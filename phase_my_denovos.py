'''
Joanna Kaplanis 
24/04/2019


Script to phase de novo SNVs using nearby het variants within 500bp 

'''

# IMPORTS ----------------------------------

import pysam
import vcf
import numpy as np
import argparse
import pandas as pd 
from collections import defaultdict

#mapping quality threshold
MAP_QUAL_TH = 20

# FUNCTIONS --------------------------------

def parse_args():
    ''' parse input arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-dnmfile',type = str, help = "file with de novos")
    parser.add_argument('-id',type = str, default = "", help = "id to subset to")
    parser.add_argument('-outfile',type = str, default = "phased_output.tab",help = "output file path")
    args = parser.parse_args()
    return(args)

def phase_my_dnm(vcf_ids,pos,chrom,ref,alt,vcfs,idcram,window = 500):
    '''
    Phase my de novo!
    Returns numpy array with information about phase information for de novo SNV. 
    '''
    allinfo = np.array([])
    vcf_reader = vcf.Reader(filename=vcfs[0],compressed=True)
    records = vcf_reader.fetch(chrom,pos-window,pos+window)
    for record in records:
        print(record)
        #is there a het variant within 500bp?
        if record.is_snp and record.POS != pos and record.genotype(vcf_ids[0]).is_het:
            #phase genotype
            gt_phase = get_gt_phase(record,vcfs,vcf_ids)
            # if genotype is phased then read phase in child CRAM file
            if gt_phase != "NA":
                read_phase = get_read_phase(idcram,chrom,pos,ref,alt, record.POS, str(record.REF), str(record.ALT[0]))
                #combine phase information from GT and read to get ultimate phase
                myphase = combine_phase(gt_phase,read_phase)
                #if phased then append info
                if myphase != "NA":
                    #create list of info to return
                    info = np.array([str(record.POS), str(record.REF), str(record.ALT[0]), str(read_phase[0])+"|"+ str(read_phase[1]),myphase])
                    if allinfo.size != 0:
                        wcomma = np.core.defchararray.add(allinfo,np.full(len(allinfo),',')
                        allinfo = np.core.defchararray.add(wcomma,info)
                    else:
                        allinfo = info  
    return(allinfo)

def combine_phase(gt_phase,read_phase):
    '''
    Combine information from the genotype phase of the nearby variant and the haplotype 
    counts from the read phase to get the phase of the de novo SNV.
    Returns F (father), M (mother) or NA (not able to phase)
    '''
    same_phase = -1
    #if only one combination has >0 reads
    if read_phase[0] == 0 and read_phase[1] >0:
        same_phase = 0
    elif read_phase[0] > 0 and read_phase[1] == 0:
        same_phase = 1
    #allow for one read off
    elif read_phase[0] == 1 and read_phase[1] > 4:
        same_phase = 0
    elif read_phase[0] > 4 and read_phase[1] == 1:
        same_phase = 1
    #assign to correct parent
    if same_phase == 1:
        myphase = gt_phase
    elif same_phase == 0:
        if gt_phase == "F":
            myphase = "M"
        elif gt_phase == "M":
            myphase = "F"
    else:
        myphase = "NA"
    return(myphase)
        
def get_variant_fromvcf(chrom,pos,vcf_path,sample):
    ''' Extract variant from VCF'''
    vcf_reader = vcf.Reader(filename=vcf_path,compressed=True)
    records = vcf_reader.fetch(chrom,pos-1,pos)
    variants = [record for record in records]
    if len(variants) == 0:
        record = None
    else:
        record = variants[0]
    return(record)

def get_gt_phase(record,vcfs,vcf_ids):
    '''
    Phase het variant in child when give the variant record and paths to father/mother VCFs
    Returns F (father), M (mother) or NA (not able to phase)
    ''' 
    phase = "NA"
    father_record = get_variant_fromvcf(record.CHROM,record.POS,vcfs[1],vcf_ids[1])
    mother_record = get_variant_fromvcf(record.CHROM,record.POS,vcfs[2],vcf_ids[2])
    if father_record is not None and mother_record is not None:
        if father_record.genotype(vcf_ids[1]).gt_type == 2:
            phase = "F"
        elif mother_record.genotype(vcf_ids[2]).gt_type == 2:
            phase = "M" 
        elif father_record.genotype(vcf_ids[1]).gt_type == 1 and mother_record.genotype(vcf_ids[2]).gt_type == 0:
            phase = "F"
        elif mother_record.genotype(vcf_ids[2]).gt_type == 1 and father_record.genotype(vcf_ids[1]).gt_type == 0:
            phase = "M"
    return(phase)

def read_pair_generator(bam,chrom,start,stop):
    '''
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    '''
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(chrom,start,stop):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def get_base_quality(read, position):
    '''given pysam read and position extract base quality '''
    fullpos = read.get_reference_positions(full_length = True)
    bq = read.query_qualities[fullpos.index(position) -1]
    return(bq)


def get_base(read, position):
    '''given pysam read and position extract base'''
    fullpos = read.get_reference_positions(full_length = True)
    base = read.query_sequence[fullpos.index(position) -1]
    return(base)

def get_base_combo(read1,read2,dnm_pos,var_pos):
    ''' 
    Given paired reads and positions extract base phase on reads.
    Returns a string of the haplotype combination in the order de novo mutation and nearby variant.
    '''
    com = ""
    #get positions in read
    read1pos = read1.get_reference_positions()
    read2pos = read2.get_reference_positions()
    #check if dnm and variant are in same read pair and get corresponding bases
    if dnm_pos in read1pos:
        if var_pos in read1pos:
                com = get_base(read1,dnm_pos) + get_base(read1,var_pos)
        elif var_pos in read2pos:
                com = get_base(read1,dnm_pos) + get_base(read2,var_pos)
    elif dnm_pos in read2pos:
        if var_pos in read2pos:
                com = get_base(read2,dnm_pos) + get_base(read2,var_pos)
        elif var_pos in read1pos:
                com = get_base(read2,dnm_pos) + get_base(read1,var_pos)
    return(com)

def count_phases(coms, dnm_ref,dnm_alt,var_ref,var_alt):
    '''
    Count haplotyples of counts.
    Returns a list of two numbers the first is the number of reads indicating the de novo 
    lies on the same haplotype as the variant. The second is the number of reads indicating
    the de novo lies on the opposite haplotype to the variant.
    '''
    rr = coms.count(dnm_ref+var_ref)
    ra = coms.count(dnm_ref+var_alt)
    aa = coms.count(dnm_alt+var_alt)
    ar = coms.count(dnm_alt+var_ref)
    same = rr +aa
    diff = ar + ar
    return([same,diff])
           
def get_read_phase(idcram,chrom,dnm_pos,dnm_ref,dnm_alt,var_pos,var_ref,var_alt):
    '''
    Get read phase given the CRMAM file for the child, chrom, de novo position and het variant position
    '''
    #create start and end to pull from cram file
    if dnm_pos > var_pos:
        start = var_pos
        end = dnm_pos
    else:
        start = dnm_pos
        end = var_pos
    #read in cram file
    samfile = pysam.AlignmentFile(idcram,"rc")
    #iterate through read pairs
    coms = []
    for read1, read2 in read_pair_generator(samfile,chrom,start,end):
        #filter reads on mapping quality
        if read1.mapping_quality <= MAP_QUAL_TH and read2.mapping_quality <= MAP_QUAL_TH:
            continue
        com = get_base_combo(read1,read2,dnm_pos,var_pos)
        #print combo
        if len(com) > 0:
            coms = coms + [com]
    samfile.close()
    #count haplotype combinations
    counts = count_phases(coms, dnm_ref,dnm_alt,var_ref,var_alt)
    return(counts)

def main():
    args = parse_args()
    dnms  = pd.read_csv(args.dnmfile,sep = "\t")
    with open(args.outfile,'w') as f:
        myheader = "\t".join(dnms.columns.tolist()+ ["phase_var_pos","phase_var_ref","phase_var_alt","AA_AR_read_support","phase"])+"\n"
        f.write(myheader)
        if len(args.id) > 0:
            dnms = dnms[dnms.id == args.id]
        for id in dnms.id.unique():
            idnms = dnms[dnms.id == id]
            for index,row in idnms.iterrows():
                #only phasing SNPs
                if len(row.ref) == 1 and len(row.alt) == 1:
                    #parse vcf input
                    vcfs = row.vcfs.split("|")
                    vcf_ids = row.vcf_ids.split("|")
                    #phase
                    info = phase_my_dnm(vcf_ids,int(row.pos),row.chrom,row.ref,row.alt,vcfs,row.cram)
                    #write phased dnms to file
                    if len(info) > 0:
                        myline = "\t".join(list(map(str,list(row))) + list(info)) + "\n"
                        f.write(myline)

        

# SCRIPT -----------------------------------

if __name__=='__main__':
    main()