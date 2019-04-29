'''
Script to phase de novo mutations for 10k

Joanna Kaplanis 
24/04/2019

DNMFILE = /nfs/ddd0/ks20/active_metaDNM_analyses/DDD_RUMC_GDX_denovos_filtered_31058_ntrios_2018_12_21.tab

bsub -q long -R'select[mem>20000] rusage[mem=20000]' -M20000 -o phase.out python /nfs/users/nfs_j/jk18/PhD/ddd_denovo_enrichment/phasing/phase_denovos_v2.py -dnmfile /nfs/ddd0/ks20/active_metaDNM_analyses/DDD_RUMC_GDX_denovos_filtered_31058_ntrios_2018_12_21.tab

'''

# IMPORTS ----------------------------------

import pysam
import vcf
import numpy as np
import argparse
import pandas as pd 
from collections import defaultdict

TRIO_INFO = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2017-12-15/trios_bams_vcfs.txt"
BASE_QUAL_TH = 11
MAP_QUAL_TH = 20

# FUNCTIONS --------------------------------

def parse_args():
    ''' parse arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-dnmfile',type = str, help = "file with de novos")
    parser.add_argument('-id',type = str, default = "", help = "id to subset to")
    parser.add_argument('-outfile',type = str, default = "phased_output.tab",help = "output file path")
    args = parser.parse_args()
    return(args)

def get_files(id):
    '''get ddd files '''
    idvcf = ""
    idcram = ""
    mainids = ""
    with open(TRIO_INFO,'r') as f:
        for line in f:
            fields = line.strip("\n").split("\t")
            if id == fields[2]:
                idvcf = fields[14]
                idcram = "/lustre/scratch115/projects/ddd/EGA_submission/" + fields[5] + ".cram"
                mainids = fields[5:8]
                break
    return(idvcf,idcram,mainids)

def phase_my_dnm(mainids,pos,chrom,ref,alt,idvcf,idcram,window = 500):
    ''' phase my de novo! '''
    allinfo = np.array([])
    vcf_reader = vcf.Reader(filename=idvcf,compressed=True)
    records = vcf_reader.fetch(chrom,pos-window,pos+window)
    for record in records:
        #is there a phasable het variant within 500bp?
        if record.POS != pos and record.genotype(mainids[0]).is_het and not(record.genotype(mainids[1]).is_het and record.genotype(mainids[2]).is_het):
            gt_phase = get_gt_phase(record,mainids)
            if gt_phase != "NA":
                read_phase = get_read_phase(idcram,chrom,pos,ref,alt, record.POS, str(record.REF), str(record.ALT[0]))
                myphase = combine_phase(gt_phase,read_phase)
                if myphase != "NA":
                    #create list of info to return
                    info = np.array([str(record.POS), str(record.REF), str(record.ALT[0]), str(read_phase[0])+"|"+ str(read_phase[1]),myphase])
                    if allinfo.size != 0:
                        np.core.defchararray.add(allinfo,info)
                    else:
                        allinfo = info  
    return(allinfo)

def combine_phase(gt_phase,read_phase):
    '''get ultimate phase'''
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
        

def get_gt_phase(record,mainids):
    ''' phase genotype''' 
    phase = "NA"
    if record.genotype(mainids[1]).gt_type == 2:
        phase = "F"
    elif record.genotype(mainids[2]).gt_type == 2:
        phase = "M" 
    elif record.genotype(mainids[1]).gt_type == 1 and record.genotype(mainids[2]).gt_type == 0:
        phase = "F"
    elif record.genotype(mainids[2]).gt_type == 1 and record.genotype(mainids[1]).gt_type == 0:
        phase = "M"
    return(phase)

def read_pair_generator(bam,chrom,start,stop):
    '''
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    '''
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(chrom,start,stop):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
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
    bq = read.query_qualities[position - read.reference_start -1]
    return(bq)


def get_base(read, position):
    '''given pysam read and position extract base'''
    base = read.query_sequence[position - read.reference_start -1]
    return(base)

def get_base_combo(read1,read2,dnm_pos,var_pos):
    ''' given reads and positions extract base phase on read '''
    com = ""
    #get positions in read
    read1pos = read1.get_reference_positions()
    read2pos = read2.get_reference_positions()
    #check if dnm and variant are in same read pair and get corresponding bases
    if dnm_pos in read1pos:
        if var_pos in read1pos:
            if get_base_quality(read1,dnm_pos) > BASE_QUAL_TH and get_base_quality(read1,var_pos) > BASE_QUAL_TH:
                com = get_base(read1,dnm_pos) + get_base(read1,var_pos)
        elif var_pos in read2pos:
            if get_base_quality(read1,dnm_pos) > BASE_QUAL_TH and get_base_quality(read2,var_pos) > BASE_QUAL_TH:
                com = get_base(read1,dnm_pos) + get_base(read2,var_pos)
    elif dnm_pos in read2pos:
        if var_pos in read2pos:
            if get_base_quality(read2,dnm_pos) > BASE_QUAL_TH and get_base_quality(read2,var_pos) > BASE_QUAL_TH:
                com = get_base(read2,dnm_pos) + get_base(read2,var_pos)
        elif var_pos in read1pos:
            if get_base_quality(read2,dnm_pos) > BASE_QUAL_TH and get_base_quality(read1,var_pos) > BASE_QUAL_TH:
                com = get_base(read2,dnm_pos) + get_base(read1,var_pos)
    return(com)

def count_phases(coms, dnm_ref,dnm_alt,var_ref,var_alt):
    '''
    Count haplotyples of counts
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
    Get read phase given the cram file, chrom, de novo position and phasable variant position
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
    counts = count_phases(coms, dnm_ref,dnm_alt,var_ref,var_alt)
    return(counts)

def main():
    args = parse_args()
    dnms  = pd.read_csv(args.dnmfile,sep = "\t")
    with open(args.outfile,'w') as f:
        myheader = "\t".join(dnms.columns.tolist()+ ["phase_var_pos","phase_var_ref","phase_var_alt","AA|AR_read_support","phase"])+"\n"
        f.write(myheader)
        #print header here
        if len(args.id) > 0:
            dnms = dnms[dnms.id == args.id]
        for id in dnms.id.unique():
            idnms = dnms[dnms.id == id]
            idvcf, idcram, mainids = get_files(id)
            if idvcf == "":
                continue
            for index,row in idnms.iterrows():
                info = phase_my_dnm(mainids,int(row.pos),row.chrom,row.ref,row.alt,idvcf,idcram)
                if len(info) > 0:
                    myline = "\t".join(list(map(str,list(row))) + list(info)) + "\n"
                    f.write(myline)

        

# SCRIPT -----------------------------------

if __name__=='__main__':
    main()