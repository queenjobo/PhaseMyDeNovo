'''
Joanna Kaplanis 
16/10/2025


Script to phase de novo SNVs using nearby het variants

'''

#!/usr/bin/env python3

# IMPORTS ----------------------------------
import pysam
import numpy as np
import argparse
import pandas as pd
from collections import defaultdict

# mapping quality threshold
MAP_QUAL_TH = 20

# ------------ Helpers  ------------
def _is_snp(rec):
    # Single-alt SNP only 
    return (
        len(rec.alts or []) == 1
        and len(rec.ref) == 1
        and len(rec.alts[0]) == 1
    )

def _is_het(sample_call):
    gt = sample_call.get("GT")
    # GT is a tuple like (0,1) or (1,0); treat phased/unphased the same
    return gt is not None and len(gt) == 2 and gt[0] != gt[1]

def _is_hom_alt(sample_call):
    gt = sample_call.get("GT")
    return gt is not None and len(gt) == 2 and gt[0] == gt[1] and gt[0] is not None and gt[0] > 0

def _is_hom_ref(sample_call):
    gt = sample_call.get("GT")
    return gt == (0, 0)

# FUNCTIONS --------------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-dnmfile", type=str,
        help="Tab file of de novo mutations to phase. Columns: id,chrom,pos,ref,alt,vcfs,vcf_ids,cram"
    )
    parser.add_argument(
        "-id", type=str, default="",
        help="Subset to a single individual ID."
    )
    parser.add_argument(
        "-outfile", type=str, default="my_phased_denovos.tab",
        help="Output file path."
    )
    parser.add_argument(
        "--reference", type=str, default=None,
        help="(Optional) reference FASTA for CRAM if required."
    )
    return parser.parse_args()

def phase_my_dnm(vcf_ids, pos, chrom, ref, alt, vcfs, idcram, window=500, reference=None):
    """
    Phase a de novo SNV using a nearby heterozygous variant.
    Returns a numpy array of strings to append to the output row.
    """
    allinfo = np.array([])

    # pysam VariantFile.fetch uses 0-based, half-open intervals
    start = max(0, pos - 1 - window)
    end = pos + window

    with pysam.VariantFile(vcfs[0]) as vcf_child:
        for record in vcf_child.fetch(chrom, start, end):
            if not _is_snp(record):
                continue
            if record.pos == pos:
                continue  # skip the DNM itself
            child_call = record.samples.get(vcf_ids[0])
            if child_call is None or not _is_het(child_call):
                continue

            # phase genotype against parents
            gt_phase = get_gt_phase(record, vcfs, vcf_ids)
            if gt_phase == "NA":
                continue

            # read-backed phasing in child CRAM
            read_phase = get_read_phase(
                idcram, chrom, pos, ref, alt,
                record.pos, record.ref, record.alts[0],
                reference=reference
            )
            myphase = combine_phase(gt_phase, read_phase)
            if myphase != "NA":
                info = np.array([
                    str(record.pos),
                    str(record.ref),
                    str(record.alts[0]),
                    f"{read_phase[0]}|{read_phase[1]}",
                    myphase
                ])
                if allinfo.size != 0:
                    wcomma = np.core.defchararray.add(allinfo, np.full(len(allinfo), ","))
                    allinfo = np.core.defchararray.add(wcomma, info)
                else:
                    allinfo = info
    return allinfo

def combine_phase(gt_phase, read_phase):
    """
    Combine GT phase of nearby variant with read-backed phase counts.
    Returns 'F', 'M', or 'NA'.
    """
    same_phase = -1
    # strong single-sided support
    if read_phase[0] == 0 and read_phase[1] > 1:
        same_phase = 0
    elif read_phase[0] > 1 and read_phase[1] == 0:
        same_phase = 1
    # allow one read off
    elif read_phase[0] == 1 and read_phase[1] > 4:
        same_phase = 0
    elif read_phase[0] > 4 and read_phase[1] == 1:
        same_phase = 1

    if same_phase == 1:
        return gt_phase
    elif same_phase == 0:
        return "M" if gt_phase == "F" else ("F" if gt_phase == "M" else "NA")
    else:
        return "NA"

def get_variant_fromvcf(chrom, pos, vcf_path):
    """Extract a single variant by position from a VCF."""
    with pysam.VariantFile(vcf_path) as vf:
        recs = list(vf.fetch(chrom, pos - 1, pos))  # half-open
    return recs[0] if recs else None

def get_gt_phase(record_child, vcfs, vcf_ids):
    """
    Determine whether the child's het variant matches father or mother:
      - If parent is hom-alt: assign to that parent.
      - If one parent het, the other hom-ref: assign to the het parent.
    Returns 'F', 'M', or 'NA'.
    """
    phase = "NA"
    father_record = get_variant_fromvcf(record_child.chrom, record_child.pos, vcfs[1])
    mother_record = get_variant_fromvcf(record_child.chrom, record_child.pos, vcfs[2])

    if father_record is None or mother_record is None:
        return phase

    f_call = father_record.samples.get(vcf_ids[1])
    m_call = mother_record.samples.get(vcf_ids[2])
    if f_call is None or m_call is None:
        return phase

    if _is_hom_alt(f_call):
        phase = "F"
    elif _is_hom_alt(m_call):
        phase = "M"
    elif _is_het(f_call) and _is_hom_ref(m_call):
        phase = "F"
    elif _is_het(m_call) and _is_hom_ref(f_call):
        phase = "M"
    return phase

def read_pair_generator(bam, chrom, start, stop):
    """
    Generate proper read pairs within region.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(chrom, start, stop):
        if (not read.is_proper_pair
            or read.is_secondary
            or read.is_supplementary
            or read.is_duplicate):
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                mate = read_dict[qname][1]
                if mate is not None:
                    yield read, mate
            else:
                mate = read_dict[qname][0]
                if mate is not None:
                    yield mate, read
            del read_dict[qname]

def _base_at_refpos(read, position):
    """
    Get the read base aligned to a given reference position.
    Assumes 'position' is present in read.reference_positions.
    """
    refpos = read.get_reference_positions(full_length=True)
    idx = refpos.index(position)  # 0-based index; no -1
    return read.query_sequence[idx]

def get_base_combo(read1, read2, dnm_pos, var_pos):
    """
    Return haplotype combo (dnm+var) from the pair if both positions are covered.
    """
    com = ""
    r1pos = set(read1.get_reference_positions())
    r2pos = set(read2.get_reference_positions())

    if dnm_pos in r1pos:
        if var_pos in r1pos:
            com = _base_at_refpos(read1, dnm_pos) + _base_at_refpos(read1, var_pos)
        elif var_pos in r2pos:
            com = _base_at_refpos(read1, dnm_pos) + _base_at_refpos(read2, var_pos)
    elif dnm_pos in r2pos:
        if var_pos in r2pos:
            com = _base_at_refpos(read2, dnm_pos) + _base_at_refpos(read2, var_pos)
        elif var_pos in r1pos:
            com = _base_at_refpos(read2, dnm_pos) + _base_at_refpos(read1, var_pos)
    return com

def count_phases(coms, dnm_ref, dnm_alt, var_ref, var_alt):
    """
    Return [same, diff] haplotype read counts.
      same = RR + AA
      diff = AR + RA
    """
    rr = coms.count(dnm_ref + var_ref)
    ra = coms.count(dnm_ref + var_alt)
    aa = coms.count(dnm_alt + var_alt)
    ar = coms.count(dnm_alt + var_ref)
    same = rr + aa
    diff = ar + ra
    return [same, diff]

def get_read_phase(idcram, chrom, dnm_pos, dnm_ref, dnm_alt, var_pos, var_ref, var_alt, reference=None):
    """
    Read-backed phase from child CRAM around the two positions.
    """
    start = min(dnm_pos, var_pos)
    end = max(dnm_pos, var_pos)

    if reference:
        samfile = pysam.AlignmentFile(idcram, "rc", reference_filename=reference)
    else:
        samfile = pysam.AlignmentFile(idcram, "rc")

    coms = []
    try:
        for read1, read2 in read_pair_generator(samfile, chrom, start - 1, end):
            # skip if either read is low MQ
            if (read1.mapping_quality is None or read2.mapping_quality is None
                or read1.mapping_quality <= MAP_QUAL_TH
                or read2.mapping_quality <= MAP_QUAL_TH):
                continue
            com = get_base_combo(read1, read2, dnm_pos, var_pos)
            if com:
                coms.append(com)
    finally:
        samfile.close()

    return count_phases(coms, dnm_ref, dnm_alt, var_ref, var_alt)

def main():
    args = parse_args()
    dnms = pd.read_csv(args.dnmfile, sep="\t")
    with open(args.outfile, "w") as f:
        myheader = "\t".join(
            dnms.columns.tolist()
            + ["phase_var_pos", "phase_var_ref", "phase_var_alt", "AA_AR_read_support", "phase"]
        ) + "\n"
        f.write(myheader)

        if len(args.id) > 0:
            dnms = dnms[dnms.id == args.id]

        for id_ in dnms.id.unique():
            idnms = dnms[dnms.id == id_]
            for _, row in idnms.iterrows():
                # only phasing SNP DNMs
                if len(row.ref) == 1 and len(row.alt) == 1:
                    vcfs = row.vcfs.split("|")
                    vcf_ids = row.vcf_ids.split("|")
                    info = phase_my_dnm(
                        vcf_ids, int(row.pos), row.chrom, row.ref, row.alt, vcfs, row.cram,
                        reference=args.reference
                    )
                    if info.size > 0:
                        myline = "\t".join(list(map(str, list(row))) + list(info)) + "\n"
                        f.write(myline)

if __name__ == "__main__":
    main()
