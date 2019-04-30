# PhaseMyDeNovo
Script for phasing DNMs. The method is generally to look for a nearby heterozygous variant within 1MB of the de novo SNV. Check the VCFs to see if the variant was inhertied from the mother or the father. Then check the reads in the child CRAM file to see if the de novo variant lies on the same or opposite haplotype as this variant. 
To determine phasing I currently use the following rules:


```
usage: phase_my_denovos.py [-h] [-dnmfile DNMFILE] [-id ID] [-outfile OUTFILE]

optional arguments:
  -h, --help        show this help message and exit
  -dnmfile DNMFILE  Tab file with de novo mutations to phase. Must contain the
                    following columns: id,chrom,pos,ref,alt,vcfs,vcf_ids,cram
  -id ID            String of ID to subset file to if only want to phase
                    variants from a specific individual.
  -outfile OUTFILE  Output file path.
```

### INPUT

The input file of de novos should have the following columns:

id - individual id  
chrom - chromosome  
pos - genomic position  
ref - reference allele   
alt - alternatie allele  
vcfs - path to vcfs in the order of child vcf, father vcf, mother vcf separated by "|". Example: 'child.vcf.gz|father.vcf.gz|mother.vcf.gz'  
vcf_ids - list of sample IDs as they appear in the VCFs in the order of child sample id, father sample id, mother sample id separated by "|". Example 'sample2494|sample5792|sample9897'  
cram - path to CRAM file for the child with the de novo mutation  

### OUTPUT
Note: only variants that were able to be phased are outputted.

The output file of phased de novos will have all the original columns from the input file as well as:

phase_var_pos - genomic position of the nearby variant  
phase_var_ref - the reference allele of the nearby variant  
phase_var_alt - the alternative allele of the nearby variant  
AA_AR_read_support - evidence of read support for phasing: number of reads supporting that the variant is on the same haplotype as the de novo SNV and number of reads supporting that the variant is on the opposite haplotupe as the de novo SNV separated by a "|". Example 0|12  
phase - what parent the de novo variant phases to either "F" for father or "M" for mother



