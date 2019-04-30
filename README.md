# PhaseMyDeNovo
Script for phasing DNMs

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

The input file of de novos should have the following columns:

