Simulate structural variants in a genome. Add SNPs, small indels, large
deletions and translocations. Translocations may be conserved or non-conserved
(accompanied by a deletion of the tail end of the translocated region of no
greater than 50% size of the translocation).
Usage:
python3 SimuVar.py
    -r   input fasta filename        reference fasta file (required)
    -s   int                         number of SNPs
    -i   int                         number of indels
    -d   int                         number of large deletions
    -t   int                         number of translocations
    -nc  True/False                  add at least one (default False)
    -f   output fasta filename       name of output fasta file (required)
    -v   output vcf filename         name of output vcf file (required)

Example:
python3 SimuVar.py \
    -r Ecoli.fa \
    -s 100 \
    -i 30 \
    -d 10 \
    -t 3 \
    -nc True \
    -f ref_mut.fa \
    -v ref_mut.vcf
