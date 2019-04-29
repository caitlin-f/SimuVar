Simulate structural variants in a genome. Add SNPs, small indels, large
deletions and translocations. Translocations may be conserved or non-conserved
(accompanied by a deletion of the tail end of the translocated region of no
greater than 50% size of the translocation).  

### Usage:
```
python3 SimuVar.py -r ref.fa -s 100 -i 30 -d 10 -t 3 -nc True -f ref_mut.fa -v ref_mut.vcf

Arguments:
	-h, --help 		show this help message and exit  
    -r   STR        reference fasta file (required)
    -s   INT		number of SNPs
    -i   INT		number of indels
    -d   INT		number of large deletions
    -t   INT		number of translocations
    -nc  BOOL		add at least one non-conserved translocation True/False (default False)
    -f   STR		name of output fasta file (required)
    -v   STR		name of output vcf file (required)

```