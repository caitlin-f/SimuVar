from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as ss
import random
import datetime
import os

"""
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
"""


class Genome:
    """ Representation of the simulated genome """

    def __init__(self, refseq):
        self.name = refseq.id
        self.seq = refseq.seq
        self.mut_seq = refseq.seq.tomutable()
        self.unavail_pos = [] # breakpoint or deleted positions
        self.__variants = [] # list of Variants in descending order of start pos

    def add_variant(self, variant):
        """ Add variant to list of variants in descending order of pos """
        self.__variants.append(variant)
        self.__variants.sort(reverse=True)

    def get_variants(self):
        return self.__variants

    def __str__(self):
        var_string = ""
        for var in reversed(self.__variants):
            var_string += str(var)
            var_string += '\n'
        var_string = var_string.strip('\n')
        return("{}\n{}".format(self.name, var_string))

    def __repr__(self):
        var_string = ""
        for var in self.__variants:
            var_string += str(var)
            var_string += '\n'
        var_string = var_string.strip('\n')
        return("{}\n{}".format(self.name, var_string))


class Variant:
    """ Variants to be introduced.
    Translocation origin end is the insert position.
    Translocation insert end is the origin start position """

    def __init__(self, type, start, end, size):
        self.type = type
        self.start = start
        self.end = end
        self.size = size
        self.ref = None
        self.alt = None

    def __str__(self):
        return("{} start:{}, end:{}, size:{}, ref:{}, alt:{}".format(self.type, self.start, self.end, self.size, self.ref, self.alt))

    def __repr__(self):
        return("{} start:{}, end:{}, size:{}, ref:{}, alt:{}".format(self.type, self.start, self.end, self.size, self.ref, self.alt))

    def __eq__(self, other):
        return self.start == other.start

    def __ne__(self, other):
        return self.start == other.start

    def __lt__(self, other):
        return self.start < other.start

    def __le__(self, other):
        return self.start <= other.start

    def __gt__(self, other):
        return self.start > other.start

    def __ge__(self, other):
        return self.start >= other.start

    def __hash__(self):
        return int("{}{}{}".format(hash(self.type), self.start, self.end))


def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-r', type=str, required=True, metavar='[reference.fa]', help="Reference file in fasta format")
    parser.add_argument('-s', type=int, metavar='number of SNPs', help="Number of SNPs to introduce")
    parser.add_argument('-i', type=int, metavar='number of indels', help="Number of small indels to introduce. Indel sizes randomly sampled between 1-100bp with a ratio of 9:1 for size being < 11bp")
    parser.add_argument('-d', type=int, metavar='number of large deletions', help="Number of large deletions to introduce. Deletion size randomly sampled between 100-5000bp")
    parser.add_argument('-t', type=int, metavar='number of translocations', help="number of translocations, regions will be randomly selected ranging between 500-5000bp in size")
    parser.add_argument('-nc', type=bool, default=False, metavar='add at least one non-conserved translocation', help="True/False add one non-conserved translocation. Deletion occurs from end of the translocation. Default=False")
    parser.add_argument('-f', type=str, required=True, metavar='output fasta filename', help="Output fasta filename and file path")
    parser.add_argument('-v', type=str, required=True, metavar='output vcf file of mutations', help="Output vcf filename and file path")

    return parser.parse_args()


# SNPs
def define_snps(genome, num):
    """ Define location of SNP """
    for n in range(num):
        snp_pos = get_snp_pos(genome)
        var = Variant("snp", snp_pos, snp_pos, 0)
        genome.add_variant(var)
        genome.unavail_pos.append(snp_pos)

def get_snp_pos(genome):
    """ Helper function. SNP position generator """
    snp_pos = random.randint(100, len(genome.seq)-100)
    if snp_pos in genome.unavail_pos:
        snp_pos = get_snp_pos(genome)
    return snp_pos


# small indels
def define_indels(genome, num):
    """ Define location, size and type of indel (deletion or insertion) between
    1bp-100bp in length """
    # create list of indel sizes to sample from
    # 1-100bp indels with ratio 9:1 of indels being < 11bp
    size_range = indel_sizes()
    for n in range(num):
        start_pos, end_pos, size = get_indel_pos(genome, size_range)
        var = Variant("indel", start_pos, end_pos, size)
        genome.add_variant(var)
        # add to unavail list
        for j in range(start_pos, end_pos):
            genome.unavail_pos.append(j)

def get_indel_pos(genome, size_range):
    """ Helper function. Indel type and position generator """
    start_pos = random.randint(100,len(genome.seq)-110) # positions 100bp from start or end will not be variable
    size = size_range[random.randint(0, len(size_range)-1)]
    end_pos = start_pos + random.randint(1,10)

    if random.randint(0,1) == 0: # insertion
        end_pos = start_pos+1
    else:
        end_pos = start_pos + size
        size *= -1

    unavail = False
    for n in range(start_pos, end_pos):
        if n in genome.unavail_pos:
            unavail = True
            break
    if unavail:
        start_pos, end_pos = get_del_pos(genome)
    return (start_pos, end_pos, size)

def indel_sizes():
    """ Build list of possible indel sizes to select from, with a 9:1 ratio of
    indels < 11bp in length """
    size_range = []
    for n in range(90):
        size_range.extend([1,2,3,4,5,6,7,8,9,10])
    for n in range(11,101):
        size_range.append(n)
    return size_range


# large deletions
def define_deletions(genome, num):
    """ Define location and size of large deletions 100bp-5000bp in length """
    start = []
    end = []
    for n in range(num):
        start_pos, end_pos = get_del_pos(genome)
        # add deletion Variants to genome list
        var = Variant("deletion", start_pos, end_pos, start_pos-end_pos)
        genome.add_variant(var)
        # add to unavail list
        for j in range(start_pos, end_pos):
            genome.unavail_pos.append(j)

def get_del_pos(genome):
    """ Helper function. Deletion position generator """
    start_pos = random.randint(100,len(genome.seq)-5100) # positions 100bp from start or end will not be variable
    end_pos = start_pos + random.randint(100,5000)
    unavail = False
    for n in range(start_pos, end_pos):
        if n in genome.unavail_pos:
            unavail = True
            break
    if unavail:
        start_pos, end_pos = get_del_pos(genome)
    return (start_pos, end_pos)


# translocations
def define_translocations(genome, num, nc):
    """
    Define the translocation variants plus deletions in the case on
    non-conservative translocations
    params:
        genome: Genome
        pos_list: translocation positions
        nc: True/False non-conserved translocation
    """
    start = []
    end = []
    for n in range(num):
        start_pos = random.randint(100,len(genome.seq)-5100) # positions 100bp from start or end will not be variable
        end_pos = start_pos + random.randint(500,5000)
        start.append(start_pos)
        end.append(end_pos)

    if nc: # if non-conservative translocations specified
        del_start = []
        del_end = []
        nc_pos = [p for p in range(0, len(start))]
        for n in range(0,len(start),2): # 50:50 chance that half will be non-conserved
            if not del_start or random.randint(0,1) == 0: # ensures at least 1 will be non-conserved
                length = len(nc_pos)
                pop_pos = random.randint(0,length-1)
                idx = nc_pos.pop(pop_pos)
                nc_size = random.randint(100, ((end[idx]-start[idx])//2)-1) # size between 100 and half the translocation size
                start_pos = end[idx]-nc_size
                end_pos = end[idx]
                del_start.append(start_pos)
                del_end.append(end_pos)
                end[idx] = start_pos
                # add new deletion Variants to genome list
                var = Variant("deletion", start_pos, end_pos, start_pos-end_pos)
                genome.add_variant(var)
                # add new deletions to unavail list
                for j in range(start_pos, end_pos):
                    genome.unavail_pos.append(j)

    # add translocation Variants to genome list
    for v in range(len(start)):
        pos = get_trans_pos(genome) # get new position
        # add either side of insertion point to unavail list
        genome.unavail_pos.append(pos-1)
        genome.unavail_pos.append(pos)
        genome.unavail_pos.append(pos+1)
        # add translocated region to unavail list
        for j in range(start[v], end[v]):
            genome.unavail_pos.append(j)
        # add Variant to genome's variant list
        var = Variant("translocation origin", start[v], pos, end[v]-start[v])
        genome.add_variant(var)
        var = Variant("translocation insert", pos, start[v], end[v]-start[v])
        genome.add_variant(var)

def get_trans_pos(genome):
    """ Helper function. Translocation position generator """
    pos = random.randint(100, len(genome.seq)-100) # insert position
    if pos in genome.unavail_pos:
        pos = get_trans_pos(genome)
    return pos


def mutate_seq(genome):
    """ Mutates the reference sequence with variants """
    for var in genome.get_variants():
        if var.type == "snp":
            mutate_snp(genome, var)
        elif var.type == "indel":
            mutate_indel(genome, var)
        elif var.type == "deletion":
            mutate_deletion(genome, var)
        elif var.type == "translocation origin":
            mutate_trans_orig(genome, var)
        elif var.type == "translocation insert":
            mutate_trans_ins(genome, var)


def mutate_snp(genome, var):
    """ Mutate reference with snp """
    nt_options = {'A':['T','G','C'], 'T':['A','G','C'], 'G':['A','T','C'], 'C':['A','T','G']}
    n = random.randint(0,2)
    nt = nt_options.get(genome.seq[var.start])[n]
    genome.mut_seq[var.start] = nt

    var.ref = genome.seq[var.start]
    var.alt = nt


def mutate_indel(genome, var):
    """ Mutation reference with small indel """
    if var.size > 0: # insertion
        if var.size <= 10:
            new_seq = small_insert(var)
        if var.size > 10:
            new_seq = large_insert(var)
        for nt in new_seq:
            genome.mut_seq.insert(var.start,nt)

        var.ref = "."
        var.alt = "".join(new_seq)

    else: # deletion
        for i in range(abs(var.size)):
            genome.mut_seq.pop(var.start)
        var.ref = genome.seq[var.start:var.end]
        var.alt = "."

def small_insert(var):
    """ Helper function. Homopolymer string builder for small indel insertion """
    nt_options = ['A', 'T', 'G', 'C']
    n = random.randint(0,3)
    insert = []
    for i in range(var.size):
        insert.append(nt_options[n])
    return insert

def large_insert(var):
    """ Helper function. Random sequence generator for large indels insertion"""
    nt_options = ['A', 'T', 'G', 'C']
    insert = []
    for i in range(var.size):
        n = random.randint(0,3)
        insert.append(nt_options[n])
    return insert


def mutate_deletion(genome, var):
    for i in range(abs(var.size)):
        genome.mut_seq.pop(var.start)
    var.ref = genome.seq[var.start:var.end]
    var.alt = "."


def mutate_trans_orig(genome, var):
    for i in range(var.size):
        genome.mut_seq.pop(var.start)
    var.ref = genome.seq[var.start:var.start+var.size]
    var.alt = "."


def mutate_trans_ins(genome, var):
    insert = list(genome.seq[var.end:var.end+var.size])
    for nt in insert:
        genome.mut_seq.insert(var.start,nt)
    var.ref = "."
    var.alt = "".join(insert)


def write_vcf(file, ref, genome):
    """ Write vcf of introduced mutations """
    vcf = open(file, "w")
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("##fileDate={}\n".format(datetime.datetime.today().strftime('%Y%m%d')))
    vcf.write("##source={}\n".format(os.path.basename(__file__)))
    vcf.write("##reference={}\n".format(ref))
    vcf.write("##contig=<ID={}\n".format(genome.name))
    vcf.write("##translocation_origin=END position relates to insertion position\n")
    vcf.write("##translocation_insert=END position relates to original start position\n")
    vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
    vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant detected">\n')
    vcf.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of variant region">\n')
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    for var in reversed(genome.get_variants()):
        info = "END={};SVTYPE={};LEN={}".format(var.end,var.type,var.size)
        vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        genome.name,var.start,'.',var.ref,var.alt,'.','.',info,'.'))

    vcf.close()


def main():
    args = parse_args()

    # read in reference
    for record in SeqIO.parse(args.r, "fasta"):
        ref_seq = record

    genome = Genome(ref_seq)

    # generate non-overlapping structural variants
    print("Defining variants")
    if args.s:
        define_snps(genome, args.s)
    if args.i:
        define_indels(genome, args.i)
    if args.t:
        define_translocations(genome, args.t, args.nc)
    if args.d:
        define_deletions(genome, args.d)

    # mutate reference with variants
    print("Mutating sequence")
    mutate_seq(genome)

    print("Writing output")
    new_seq = SeqRecord(genome.mut_seq)
    new_seq.id = genome.name
    new_seq.description = "| SimuVar mutated sequence"
    SeqIO.write(new_seq, args.f, "fasta")

    write_vcf(args.v, args.r, genome)


if __name__ == '__main__':
    main()
