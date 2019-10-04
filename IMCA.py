#!/usr/bin/python
#import sys
#import os
import copy
import argparse
import collections
import pysam
import mappy

# operations defined in CIGAR strings
# 0 M match/mismatch
# 1 I insertion
# 2 D deletion
# 3 N intron
# 4 S soft clip
# 5 H hard clip
# 6 P padding (?)
# 7 = match
# 8 X mismatch
OPERATIONS_QUERY_CONSUMING = set((0, 1, 4, 7, 8))
OPERATIONS_REFERENCE_CONSUMING = set((0, 2, 3, 7, 8))
OPERATIONS_INSERTION_LIKE = OPERATIONS_QUERY_CONSUMING - OPERATIONS_REFERENCE_CONSUMING
OPERATIONS_DELETION_LIKE = OPERATIONS_REFERENCE_CONSUMING - OPERATIONS_QUERY_CONSUMING

arguments = None
reference_for_reads = None
reference_for_contigs = None
contigs_ref = None
contigs_fa = None
reads = None
contig_mappings = collections.defaultdict(list)
FLAGS = {}
header = {}

# class mappy::Alignment:
# self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand, seg_id, cs_str, MD_str

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq',
                        action='store', type=str,
                        help='.fastq file with reads')
    parser.add_argument('-c', '--contigs',
                        action='store', type=str,
                        help='fasta file with contigs')
    parser.add_argument('-r', '--reference',
                        action='store', type=str,
                        help='fasta file with reference')
    parser.add_argument('-o', '--output',
                        action='store', type=str, default='merged.sam',
                        help='output file (defaults to merged.sam)')
    parser.add_argument('-t', '--transfer-from-reference',
                        action='store_true',
                        help='transfer reads mapped to the reference?'
                        ' By default I won\'t.')
    parser.add_argument('-k', '--keep-unmapped-contigs',
                        action='store_true',
                        help='Add unmapped contigs to the references in the output?'
                        ' By default I won\'t.')
    parser.add_argument('-a', '--ambigous-mappings',
                        action='store', type=str, default='best',
                        help='What to do with ambigous mappings?\n'
                        ' - best: keep the best mapping (default).'
                        '         if there is none, choose randomly.\n\n'
                        ' - best-rm: keep the best mapping.'
                        '         if there is none, remove this contig/read.\n'
                        '- rm: remove the ambigously mapped reads/contigs.\n'
                        '- keep: keep all the best mappings.\n')
    arguments = parser.parse_args()
    return arguments

def add_header(sam_file):
    # TODO
    pass

def calculate_flag(alignment):
    """
    From SAM Format Specification, version from 20.06.2019
    https://samtools.github.io/hts-specs/SAMv1.pdf:
    1    0x1   template having multiple segments in sequencing
    2    0x2   each segment properly aligned according to the aligner
    4    0x4   segment unmapped
    8    0x8   next segment in the template unmapped
    16   0x10  SEQ being reverse complemented
    32   0x20  SEQ of the next segment in the template being reverse complemented
    64   0x40  the first segment in the template
    128  0x80  the last segment in the template
    256  0x100 secondary alignment
    512  0x200 not passing filters, such as platform/vendor quality controls
    1024 0x400 PCR or optical duplicate
    2048 0x800 supplementary alignment
    """
    flag = 0
    if alignment.strand = -1:
        flag += 16
    if not alignment.is_primary:
        flag += 256
    return str(flag)

def alignment2sam(alignment, read_name, sequence):
    """
    From SAM Format Specification, version from 20.06.2019
    https://samtools.github.io/hts-specs/SAMv1.pdf:
    1  QNAME String [!-?A-~]{1,254}             Query template NAME
    2  FLAG  Int    [0, 2^16−1]                 bitwise FLAG
    3  RNAME String \*|[:rname:∧*=][:rname:]*   Reference sequence NAME
    4  POS   Int    [0,231−1]                   1-based leftmost mapping POSition
    5  MAPQ  Int    [0,28−1]                    MAPping Quality
    6  CIGAR String \*|([0-9]+[MIDNSHPX=])+     CIGAR string
    7  RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
    8  PNEXT Int    [0,231−1]                   Position of the mate/next read
    9  TLEN  Int    [−231+ 1,231−1]             observed Template LENgth
    10 SEQ   String \*|[A-Za-z=.]+              segment SEQuence
    11 QUAL  String [!-~]+                      ASCII of Phred-scaled base QUALity+33
    """
    # TODO: fields 7, 8, 9, potentially 11
    flag = calculate_flag(alignment)
    return '\t'.join([read_name, flag,
                      alignment.ctg, str(alignment.r_st),
                      str(alignment.mapq), alignment.cigar_str,
                      '*', '0', '0', sequence, '*'])

def write_alignment(alignment, output, read_name, sequence):
    output.write(alignment2sam(alignment, read_name, sequence))
    output.write('\n')

def proba():
    print(list(contig_mappings.items())[0])

def map_contigs(filename):
    for contig in mappy.fastx_read(filename):
        contig_name, contig_seq, contig_quality = contig
        contig_mappings[contig_name] = reference_for_contigs.map(contig_seq)

def choose_mapping(mappings):
    # TODO: it will be chosen based on FLAGS
    return mappings[0]

def find_contig_mappings(mappings2contigs):
    return [contig_mappings[mapping.ctg] for mapping in mappings2contigs]

def calculate_offset(cigar, start, stop):
    position_on_contig = start
    offset = - start
    for count, operation in cigar:
        print("Count, operation:")
        print(count, operation)
        end_of_operation = position_on_contig
        if operation in OPERATIONS_QUERY_CONSUMING:
            end_of_operation += count
        print("end of operation:")
        print(end_of_operation)
        if end_of_operation >= stop:
            print("koniec")
            break
        if operation in OPERATIONS_REFERENCE_CONSUMING:
            offset += count
        if operation in OPERATIONS_QUERY_CONSUMING:
            position_on_contig += count
            offset -= count
    print(position_on_contig, offset, stop)
    return offset
 
def transfer_mapping(read2contig, contig2ref):
    # TODO:
    # !reverse mapping
    # cigar? Do we keep the original one? Or merge with contig2ref?
    # Mapping quality - how do we want to calculate it?
    # mlen, blen attributes
    print("Read2contig:")
    print(read2contig)
    print("Contig2ref:")
    print(contig2ref)
    start_offset = calculate_offset(contig2ref.cigar, contig2ref.q_st, read2contig.r_st)
    start = read2contig.r_st + contig2ref.r_st + start_offset
    end_offset = calculate_offset(contig2ref.cigar, contig2ref.q_st, read2contig.r_en)
    #end = start + read2contig.r_en - read2contig.r_st
    end = read2contig.r_en + contig2ref.r_st + end_offset
    # ^ na pewno mozna zrobic lepiej mapowanie koncowego koordynatu
    print("start, end:")
    print(start, end)
    read2ref = mappy.Alignment(contig2ref.ctg, contig2ref.ctg_len, # chromosome
                               start, end,                         # start, end on reference
                               1,                                  # strand - todo
                               read2contig.q_st, read2contig.q_en, # start, end on query
                               read2contig.mapq,                   # mapping quality - potentially todo
                               read2contig.cigar,
                               read2contig.is_primary,
                               read2contig.mlen, read2contig.blen,
                               read2contig.NM,
                               read2contig.trans_strand,
                               read2contig.read_num - 1,
                               read2contig.cs, read2contig.MD)
    # self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand, seg_id, cs_str, MD_str
    # doesn't work for inverse mapping
    print(read2ref)
    return read2ref

def get_transferred_mapping(read2contig, contig2ref):
    # TODO: proper mapping choosing
    if contig2ref == []:
        return None
    read2contig = read2contig[0]
    contig2ref = list(contig2ref[0])[0]
    mapping = transfer_mapping(read2contig, contig2ref)
    return mapping
 
def main():
    global arguments
    global reference_for_reads
    global reference_for_contigs
    global contigs_ref
    global contigs_fa
    global reads
    global contig_mappings
    arguments = parse_arguments()
    reference_for_reads = mappy.Aligner(arguments.reference, preset='sr')
    # tu bedzie inny preset, ale teraz mam krotkie contigi
    reference_for_contigs = mappy.Aligner(arguments.reference, preset='sr')
    contigs_ref = mappy.Aligner(arguments.contigs, preset='sr')
    reads = mappy.fastx_read(arguments.fastq)
    output = open(arguments.output, 'w')
    map_contigs(arguments.contigs)
    FLAGS['transfer_from_reference'] = arguments.transfer_from_reference
    for read in reads:
        proba()
        name, seq, quality = read
        mappings2reference = list(reference_for_reads.map(seq))
        print(mappings2reference)
        if mappings2reference != [] and not FLAGS['transfer_from_reference']:
            # read mapped 2 ref and we want to keep this mapping
            mapping = choose_mapping(mappings2reference)
            write_alignment(mapping, output, name, seq)
            continue
        mappings2contigs = list(contigs_ref.map(seq))
        if mappings2reference == []:
            if mappings2contigs == []:
                # read not mapped anywhere
                # currently I don't write unmapped reads
                continue
            else:
                # read mapped only 2 contigs
                # we want to transfer mapping
                mappings = find_contig_mappings(mappings2contigs)
                mapping = choose_mapping(mappings2contigs)
                write_alignment(mapping, output. name, seq)
                continue
        else:
            if mappings2contigs == []:
                # read mapped only 2 reference
                mapping = choose_mapping(mappings2reference)
                write_alignment(mapping, output, name, seq)
                continue
            else:
                # read mapped 2 reference and 2 contigs
                # we want to transfer mapping
                # TODO: it should be in a method
                mappings = find_contig_mappings(mappings2contigs)
                mapping = get_transferred_mapping(mappings2contigs, mappings)
                if mapping is None:
                    # read mapped 2 contig, but contig not mapped to ref
                    continue
                mapping = choose_mapping(mappings2contigs)
                write_alignment(mapping, output, name, seq)
                continue
        #write_alignment(read, output)
    output.close()
    add_header(output)

if __name__ == '__main__':
    main()
