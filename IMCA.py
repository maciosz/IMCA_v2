#!/usr/bin/python
#import sys
#import os
import argparse
import collections
import pysam
import mappy

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

def alignment2sam(alignment, read_name, sequence):
    flag = 0
    return '\t'.join([read_name, flag,
                      aligment.ctg, alignment.r_st,
                      alignment.mapq, alignment.cigar_str,
                      '*', '0', '0', sequence, '*'])

def main():
    arguments = parse_arguments()
    reference = mappy.Aligner(arguments.reference)
    contigs_ref = mappy.Aligner(arguments.contigs)
    contigs_fa = mappy.fastx_read(arguments.contigs)
    reads = mappy.fastx_read(arguments.fastq)
    contig_mappings = collections.defaultdict(list)
    for contig in contigs_fa:
        contig_name, contig_seq, contig_quality = contig
        contig_mappings[contig_name] = reference.map(contig_seq)
    for read in reads:
        name, seq, quality = read
        mappings2reference = list(reference.map(seq))
        mappings2contigs = list(contigs_ref.map(seq))
        if mappings2reference == []:
            if mappings2contigs == []:
                # read not mapped anywhere
                pass
            else:
                # read mapped only 2 contigs
        else:
            if mappings2contigs == []:
                # read mapped only 2 reference
                pass
            else:
                # read mapped 2 reference and 2 contigs
                if arguments.transfer_from_reference:
                    # we want to transfer mapping
                    pass
                else:
                    # we ignore contig mappings
                    pass
        for mapping in mapping2contigs:
            contig_mapping = contig_mappings[mapping.ctg]
            if contig_mapping == []:
                # read mapped 2 contig, contig not mapped 2 reference
                pass
            else:
                # read mapped 2 contig, contig mapped 2 reference
                pass

if __name__ == '__main__':
    main()
