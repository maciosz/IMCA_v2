#!/usr/bin/python
#import sys
#import os
import copy
import argparse
import collections
import pysam
import mappy
from mapped_segment import MappedSegment

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
reads2reference = None
reads2contigs = None
contigs2reference = None
reference_for_reads = None
reference_for_contigs = None
contigs_as_reference = None
all_reads = set()
output = None
# class mappy::Alignment:
# self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand, seg_id, cs_str, MD_str

def parse_arguments():
    parser = argparse.ArgumentParser(description=
            'You should provide reads and contigs.'
            ' If they are already mapped and in sam/bam format,'
            ' you should also provide reads2contigs.'
            ' If they are in fasta/fastq format and you want me to map them,'
            ' you should also provide reference.')
    parser.add_argument('-r', '--reads',
                        action='store', type=str,
                        help='reads; either bam/sam file with reads mapped to reference or fasta/fastq file to be mapped to reference with minimap.')
    parser.add_argument('-c', '--contigs',
                        action='store', type=str,
                        help='contigs; either bam/sam file with contigs mapped to reference or fasta/fastq file to be mapped to reference with minimap.')
    parser.add_argument('--reference',
                        action='store', type=str,
                        help='fasta file with reference. ')
    parser.add_argument('--reads2contigs',
                        action='store', type=str,
                        help='sam/bam file with reads mapped to contigs.')
    parser.add_argument('-o', '--output',
                        action='store', type=str, default='merged.sam',
                        help='output file (defaults to merged.sam)')
    #parser.add_argument('-t', '--transfer-from-reference',
    #                    action='store_true',
    #                    help='transfer reads mapped to the reference?'
    #                    ' By default I won\'t.')
    parser.add_argument('--mapped2ref',
                        action='store', type=str,
                        default='condition',
                        help='Read is mapped to reference. What should I do with it?'
                             ' 1. leave - leave it alone, don\'t touch it.'
                             ' 2. transfer - if it can be transferred via contigs, do it.'
                             ' 3. check - check if it should be tranferred.'
                             ' If read meets some condition, it will be transferred,'
                             ' otherwise it will be left alone.'
                             ' You can specify the condition with "--condition" argument.')
    parser.add_argument('--condition',
                        action='store', type=str, nargs='+',
                        help='So, you specified "check" in "mapped2ref" argument.'
                             ' What condition should read meet to be transferred?'
                             ' 1. ambiguous - have more than one mapping.'
                             ' 2. quality [x] - have mapping quality below x.'
                             ' 3. cigar? Sth depending on CIGAR string. Obviously not implemented yet.'
                             '! You can specify more than one option here!')
    parser.add_argument('-k', '--keep-unmapped-contigs',
                        action='store_true',
                        help='Add unmapped contigs to the references in the output?'
                        ' By default I won\'t.')
    parser.add_argument('--test',
                        action='store', type=str)
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
    check_arguments(arguments)
    return arguments

def check_arguments(args):
    # TODO
    pass

def proba():
    print(list(contig_mappings.items())[0])

def guess_fileformat(filename):
    if filename.endswith("bam") or filename.endswith("sam"):
        return "bam"
    elif filename.endswith("fasta") or filename.endswith("fa"):
        return "fasta"
    else:
        raise ValueError("I can't guess format of %s."
                " Currently I support sam, bam and fasta files"
                " and I expect .sam, .bam, .fa or .fasta extension."
                " I can't guess without the proper name, sorry."
                % filename)

def map_with_mappy(query, subject):
    """
    Here I could already choose best mappings or sth like that.
    """
    result = collections.defaultdict(list)
    for read in query:
        name, seq, quality = read
        mappings = subject.map(seq)
        for mapping in mappings:
            mapping = MappedSegment(mapping, "mappy")
            mapping.name = name
            result[name].append(mapping)
    return result

def map_reads2reference(fasta):
    #return [reference.map(read) for read in fasta]
    read_in_reference_for_reads()
    if reference_for_reads is None:
        raise ValueError("Reference not read in."
                         " I can't map to it if it's not read in. Duh.")
    return map_with_mappy(fasta, reference_for_reads)

def map_contigs2reference(fasta):
    #return [reference.map(read) for read in fasta]
    read_in_reference_for_contigs()
    if reference_for_contigs is None:
        raise ValueError("Reference not read in."
                         " I can't map to it if it's not read in. Duh.")
    return map_with_mappy(fasta, reference_for_contigs)

def map2contigs(fasta):
    #return [contigs_as_reference.map(read) for read in fasta]
    if contigs_as_reference is None:
        raise ValueError("Contigs not read in as reference."
                         " I can't map to them if they're not read in as reference.")
    return map_with_mappy(fasta, contigs_as_reference)

def read_in_bam(filename):
    bam = pysam.AlignmentFile(filename)
    #bam = change_bam2list(bam)
    bam = change_bam2dict(bam)
    return bam

def read_in_fasta_as_query(filename):
    return mappy.fastx_read(filename)

def read_in_fasta_as_reference(filename, preset="sr"):
    return mappy.Aligner(filename, preset)

def change_bam2list(bam):
    segments_list = []
    for read in bam:
        read = MappedSegment(read, "pysam")
        segments_list.append(read)
    return segments_list

def change_bam2dict(bam):
    segments_dict = collections.defaultdict(list)
    for read in bam:
        read = MappedSegment(read, "pysam")
        segments_dict[read.name].append(read)
    return segments_dict

def read_in_reference_for_reads():
    global reference_for_reads
    if reference_for_reads is None:
        #if hasattr(arguments, 'reference'):
        if not arguments.reference is None:
            reference_for_reads = read_in_fasta_as_reference(arguments.reference)

def read_in_reference_for_contigs():
    global reference_for_contigs
    if reference_for_contigs is None:
        #if hasattr(arguments, 'reference'):
        if not arguments.reference is None:
            #reference_for_contigs = read_in_fasta_as_reference(arguments.reference, "asm10")
            reference_for_contigs = read_in_fasta_as_reference(arguments.reference, "sr")

def read_in_contigs_as_reference():
    global contigs_as_reference
    contigs_format = guess_fileformat(arguments.contigs)
    if contigs_format == "fasta":
        contigs_as_reference = mappy.Aligner(arguments.contigs)

#def read_in_file(filename):
#    fileformat = guess_fileformat(filename)
#    if fileformat == "bam":
#        return read_in_bam(filename)
#    elif fileformat == "fasta":
#        fasta = read_in_fasta_as_query(filename)
#        return map2reference(fasta)

def read_in_contigs2reference():
    global contigs2reference
    filename = arguments.contigs
    #contigs2reference = read_in_file(filename)
    fileformat = guess_fileformat(filename)
    if fileformat == "bam":
        contigs2reference = read_in_bam(filename)
    elif fileformat == "fasta":
        fasta = read_in_fasta_as_query(filename)
        contigs2reference = map_contigs2reference(fasta)

def read_in_reads2reference():
    global reads2reference
    filename = arguments.reads
    #reads2reference = read_in_file(filename)
    fileformat = guess_fileformat(filename)
    if fileformat == "bam":
        reads2reference = read_in_bam(filename)
    elif fileformat == "fasta":
        fasta = read_in_fasta_as_query(filename)
        reads2reference =  map_reads2reference(fasta)

def read_in_reads2contigs():
    global reads2contigs
    #if hasattr(arguments, 'reads2contigs'):
    if not arguments.reads2contigs is None:
        reads2contigs = read_in_bam(arguments.reads2contigs)
    else:
        reads = read_in_fasta_as_query(arguments.reads)
        read_in_contigs_as_reference()
        reads2contigs = map2contigs(reads)

def get_all_reads():
    global all_reads
    for read in reads2contigs.keys():
        all_reads.add(read)
    for read in reads2reference.keys():
        all_reads.add(read)
        
def main():
    global arguments
    global reads2contigs
    global reads2reference
    global contigs2reference
    global output
    arguments = parse_arguments()
    read_in_contigs2reference()
    print(contigs2reference)
    read_in_reads2reference()
    print(reads2reference)
    read_in_reads2contigs()
    print(reads2contigs)
    get_all_reads()
    
    #output = pysam.AlignmentFile(arguments.output, "w")
    output = open(arguments.output, "w")
    for read_name, mappings in reads2reference.items():
        if arguments.mapped2ref == "leave":
            for mapping in mappings:
                mapping.write_alignment(output)
    for read_name, mappings in reads2contigs.items():
        print(read_name)
        for mapping in mappings:
            print(mapping)
            contig_name = mapping.reference
            print(contig_name)
            contig_mappings = contigs2reference[contig_name]
            for contig_mapping in contig_mappings:
                print(contig_mapping)
                # transfer in place
                mapping.transfer_mapping(contig_mapping)
                # for testing purposes I write them all to file
                mapping.write_alignment(output)
        # here we decide what we want to do with these mappings
        # potentially we write one / some / all to output

if __name__ == '__main__':
    main()
