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

class MappedSegment:
    """
    Class to store mapped segment, like read or contig,
    unifying interface of mappy and pysam.

    Has the following attributes:

    ===========================================================================
    name            | mappy name |      pysam name       |relation between them
    ---------------------------------------------------------------------------
    strand          |   strand   |      is_reverse       |  ~
    is_primary      | is_primary |     is_secondary      |  !=
    reference       |    ctg     |    reference_name     |  =
    start_on_query  |    q_st    | query_alignment_start |  ?
    end_on_query    |    q_en    |  query_alignment_end  |  ?
    start_on_ref    |    r_st    |    reference_start    |  ?
    end_on_ref      |    r_en    |     reference_end     |  ?
    mapping_quality |    mapq    |    mapping_quality    |  =
    cigar           | cigar_str  |      cigarstring      |  =
    cigar_tuples    |   cigar    |      cigartuples      | reversed
    is_mapped       |     -      |      !is_unmapped     |
    name
    sequence
    ===========================================================================
    ctg_len? uzywam tego w nowej wersji IMCA, ale chyba tylko dlatego ze mappy tego potrzebowal
        wiec jesli nie bede operowac na jego obiektach to chyba spoko

    ------
    Attributes of mappy.Alignment:
    ------
    ctg: name of the reference sequence the query is mapped to
    ctg_len: total length of the reference sequence
    r_st and r_en: start and end positions on the reference
    q_st and q_en: start and end positions on the query
    strand: +1 if on the forward strand; -1 if on the reverse strand
    mapq: mapping quality
    blen: length of the alignment, including both alignment matches and gaps but excluding ambiguous bases.
    mlen: length of the matching bases in the alignment, excluding ambiguous base matches.
    NM: number of mismatches, gaps and ambiguous poistions in the alignment
    trans_strand: transcript strand. +1 if on the forward strand; -1 if on the reverse strand; 0 if unknown
    is_primary: if the alignment is primary (typically the best and the first to generate)
    read_num: read number that the alignment corresponds to; 1 for the first read and 2 for the second read
    cigar_str: CIGAR string
    cigar: CIGAR returned as an array of shape (n_cigar,2). The two numbers give the length and the operator of each CIGAR operation.
    MD: the MD tag as in the SAM format. It is an empty string unless the MD argument is applied when calling mappy.Aligner.map().
    cs: the cs tag.

    ------
    Attributes of pysam.AlignedSegment:
    ------
    bin: properties bin
    cigarstring: the cigar alignment as a string.
    cigartuples: the cigar alignment. The alignment is returned as a list of tuples of (operation, length).
    flag: properties flag
    is_duplicate: true if optical or PCR duplicate
    is_paired: true if read is paired in sequencing
    is_proper_pair: true if read is mapped in a proper pair
    is_qcfail: true if QC failure
    is_read1: true if this is read1
    is_read2: true if this is read2
    is_reverse: true if read is mapped to reverse strand
    is_secondary: true if not primary alignment
    is_supplementary: true if this is a supplementary alignment
    is_unmapped: true if read itself is unmapped
    mapping_quality: mapping quality
    mate_is_reverse: true is read is mapped to reverse strand
    mate_is_unmapped: true if the mate is unmapped
    next_reference_id: the reference id of the mate/next read.
    next_reference_name: reference name of the mate/next read (None if no AlignmentFile is associated)
    next_reference_start: the position of the mate/next read.
    query_alignment_end: end index of the aligned query portion of the sequence (0-based, exclusive) This the index just past the last base in seq that is not soft-clipped.
    query_alignment_length: length of the aligned query sequence. This is equal to qend - qstart
    query_alignment_qualities: aligned query sequence quality values (None if not present). These are the quality values that correspond to query, that is, they exclude qualities of soft clipped bases. This is equal to qual[qstart:qend]. Quality scores are returned as a python array of unsigned chars. Note that this is not the ASCII-encoded value typically seen in FASTQ or SAM formatted files. Thus, no offset of 33 needs to be subtracted.
    query_alignment_sequence: aligned portion of the read. This is a substring of seq that excludes flanking bases that were soft clipped (None if not present). It is equal to seq[qstart:qend].
    SAM/BAM files may include extra flanking bases that are not part of the alignment. These bases may be the result of the Smith-Waterman or other algorithms, which may not require alignments that begin at the first residue or end at the last. In addition, extra sequencing adapters, multiplex identifiers, and low-quality bases that were not considered for alignment may have been retained.
    query_alignment_start: start index of the aligned query portion of the sequence (0-based, inclusive). This the index of the first base in seq that is not soft-clipped.
    query_length: the length of the query/read. This value corresponds to the length of the sequence supplied in the BAM/SAM file. The length of a query is 0 if there is no sequence in the BAM/SAM file. In those cases, the read length can be inferred from the CIGAR alignment, see pysam.AlignedSegment.infer_query_length(). The length includes soft-clipped bases and is equal to len(query_sequence). This property is read-only but can be set by providing a sequence. Returns 0 if not available.
    query_name: the query template name (None if not present)
    query_qualities: read sequence base qualities, including soft clipped bases (None if not present). Quality scores are returned as a python array of unsigned chars.
    Note that this is not the ASCII-encoded value typically seen in FASTQ or SAM formatted files. Thus, no offset of 33 needs to be subtracted. Note that to set quality scores 
    the sequence has to be set beforehand as this will determine the expected length of the quality score array. This method raises a ValueError if the length of the quality scores
    and the sequence are not the same.
    query_sequence: read sequence bases, including soft clipped bases (None if not present). The sequence is returned as it is stored in the BAM file. Some mappers might have stored a reverse complement of the original read sequence.
    reference_end: aligned reference position of the read on the reference genome. reference_end points to one past the last aligned residue. Returns None if not available (read is unmapped or no cigar alignment present).
    reference_id: reference ID This field contains the index of the reference sequence in the sequence dictionary. To obtain the name of the reference sequence, use get_reference_name()
    reference_length: aligned length of the read on the reference genome. This is equal to aend - pos. Returns None if not available.
    reference_name: reference name
    reference_start: 0-based leftmost coordinate
    tags: deprecated, use get_tags() instead
    template_length: the observed query template length
    """

    def __init__(self, segment, which_format):
        if which_format == "mappy":
            self.read_from_mappy(segment)
        elif which_format == "pysam":
            self.read_from_pysam(segment)
        elif which_format == "new":
            self.create_from_scratch(segment)
        else:
            raise ValueError("Unknown format: %s;"
                             " accepted formats are mappy, pysam,"
                             " and dict of attributes." % which_format)

    def read_from_mappy(self, mappy_alignment):
        if mappy_alignment.strand == 1:
            self.strand = '+'
        elif mappy_alignment.strand == -1:
            self.strand = '-'
        else:
            self.strand = '.'
        self.is_primary = mappy_alignment.is_primary
        self.reference = mappy_alignment.ctg
        # coordinates: to check
        self.start_on_query = mappy_alignment.q_st
        self.end_on_query = mappy_alignment.q_en
        self.start_on_ref = mappy_alignment.r_st
        self.end_on_ref = mappy_alignment.r_en
        self.mapping_quality = mappy_alignment.mapq
        self.cigar = mappy_alignment.cigar_str
        self.cigar_tuples = mappy_alignment.cigar
        self.is_mapped = True
        self.name = "."

    def read_from_pysam(self, pysam_aligned_segment):
        if pysam_aligned_segment.is_reverse== True:
            self.strand = '-'
        else:
            self.strand = '+'
        self.is_primary = not pysam_aligned_segment.is_secondary
        self.reference = pysam_aligned_segment.reference_name
        # coordinates: to check
        self.start_on_ref = pysam_aligned_segment.query_alignment_start
        self.end_on_ref = pysam_aligned_segment.query_alignment_end
        self.start_on_query = pysam_aligned_segment.query_alignment_start
        self.end_on_query = pysam_aligned_segment.query_alignment_end
        self.mapping_quality = pysam_aligned_segment.mapping_quality
        self.cigar = pysam_aligned_segment.cigarstring
        self.cigar_tuples = [(i[1], i[0]) for i in pysam_aligned_segment.cigartuples]
        self.is_mapped = not pysam_aligned_segment.is_unmapped
        self.name = pysam_aligned_segment.query_name

    def create_from_scratch(self, attributes):
        self.strand = attributes['strand']
        self.is_primary = attributes['is_primary']
        self.reference = attributes['reference']
        self.start_on_query = attributes['start_on_query']
        self.end_on_query = attributes['end_on_query']
        self.start_on_ref = attributes['start_on_ref']
        self.end_on_ref = attributes['end_on_ref']
        self.mapping_quality = attributes['mapping_quality']
        self.cigar = attributes['cigar']
        self.cigar_tuples = attributes['cigar_tuples']
        self.is_mapped = attributes['is_mapped']
        self.name = attributes['name']

    def calculate_flag(self):
        """
        From SAM Format Specification
        https://samtools.github.io/hts-specs/SAMv1.pdf
        version from 20.06.2019:
        =============================================================================
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
        =============================================================================
        """
        flag = 0
        if self.strand == '-':
            flag += 16
        if not self.is_primary:
            flag += 256
        return str(flag)

    def alignment2sam(self):
        """
        Return a line from SAM file.

        From SAM Format Specification
        https://samtools.github.io/hts-specs/SAMv1.pdf
        version from 20.06.2019:
        =================================================================================
        1  QNAME String [!-?A-~]{1,254}             Query template NAME
        2  FLAG  Int 
        [0, 2^16-1] 
        bitwise
        FLAG
        3  RNAME String \*|[:rname:^*=][:rname:]*   Reference sequence NAME
        4  POS   Int    [0/2^31-1]                  1-based leftmost mapping POSition
        5  MAPQ  Int   
        [0,2^8-1]       
        MAPping
        Quality
        6  CIGAR String \*|([0-9]+[MIDNSHPX=])+     CIGAR string
        7  RNEXT String \*|=|[:rname:^*=][:rname:]* Reference name of the mate/next read
        8  PNEXT Int    [0,2^31-1]                  Position of the mate/next read
        9  TLEN  Int    [-2^31+1, 2^31-1]           observed Template LENgth
        10 SEQ   String \*|[A-Za-z=.]+              segment SEQuence
        11 QUAL  String [!-~]+                      ASCII of Phred-scaled base QUALity+33
        =================================================================================
        """
        # TODO: fields 7, 8, 9, potentially 11
        flag = self.calculate_flag()
        if hasattr(self, "sequence"):
            sequence = self.sequence
        else:
            sequence = "."
        return '\t'.join([self.name, flag,
                          self.reference, str(self.start_on_ref),
                          str(self.mapping_quality), self.cigar,
                          '*', '0', '0', sequence, '*'])

    def __str__(self):
        return self.alignment2sam()

    def write_alignment(self, output):
        """
        Write an alignment to file output in SAM format.
        """
        output.write(self.alignment2sam())
        output.write('\n')

    def calculate_offset(self, coordinate):
        """
        Given a CIGAR describing mapping of seq A to seq B,
        calculate the offset of coordinates
        needed to transfer coordinate from A to B,
        assuming mapping of A starts at start
        (that is, start-th coordinate of A is the beggining of CIGAR description).

        For example:
        If there are no indels in CIGAR,
        offset should be equal zero.

        If there is one deletion, offset will be positive:

        CIGAR: 10M5D20M ([[10, 0], [5, 2], [20, 0]])
        start: 3
        coordinate: 30

        Offset: +5

        If there is one insertion, offset will be negative.


              3M  3I    9M                 CIGAR
           ________________x____           seqA
             |   |  /          /
             |   | /          /            mapping
        _____|___|/__________/_________    seqB

        In this case:
        CIGAR: 3M3I9M
        start: 2
        coordinate: 15

        Offset: -3

        It means, if we want to transfer cooridnate 15 (denoted by x here) to seqB
        we should get the start coordinate of mapping seqA to seqB,
        add 15 and decrease by 3.
        """
        position_on_contig = self.start_on_query
        offset = - self.start_on_query
        for count, operation in self.cigar_tuples:
            #print("Count, operation:")
            #print(count, operation)
            end_of_operation = position_on_contig
            if operation in OPERATIONS_QUERY_CONSUMING:
                end_of_operation += count
            #print("end of operation:")
            #print(end_of_operation)
            if end_of_operation >= coordinate:
                #print("koniec")
                break
            if operation in OPERATIONS_REFERENCE_CONSUMING:
                offset += count
            if operation in OPERATIONS_QUERY_CONSUMING:
                position_on_contig += count
                offset -= count
        #print(position_on_contig, offset, coordinate)
        return offset
 
    def transfer_mapping(self, contig2ref):
        """
        contig2ref - MappedSegment with contig mapped to reference

        Transfers mapping in place.
        """
        # TODO:
        # check reverse mapping -> 
        #  result is slightly different when I transfer through contigs.fa and contigs_rev.fa,
        #  it probably shouldn't (98, 162 -> 93, 163; 329, 399 -> 331, 399)
        # cigar? Do we keep the original one? Or merge with contig2ref?
        # Mapping quality - how do we want to calculate it?
        #  we could keep the read2contig one
        # primary mappings?

        #print("Read2contig:")
        #print(self)
        #print("Contig2ref:")
        #print(contig2ref)

        if contig2ref.strand == '+':
            start_offset = contig2ref.calculate_offset(self.start_on_ref)
            end_offset = contig2ref.calculate_offset(self.end_on_ref)
            start = self.start_on_ref + contig2ref.start_on_ref + start_offset
            end = self.end_on_ref + contig2ref.start_on_ref + end_offset
        else:
            contig_length = sum([operation[0] for operation in contig2ref.cigar_tuples if operation[1] in OPERATIONS_QUERY_CONSUMING]) # to chyba nie jest poprawne, ale tymczasowo niech bedzie
            # TODO dokonczyc poprawianie tego else'a, teraz nie zadziala
            # nie wiem czy nie musze zsecodowac sprawdzania nici na calculate_offset
            cigar = contig2ref.cigar_tuples
            cigar.reverse()
            start_offset = calculate_offset(cigar,
                                            contig_length - contig2ref.q_en,
                                            read2contig.r_st)
            end_offset = calculate_offset(cigar,
                                          contig_length - contig2ref.q_en,
                                          contig_length - read2contig.r_st)
            start = (contig_length - read2contig.r_en) + \
                    contig2ref.r_st + end_offset
            end = (contig_length - read2contig.r_st) + \
                  contig2ref.r_st + start_offset
        #end = start + read2contig.r_en - read2contig.r_st
        # ^ na pewno mozna zrobic lepiej mapowanie koncowego koordynatu
        # tylko ze ja chyba juz nie potrzebuje atrybutu end
        #  byl mi potrzebny jak operowalam na mappy.Alignment
        print("start, end:")
        print(start, end)
        if self.strand == contig2ref.strand:
            self.strand = '+'
        else:
            self.strand = '-'
        self.reference = contig2ref.reference
        self.start_on_ref = start
        self.end_on_ref = end
        print(self)
