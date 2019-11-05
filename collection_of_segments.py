import collections
from mapped_segment import MappedSegment

class CollectionOfSegments(object):

    def __init__(self):
        self.collection = collections.defaultdict(list)

    def read_in_from_bam(self, filename):
        import pysam
        bamfile = pysam.AlignmentFile(filename)
        for read in bamfile:
            read = MappedSegment(read, "pysam")
            self.add_read(read)

    def map_reads(self, reads, aligner):
        """
        reads - generator of reads, as created by mappy
        aligner - mappy.Aligner object
        """
        for read in reads:
            name, seq, quality = read
            mappings = aligner.map(seq)
            for mapping in mappings:
                mapping = mappedSegment(mapping, "mappy")
                mapping.name = name
                self.add_mapping(mapping)

    def is_empty(self):
        return len(self.collection) == 0

    def add_mapping(self, mapping):
        self.collection[mapping.name].append(mapping)

    def write_sam(self, output_name):
        #bam = pysam.AlignmentFile(outputname, "wb")
        sam = open(output_name, "w")
        self.write_header(sam)
        for mappings in self.collection.values():
            for mapping in mappings:
                mapping.write_alignment(sam)
        sam.close()

    def write_header(self, output):
        """
        To write header I would need to have chromosome lengths, not only names.
        So I'm not sure how exactly to go about this.
        """
        pass
