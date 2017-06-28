"""
A simple class that is meant to mimic a Pysam record
and add some extra attributes so it can later be
written to a sam file
"""
import pysam
from taggd.io.record import Record

FASTQ, FASTA, SAM, BAM = range(4)

class SAMRecord(Record):
    """
    Holds a SAM record.
    Follows pysam's representation.
    """

    def __init__(self, object alseq):
        """
        Constructor. 
        Represents a pysam.AlignedSequence as a Record.
        """
        Record.__init__(self)
        assert isinstance(alseq, pysam.AlignedSegment), 'ERROR: this is not a pysam.AlignedSegment'
        self.input_format = SAM
        self.annotation = str(alseq.query_name)
        self.sequence = str(alseq.query_sequence)
        self.attributes["flag"] = alseq.flag
        self.attributes["reference_id"] = alseq.reference_id
        self.attributes["reference_start"] = alseq.reference_start
        self.attributes["mapping_quality"] = alseq.mapping_quality
        self.attributes["cigar"] = alseq.cigar
        self.attributes["next_reference_id"] = alseq.next_reference_id
        self.attributes["reference_start"] = alseq.reference_start
        self.attributes["template_length"] = alseq.template_length
        self.attributes["query_qualities"] = alseq.query_qualities
        self.attributes["tags"] = alseq.tags
        self.attributes["taggdtags"] = None

    def add_tags(self, list added):
        """
        Appends extra tags for taggd.
        :param added a list of tag tuples (name,value)
        """
        self.attributes["taggdtags"] = added

    @property
    def taggdtags_str(self):
        cdef str k
        cdef object v
        cdef str taggdtags_str = ' '.join(["{}:{}".format(k,v) for k,v in self.attributes["taggdtags"]]) if self.attributes["taggdtags"] else ""
        return taggdtags_str

    def unwrap(self, output_format=None):
        """
        Returns a pysam alignment.
        """
        cdef object a = pysam.AlignedSegment()
        if output_format in [None, SAM, BAM]:
            a.query_name = self.annotation
            a.flag = self.attributes["flag"]
            a.reference_id = self.attributes["reference_id"]
            a.reference_start = self.attributes["reference_start"]
            a.mapping_quality = self.attributes["mapping_quality"]
            a.cigar = self.attributes["cigar"]
            a.next_reference_id = self.attributes["reference_id"]
            a.next_reference_start = self.attributes["reference_start"]
            a.template_length = self.attributes["template_length"]
            a.query_sequence = self.sequence
            a.query_qualities = self.attributes["query_qualities"]
            a.tags = self.attributes["tags"]
            if self.attributes["taggdtags"] != None:
                    a.tags += self.attributes["taggdtags"]
            return a
        elif output_format == FASTQ: return self.unwrap_fastq()
        elif output_format == FASTA: return self.unwrap_fasta()

    
    def unwrap_fasta(self):
        """
        Returns a (annotation, sequence) tuple
        """
        return ("{} {}".format(self.annotation, self.taggdtags_str), self.sequence)
        
    def unwrap_fastq(self):
        """
        Returns a (annotation, sequence, quality) tuple
        """
        return ("{} {}".format(self.annotation, self.taggdtags_str), self.sequence, self.attributes["query_qualities"])