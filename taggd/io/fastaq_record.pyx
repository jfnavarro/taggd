""" 
This class inherits from Record and it represents
a FASTA or FASTQ record
"""
from taggd.io.record import Record
import pysam

FASTQ, FASTA, SAM, BAM = range(4)

class FASTAQRecord(Record):
    """
    Holds a FASTA or FASTQ record. 
    Not cdef-cythonized so as to be threading compatible.
    """

    def __init__(self, object fqr):
        """
        Constructor. Represents a (annotation, seq[, qual]) tuple as a Record.
        """
        Record.__init__(self)
        
        if isinstance(fqr, tuple):
            if len(fqr) == 2: self._init_fasta(fqr)
            elif len(fqr) == 3: self._init_fastq(fqr)

        self.attributes["taggdtags"] = None

    def _init_fasta(self, fqr):
        """
        Format specific Constructor.
        Represents a (annotation, seq) tuple as record.
        """
        self.input_format = FASTA
        self.annotation = fqr[0]
        self.sequence = fqr[1]

    def _init_fastq(self, fqr):
        """
        Format specific Constructor.
        Represents a (annotation, seq, qual) tuple as record.
        """
        self._init_fasta(fqr)
        self.input_format = FASTQ
        self.attributes["query_qualities"] = fqr[2]

    def add_tags(self, list added):
        """
        Appends extra tags for taggd.
        :param added a list of tag tuples (name,value)
        """
        self.attributes["taggdtags"] = added

    @property
    def taggdtags_str(self):
        """
        Propery that returns a string representation of the self.attributes["taggdtags"] list
        """
        cdef str k
        cdef object v
        cdef str taggdtags_str = ' '.join(["{}:{}".format(k,v) for k,v in self.attributes["taggdtags"]]) if self.attributes["taggdtags"] else ""
        return taggdtags_str
    
    def unwrap(self, output_format=None):
        """
        Returns a (annotation, sequence) tuple, a (annotation, sequence, quality) tuple or a pysam.AlignedSegment depending on the output_format variable
        """
        
        if output_format == None: output_format = self.input_format
        if "query_qualities" not in self.attributes: output_format = FASTA
        
        if output_format == FASTQ:   return self.unwrap_fastq()
        elif output_format == FASTA: return self.unwrap_fasta()
        elif output_format == BAM or output_format == SAM:
            return self.unwrap_sam()
    
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
    
    def unwrap_sam(self):
        """
        Returns a pysam alignment.
        """
        cdef object a = pysam.AlignedSegment()
        a.query_name = self.annotation
        a.query_sequence = self.sequence

        if "query_qualities" in self.attributes: a.query_qualities = self.attributes["query_qualities"]
        
        if self.attributes["taggdtags"] != None:
                a.tags += self.attributes["taggdtags"]
        
        return a

    def __str__(self):
        """
        String representation
        """
        cdef str fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
        cdef str fa_format = '>{header_comments}\n{sequence}\n'
        if   self.input_format == FASTQ: return fq_format.format(header_comments=self.annotation, sequence=self.sequence, quality=self.attributes["query_qualities"])
        elif self.input_format == FASTA: return fa_format.format(header_comments=self.annotation, sequence=self.sequence)
        elif self.input_format == SAM:   return self.unwrap_sam().__str__()