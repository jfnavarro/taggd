""" 
This class inherits from Record and it represents
a FASTA or FASTQ record
"""
from taggd.io.record import Record

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
        self.annotation = fqr[0]
        self.sequence = fqr[1]
        if len(fqr) == 3: self.attributes["quality"] = fqr[2]
        self.attributes["taggdtags"] = None

    def add_tags(self, list added):
        """
        Appends extra tags for taggd.
        :param added a list of tag tuples (name,value)
        """
        self.attributes["taggdtags"] = added

    def unwrap(self, force_fasta=False):
        """
        Returns (annotation, sequence) or (annotation, sequence, quality) if quality is available
        """
        cdef str k
        cdef object v
        cdef str taggdtags_str = ' '.join(["{}:{}".format(k,v) for k,v in self.attributes["taggdtags"]]) if self.attributes["taggdtags"] else ""
        if "quality" in self.attributes and not force_fasta:
            return ("{} {}".format(self.annotation, taggdtags_str), 
                    self.sequence, self.attributes["quality"])
        else:
            return ("{} {}".format(self.annotation, taggdtags_str), self.sequence)

    def __str__(self):
        """
        String representation
        """
        cdef str fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
        cdef str fa_format = '>{header_comments}\n{sequence}\n'
        if "quality" in self.attributes:
            return fq_format.format(header_comments=self.annotation, 
                                    sequence=self.sequence, quality=self.attributes["quality"])
        else:
            return fa_format.format(header_comments=self.annotation, sequence=self.sequence)