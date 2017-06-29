""" 
This class represents
a FASTA or FASTQ record or pysam.AlignedSegment
"""
import pysam

FASTQ, FASTA, SAM, BAM, FASTQPE = range(5)

class Record(object):
    """
    Holds a FASTA or FASTQ record. 
    Not cdef-cythonized so as to be threading compatible.
    """

    def __init__(self, object req):
        """
        Constructor. Represents a (annotation, seq[, qual]) tuple or a pysam.AlignedSegment as a Record.
        """
        self.annotation = None
        self.sequence = None
        self.attributes = dict()
        
        if isinstance(req, tuple):
            if len(req) == 2:
                if   isinstance(req[0], str):   self._init_fasta(req)
                elif isinstance(req[0], tuple): self._init_paired_fastq(req)
            elif len(req) == 3: self._init_fastq(req)

        elif isinstance(req, pysam.AlignedSegment):
            self._init_sam(req)
        
        self.attributes["taggdtags"] = None

    def _init_fasta(self, req):
        """
        Format specific Constructor.
        Represents a (annotation, seq) tuple as record.
        """
        self.input_format = FASTA
        self.annotation = req[0]
        self.sequence = req[1]

    def _init_fastq(self, req):
        """
        Format specific Constructor.
        Represents a (annotation, seq, qual) tuple as record.
        """
        self._init_fasta(req)
        self.input_format = FASTQ
        self.attributes["query_qualities"] = req[2]

    def _init_paired_fastq(self, req):
        """
        Format specific Constructor.
        Represents a ((annotation, seq, qual),(annotation, seq, qual)) nested tuple as record.
        """
        self._init_fastq( req[0] )
        self.input_format = FASTQPE
        self.attributes["second_read_annotation"]      = req[1][0]
        self.attributes["second_read_sequence"]        = req[1][1]
        self.attributes["second_read_query_qualities"] = req[1][2]
        assert self.annotation.split(' ')[0] == self.attributes["second_read_annotation"].split(' ')[0], "ERROR: order of fastq headers don't match in paired infiles!\n\tINFO: "+self.annotation.split(' ')[0]+" != "+self.attributes["second_read_annotation"].split(' ')[0]+"\n"

    def _init_sam(self, req):
        """
        Format specific Constructor.
        Represents a pysam.AlignedSegment as record.
        """
        self.input_format = SAM
        self.annotation = str(req.query_name)
        self.sequence = str(req.query_sequence)
        self.attributes["flag"] = req.flag
        self.attributes["reference_id"] = req.reference_id
        self.attributes["reference_start"] = req.reference_start
        self.attributes["mapping_quality"] = req.mapping_quality
        self.attributes["cigar"] = req.cigar
        self.attributes["next_reference_id"] = req.next_reference_id
        self.attributes["reference_start"] = req.reference_start
        self.attributes["template_length"] = req.template_length
        self.attributes["query_qualities"] = req.query_qualities
        self.attributes["tags"] = req.tags

    def add_tags(self, list added):
        """
        Appends extra tags for taggd.
        :param added a list of tag tuples (name,value)
        """
        self.attributes["taggdtags"] = added

    def get_umi_tag(self, tuple umi_coordinates):
        """
        Returns a tuple with a umi tag ("B3:Z",umi_sequence)
        """
        cdef int umi_start_position
        cdef int umi_end_position
        cdef str umi_sequence
        cdef str name = "B3:Z"
        
        umi_start_position, umi_end_position = umi_coordinates
        
        if len(self.sequence) < umi_end_position: raise ValueError("Read length is shorter than the umi end position, try rerunning with other or no umi coordinates.")
        else:
            umi_sequence = self.sequence[umi_start_position:umi_end_position]
            return (name,umi_sequence)

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
        if self.input_format == FASTQPE:
            return (
                    ("{} {}".format(self.annotation, self.taggdtags_str), self.sequence),
                    ("{} {}".format(self.attributes["second_read_annotation"], self.taggdtags_str), self.attributes["second_read_sequence"])
                   )
        else:
            return (("{} {}".format(self.annotation, self.taggdtags_str), self.sequence), )
        
    def unwrap_fastq(self):
        """
        Returns a (annotation, sequence, quality) tuple
        """
        if self.input_format == FASTQPE:
            return (
                    ("{} {}".format(self.annotation, self.taggdtags_str), self.sequence, self.attributes["query_qualities"]),
                    ("{} {}".format(self.attributes["second_read_annotation"], self.taggdtags_str), self.attributes["second_read_sequence"], self.attributes["second_read_query_qualities"] )
                   )
        else:
            return (("{} {}".format(self.annotation, self.taggdtags_str), self.sequence, self.attributes["query_qualities"]), )
    
    def unwrap_sam(self):
        """
        Returns a pysam.AlignedSegment.
        """
        cdef object a = pysam.AlignedSegment()
        cdef object b = pysam.AlignedSegment() # only used if input format is FASTQPE
        a.query_name = self.annotation
        a.query_sequence = self.sequence

        if "flag"            in self.attributes: a.flag                 = self.attributes["flag"]
        if "reference_id"    in self.attributes: a.reference_id         = self.attributes["reference_id"] 
        if "reference_start" in self.attributes: a.reference_start      = self.attributes["reference_start"] 
        if "mapping_quality" in self.attributes: a.mapping_quality      = self.attributes["mapping_quality"]
        if "cigar"           in self.attributes: a.cigar                = self.attributes["cigar"]
        if "reference_id"    in self.attributes: a.next_reference_id    = self.attributes["reference_id"] 
        if "reference_start" in self.attributes: a.next_reference_start = self.attributes["reference_start"]
        if "template_length" in self.attributes: a.template_length      = self.attributes["template_length"]
        if "query_qualities" in self.attributes: a.query_qualities      = self.attributes["query_qualities"]
        if "tags"            in self.attributes: a.tags                 = self.attributes["tags"]
            
        if self.attributes["taggdtags"] != None:
                a.tags += self.attributes["taggdtags"]
        
        if self.input_format == FASTQPE:
            b.query_name      = self.attributes["second_read_annotation"]
            b.query_sequence  = self.attributes["second_read_sequence"]
            b.query_qualities = self.attributes["second_read_query_qualities"]
            if self.attributes["taggdtags"] != None: b.tags += self.attributes["taggdtags"]
            return (a,b)
        else:
            return (a, )

    def __str__(self):
        """
        String representation
        """
        cdef str fq_format = '@{header_comments}\n{sequence}\n+\n{quality}\n'
        cdef str fa_format = '>{header_comments}\n{sequence}\n'
        if   self.input_format == FASTQ: return fq_format.format(header_comments=self.annotation, sequence=self.sequence, quality=self.attributes["query_qualities"])
        elif self.input_format == FASTA: return fa_format.format(header_comments=self.annotation, sequence=self.sequence)
        elif self.input_format == SAM:   return self.unwrap_sam().__str__()