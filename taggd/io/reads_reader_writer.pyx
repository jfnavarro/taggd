"""
Interface for writing/reading FASTQ/FASTA and SAM files
The idea is to write always the abstract class called Record
"""
import pysam as ps
import sys
import os
from cpython cimport bool
import taggd.io.fastq_utils as fu
from taggd.io.record import *

# Global variables for match types.
cdef int FASTQ, FASTA, SAM, BAM

class ReadsReaderWriter():
    """
    Provides opening of a reads file, and writing to the same format to one or more
    specified files.
    """

    def __init__(self, str reads_infile_name):
        """
        Constructor
        """

        self.file_type = -1
        self.output_file_format = -1
        self.infile_name = reads_infile_name
        self.infile = None
        self.infile_header = None
        self.second_fastq_filename = None

        # Supported file types
        global FASTQ
        global FASTA
        global SAM
        global BAM
        FASTQ, FASTA, SAM, BAM = range(4)

        # Extract the file type
        cdef str suffix = os.path.splitext(self.infile_name)[1].lower()
        if suffix == ".fa" or suffix == ".fasta":
            self.file_type = FASTA
        elif suffix == ".fq" or suffix == ".fastq":
            self.file_type = FASTQ
        elif suffix == ".sam":
            self.file_type = SAM
        elif suffix == ".bam":
            self.file_type = BAM
        else:
            raise ValueError("Unsupported reads file format!")

        # Read header.
        if self.file_type == SAM or self.file_type == BAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", 
                                           check_header=True, check_sq=False)
            self.infile_header = self.infile.header
            self.infile.close()
        elif self.file_type == FASTA or self.file_type == FASTQ:
            self.infile_header = { 'HD': {'VN': '1', 'SO':'unsorted'} }
            #self.infile_header = { 'HD': {'VN': '1', 'GO':'query'} }

    def set_second_fastq_filename(self, second_fastq_filename):
        assert self.file_type == FASTQ, 'ERROR: paried infiles are only allowed when using FASTQ format.'
        self.second_fastq_filename = second_fastq_filename
        
    def reader_open(self):
        """
        Opens the reads file using appropriate format.
        """
        # Ensure to close
        self.reader_close()
        
        # Open file.
        if self.file_type == FASTA or self.file_type == FASTQ:
            if self.second_fastq_filename and self.file_type == FASTQ:
                import itertools
                self._infile_1 = open(self.infile_name, "r")
                self._infile_2 = open(self.second_fastq_filename, "r")
                self.infile = itertools.izip_longest(
                        fu.readfq(self._infile_1),
                        fu.readfq(self._infile_2)
                    )
            else:
                self.infile = fu.readfq(open(self.infile_name, "r"))
        elif self.file_type == SAM:
            self.infile = ps.AlignmentFile(self.infile_name, "r", check_header=True, check_sq=False)
        elif self.file_type == BAM:
            self.infile = ps.AlignmentFile(self.infile_name, "rb", check_header=True, check_sq=False)
        else:
            raise ValueError("Unsupported reads file format!")

        cdef object rec

        for orig in self.infile:
            rec = Record(orig)
            yield rec

    def reader_close(self):
        """
        Closes the input file handler.
        """
        import itertools
        if self.infile != None:
            if isinstance( self.infile, itertools.izip_longest ):
                self._infile_1.close()
                self._infile_2.close()
            else:
                self.infile.close()
                self.infile = None

    def __exit__(self, type, value, tb):
        """
        Always close the input file.
        """
        self.reader_close()

    def get_format(self):
        """
        Returns the file format.
        """
        if self.file_type == FASTA:
            return "fa"
        if self.file_type == FASTQ:
            return "fq"
        if self.file_type == SAM:
            return "sam"
        if self.file_type == BAM:
            return "bam"
        return None

    def get_writer(self, str outfile_name):
        """
        Returns a writer.
        :param outfile_name the name of the file to create
        :return the file handler so records can be written on it.
        """
        # Remove if exists
        if os.path.exists(outfile_name):
            os.remove(outfile_name)

        # get the output file type
        cdef str suffix = outfile_name.split(".")[-1]
        if suffix == "fa":    self.output_file_format = FASTA
        elif suffix == "fq":  self.output_file_format = FASTQ
        elif suffix == "sam": self.output_file_format = SAM
        elif suffix == "bam": self.output_file_format = BAM
        else: raise ValueError("Unsupported output file format!")
        if self.file_type == FASTA and self.output_file_format != FASTA: raise ValueError("Unsupported output file format!\n            No quality values present.\n            FASTA is the only output format allowed for FASTA input.")

        # Open the file and returns the handler
        if self.output_file_format == FASTA or self.output_file_format == FASTQ:
            return open(outfile_name, "w")
        elif self.output_file_format == SAM:
            return ps.AlignmentFile(outfile_name, "wh", header=self.infile_header)
        elif self.output_file_format == BAM:
            return ps.AlignmentFile(outfile_name, "wb", header=self.infile_header)
        else:
            raise ValueError("Unknown file format for writer")

    def write_record(self, outfile, record):
        """
        Writes a record in the filename descriptor given.
        Important out_handler must be a descriptor opened (using get_writer())
        :param out_handler the outpuf file handler
        :param record the Record object to write
        """
        #TODO record could not be the same type of the file handler(check this)
        for read in record.unwrap(output_format=self.output_file_format):
            if   self.output_file_format == FASTA: fu.writefa_record(outfile, read)
            elif self.output_file_format == FASTQ: fu.writefq_record(outfile, read)
            elif self.output_file_format in [SAM, BAM]: outfile.write(read)
            else:
                raise ValueError("Unknown file format for record")