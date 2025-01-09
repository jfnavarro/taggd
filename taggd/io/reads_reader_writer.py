"""
Interface for writing/reading FASTQ/FASTA and SAM/BAM files
The idea is to write always the abstract class called Record
"""
import pysam
import os
from typing import Generator, Union
from taggd.io.sam_record import *
from taggd.io.fasta_record import *
from taggd.io.fastq_record import *
import dnaio


class ReadsReaderWriter:
    """
    Provides functionality for reading and writing reads files in various formats.
    Supports FASTA, FASTQ, SAM, and BAM formats.
    """

    FASTQ, FASTA, SAM, BAM = range(4)

    def __init__(self, reads_infile_name: str) -> None:
        """
        Initializes the reader/writer with the input file.

        Args:
            reads_infile_name (str): Path to the input reads file.
        """
        self.file_type = -1
        self.infile_name = reads_infile_name
        self.infile = None
        self.infile_header = None

        # Determine file type based on file extension
        suffix = os.path.splitext(self.infile_name)[1].lower()
        if suffix in [".fa", ".fasta"]:
            self.file_type = self.FASTA
        elif suffix in [".fq", ".fastq"]:
            self.file_type = self.FASTQ
        elif suffix == ".sam":
            self.file_type = self.SAM
        elif suffix == ".bam":
            self.file_type = self.BAM
        else:
            raise ValueError("Unsupported reads file format!")

        # Read header for SAM/BAM files
        if self.file_type in {self.SAM, self.BAM}:
            with pysam.AlignmentFile(self.infile_name, "r") as infile:
                self.infile_header = infile.header

    def reader_open(
        self,
    ) -> Generator[Union[dnaio.Record, pysam.AlignedSegment], None, None]:
        """
        Opens the reads file for reading and yields records.

        Yields:
            Records from the input file.
        """
        self.reader_close()

        if self.file_type == self.FASTA:
            self.infile = dnaio.FastaReader(self.infile_name)
        elif self.file_type == self.FASTQ:
            self.infile = dnaio.FastqReader(self.infile_name)
        elif self.file_type == self.SAM:
            self.infile = pysam.AlignmentFile(self.infile_name, "r")
        elif self.file_type == self.BAM:
            self.infile = pysam.AlignmentFile(self.infile_name, "rb")
        else:
            raise ValueError("Unsupported reads file format!")

        for record in self.infile:
            yield record

    def reader_close(self) -> None:
        """
        Closes the input file handler if it is open.
        """
        if self.infile is not None:
            self.infile.close()
            self.infile = None

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """
        Ensures the input file is closed on exit.
        """
        self.reader_close()

    def get_format(self) -> str:
        """
        Returns the format of the input file.

        Returns:
            The format as a string ("fa", "fq", "sam", or "bam").
        """
        if self.file_type == self.FASTA:
            return "fa"
        if self.file_type == self.FASTQ:
            return "fq"
        if self.file_type == self.SAM:
            return "sam"
        return "bam"

    def get_writer(
        self, outfile_name: str
    ) -> Union[dnaio.FastaWriter, dnaio.FastqWriter, pysam.AlignmentFile]:
        """
        Returns a writer for the specified output file.

        Args:
            outfile_name: Path to the output file.

        Returns:
            The writer object.
        """
        if os.path.exists(outfile_name):
            os.remove(outfile_name)

        if self.file_type == self.FASTA:
            return dnaio.FastaWriter(outfile_name)
        elif self.file_type == self.FASTQ:
            return dnaio.FastqWriter(outfile_name)
        elif self.file_type == self.SAM:
            return pysam.AlignmentFile(outfile_name, "wh", header=self.infile_header)
        elif self.file_type == self.BAM:
            return pysam.AlignmentFile(outfile_name, "wb", header=self.infile_header)
        else:
            raise ValueError("Unknown file format for writer")
