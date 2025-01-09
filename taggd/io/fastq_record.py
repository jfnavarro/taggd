"""
This class inherits from Record and it represents
a FASTQ record
"""
from typing import List, Tuple
from taggd.io.record import Record


class FASTQRecord(Record):
    """
    Holds a FASTQ record.
    """

    def __init__(self, header: str, sequence: str, qualities: str):
        """
        Initializes a FASTQRecord from the main attributes.

        Args:
            header: The FASTQ header.
            sequence: The FASTQ sequence.
            qualities: The FASTQ qualities.
        """
        super().__init__()
        self.header = header
        self.sequence = sequence
        self.qualities = qualities
        self.taggdtags: str = ""

    def add_tags(self, added: List[Tuple[str, str]]) -> None:
        """
        Appends tags for extra information.

        Args:
            added: A list of tag tuples (name, value).
        """
        self.taggdtags = " ".join([f"{k}:{v}" for k, v in added])

    def unwrap(self) -> Tuple[str, str, str]:
        """
        Returns the FASTQ record as a tuple (header, sequence, quality).

        Returns:
            The FASTQ record with tags included in the header.
        """
        return (f"{self.header} {self.taggdtags}", self.sequence, self.qualities)

    def __str__(self) -> str:
        """
        Returns the string representation of the FASTQ record in standard FASTQ format.

        Returns:
            The FASTQ record as a formatted string.
        """
        fq_format = "@{header_comments}\n{sequence}\n+\n{quality}\n"
        return fq_format.format(
            header_comments=self.annotation,
            sequence=self.sequence,
            quality=self.qualities,
        )
