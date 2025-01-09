"""
This class inherits from Record and it represents
a FASTA record
"""
from typing import List, Tuple
from taggd.io.record import Record


class FASTARecord(Record):
    """
    Holds a FASTA record.
    """

    def __init__(self, header: str, sequence: str):
        """
        Initializes a FASTARecord from its main attributes.

        Args:
            header: The FASTA header.
            sequence: The FASTA sequence.
        """
        super().__init__()
        self.header = header
        self.sequence = sequence
        self.taggdtags: str = ""

    def add_tags(self, added: List[Tuple[str, str]]) -> None:
        """
        Appends tags for extra information.

        Args:
            added: A list of tag tuples (name, value).
        """
        self.taggdtags = " ".join([f"{k}:{v}" for k, v in added])

    def unwrap(self) -> Tuple[str, str]:
        """
        Returns the FASTA record as a tuple (annotation, sequence).

        Returns:
            The FASTA record with tags included in the annotation.
        """
        return (f"{self.header} {self.taggdtags}", self.sequence)

    def __str__(self) -> str:
        """
        Returns the string representation of the FASTA record in standard FASTA format.

        Returns:
            The FASTA record as a formatted string.
        """
        fa_format = ">{header_comments}\n{sequence}\n"
        return fa_format.format(header_comments=self.header, sequence=self.sequence)
