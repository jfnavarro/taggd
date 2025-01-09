"""
Base classes for different records that will be written to files
"""

from typing import Any, Dict, List, Optional, Tuple

class Record:
    """
    Represents a generic record for various file formats.

    This class serves as a universal object for representing records, with utility
    classes handling specific file format conversions.
    """

    def __init__(self) -> None:
        """
        Initializes a Record instance.
        """
        self.annotation = None
        self.sequence = None
        self.attributes = {}

    def unwrap(self) -> Optional[Any]:
        """
        Converts the record to its original representation.

        Returns:
            The original representation of the record, if applicable.
        """
        return None

    def add_tags(self, added: List[Tuple[str, Any]]) -> None:
        """
        Adds tags to the record as key-value pairs.

        Args:
            added: A list of tuples where each tuple contains
            a tag name and its corresponding value.
        """
        for key, value in added:
            self.attributes[key] = value

    def __str__(self) -> str:
        """
        Returns a string representation of the record.

        Returns:
            A tab-separated string of the annotation, sequence, and attributes.
        """
        return f"{str(self.annotation)}\t{str(self.sequence)}\t{str(self.attributes)}"
