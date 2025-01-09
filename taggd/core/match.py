"""
This represent a barcode match with all
the necessary information that can be used to
write to a file later
"""
# Match type constants
KILL = -1
UNMATCHED = 0
MATCHED_PERFECTLY = 1
MATCHED_UNAMBIGUOUSLY = 2
MATCHED_AMBIGUOUSLY = 3


def match_type_to_str(match_type: int) -> str:
    """
    Converts a match type integer to its corresponding string representation.

    Args:
        match_type (int): The integer representation of the match type.

    Returns:
        str: The string representation of the match type.
    """
    if match_type == UNMATCHED:
        return "UNMATCHED"
    elif match_type == MATCHED_PERFECTLY:
        return "MATCHED_PERFECTLY"
    elif match_type == MATCHED_UNAMBIGUOUSLY:
        return "MATCHED_UNAMBIGUOUSLY"
    elif match_type == MATCHED_AMBIGUOUSLY:
        return "MATCHED_AMBIGUOUSLY"
    return "KILL"


def get_match_header() -> str:
    """
    Returns a header for match-to-string conversion.

    Returns:
        str: The header string.
    """
    return "#Annotation\tMatch_result\tBarcode\tEdit_distance"


class Match:
    """
    Represents a match, storing relevant details about the record, match type,
    barcode, and edit distance.
    """

    def __init__(self, match_type: int, barcode: str, edit_distance: int = 1) -> None:
        """
        Initializes a Match instance.

        Args:
            match_type: The type of the match (e.g., MATCHED_PERFECTLY, UNMATCHED).
            barcode: The barcode associated with the match.
            edit_distance: The edit distance of the match. Defaults to 1.
        """
        self.match_type = match_type
        self.barcode = barcode
        self.edit_distance = edit_distance

    def __str__(self) -> str:
        """
        Converts the Match instance to a string representation.

        Returns:
            The string representation of the match in tab-separated format.
        """
        return (
            f"{match_type_to_str(self.match_type)}\t"
            f"{self.barcode}\t"
            f"{self.edit_distance}"
        )
