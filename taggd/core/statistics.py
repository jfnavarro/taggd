"""
A wrapper class to store some statistics
"""


class Statistics(object):
    """
    Shorthand for a statistics object.
    """

    def __init__(self, max_edit_distance):
        """
        Constructor
        """
        self.time = 0
        self.total_reads = 0
        self.total_reads_wr = 0
        self.perfect_matches = 0
        self.imperfect_unambiguous_matches = 0
        self.imperfect_ambiguous_matches = 0  # Non-unique
        self.unmatched = 0
        self.edit_distance_counts = [0] * (max_edit_distance + 1)

    def __str__(self):
        """
        String representation
        """
        mat = self.perfect_matches + self.imperfect_unambiguous_matches
        amb = (
            self.total_reads
            - self.perfect_matches
            - self.imperfect_unambiguous_matches
            - self.unmatched
        )

        return (
            f"# Total execution time in secs: {self.time}\n"
            f"# Total reads: {self.total_reads}\n"
            f"# Total reads written: {self.total_reads_wr}\n"
            f"# Matches: {mat}   [{mat * 100.0 / self.total_reads:.2f}%]\n"
            f"#   - Perfect matches: {self.perfect_matches}   "
            f"[{self.perfect_matches * 100.0 / self.total_reads:.2f}%]\n"
            f"#   - Imperfect matches: {self.imperfect_unambiguous_matches}   "
            f"[{self.imperfect_unambiguous_matches * 100.0 / self.total_reads:.2f}%]\n"
            f"# Ambiguous matches: {amb}   [{amb * 100.0 / self.total_reads:.2f}%]\n"
            f"#   - Non-unique ambiguous matches: {self.imperfect_ambiguous_matches}\n"
            f"# Unmatched: {self.unmatched}   "
            f"[{self.unmatched * 100.0 / self.total_reads:.2f}%]\n"
            f"# Matched edit distance counts for 0,1,...: {self.edit_distance_counts}"
        )

    def __iadd__(self, other):
        """
        Overrided += operator
        """
        self.time += other.time
        self.total_reads += other.total_reads
        self.total_reads_wr += other.total_reads_wr
        self.perfect_matches += other.perfect_matches
        self.imperfect_unambiguous_matches += other.imperfect_unambiguous_matches
        self.imperfect_ambiguous_matches += other.imperfect_ambiguous_matches
        self.unmatched += other.unmatched
        self.edit_distance_counts = [
            x + y for x, y in zip(self.edit_distance_counts, other.edit_distance_counts)
        ]
        return self
