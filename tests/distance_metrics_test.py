import pytest
from taggd.misc.distance_metrics import (
    py_hamming_distance,
    py_levenshtein_distance,
    py_subglobal_distance,
)


def test_hamming_distance_basic():
    assert py_hamming_distance("ACGT", "ACGA") == 1
    assert py_hamming_distance("ACGT", "AGGT") == 1
    assert py_hamming_distance("ACGT", "TTTT") == 3


def test_hamming_distance_exact_match():
    assert py_hamming_distance("ACGT", "ACGT") == 0


def test_hamming_distance_limit_exceeded():
    assert py_hamming_distance("ACGT", "TTTT", limit=2) == 3  # Exceeds limit
    assert py_hamming_distance("ACGT", "AGGT", limit=2) == 1  # Within limit


def test_hamming_distance_empty_sequences():
    assert py_hamming_distance("", "") == 0


def test_hamming_distance_unequal_length():
    with pytest.raises(ValueError):
        py_hamming_distance("ACGT", "ACG")


def test_levenshtein_distance_basic():
    assert py_levenshtein_distance("kitten", "sitting") == 3  # 2 subs, 1 insertion
    assert py_levenshtein_distance("flaw", "lawn") == 2  # 2 subs
    assert py_levenshtein_distance("abcdef", "azced") == 3  # 3 subs


def test_levenshtein_distance_exact_match():
    assert py_levenshtein_distance("ACGT", "ACGT") == 0


def test_levenshtein_distance_limit_exceeded():
    assert py_levenshtein_distance("kitten", "sitting", limit=2) == 3  # Exceeds limit
    assert py_levenshtein_distance("kitten", "kitten", limit=2) == 0  # Within limit


def test_levenshtein_distance_empty_sequences():
    assert py_levenshtein_distance("", "") == 0


def test_levenshtein_distance_one_empty_sequence():
    assert py_levenshtein_distance("kitten", "") == len("kitten")
    assert py_levenshtein_distance("", "sitting") == len("sitting")


def test_subglobal_distance_basic():
    assert py_subglobal_distance("ACGTACGT", "ACGT") == 0  # Perfect alignment
    assert py_subglobal_distance("ACGTACGT", "CGTT") == 1  # 1 substitution
    assert py_subglobal_distance("ACGTACGT", "TTTT") == 3  # 3 substitutions


def test_subglobal_distance_indels():
    # Insertions and deletions
    assert py_subglobal_distance("ACGTACGT", "ACGTA") == 0  # Perfect alignment
    assert py_subglobal_distance("ACGTACGT", "ACGCAC") == 1  # 1 insertion


def test_subglobal_distance_exact_match():
    assert py_subglobal_distance("ACGT", "ACGT") == 0


def test_subglobal_distance_s2_longer_than_s1():
    with pytest.raises(ValueError):
        py_subglobal_distance("ACGT", "ACGTAC")


def test_subglobal_distance_empty_sequences():
    with pytest.raises(ValueError):
        py_subglobal_distance("", "ACGT")
