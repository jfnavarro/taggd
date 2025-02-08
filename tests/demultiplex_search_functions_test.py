import pytest
from taggd.core.demultiplex_search_functions import (
    py_get_candidates,
    py_get_distances,
    py_get_top_hits,
)
import taggd.constants as constants


def test_get_candidates_basic():
    kmer2seq = {"AAA": {"barcode1": [0, 3]}, "AAT": {"barcode2": [2]}}
    read_barcode = "AAAAT"
    result = py_get_candidates(read_barcode, 3, 1, kmer2seq, False, 0, 0, 2)
    assert result == ["barcode1", "barcode2"]


def test_get_candidates_mixed():
    kmer2seq = {"AAA": {"barcode1": [0, 3]}, "CCT": {"barcode2": [2]}}
    read_barcode = "AAAAT"
    result = py_get_candidates(read_barcode, 3, 1, kmer2seq, False, 0, 0, 2)
    assert result == ["barcode1"]


def test_get_candidates_no_matches():
    kmer2seq = {"CCC": {"barcode1": [0, 3]}}
    read_barcode = "AAAAT"
    result = py_get_candidates(read_barcode, 3, 1, kmer2seq, False, 0, 0, 2)
    assert result == []


def test_get_candidates_empty_inputs():
    result = py_get_candidates("", 3, 1, {}, False, 0, 0, 2)
    assert result == []


def test_get_candidates_with_overhang():
    kmer2seq = {"AAA": {"barcode1": [0]}, "AAT": {"barcode2": [2]}}
    read_barcode = "AAAAT"
    result = py_get_candidates(read_barcode, 3, 1, kmer2seq, False, 1, 1, 2)
    assert result == ["barcode1", "barcode2"]


def test_get_distances_basic():
    read_barcode = "ACGT"
    candidates = ["ACGA", "AGCT", "CCCC"]
    result = py_get_distances(read_barcode, candidates, constants.HAMMING, 2)
    assert result == [("ACGA", 1), ("AGCT", 2)]


def test_get_distances_with_levenshtein():
    read_barcode = "ACGT"
    candidates = ["ACG", "ACGTG", "TTTT"]
    result = py_get_distances(read_barcode, candidates, constants.LEVENSHTEIN, 2)
    expected = [
        ("ACG", 1),  # One deletion
        ("ACGTG", 1),  # One insertion
    ]
    assert result == expected


def test_get_distances_with_subglobal():
    read_barcode = "ACGT"
    candidates = ["ACGA", "AGCT", "TTTT"]
    result = py_get_distances(read_barcode, candidates, constants.SUBGLOBAL, 2)
    expected = [
        ("ACGA", 1),  # One mismatch
        ("AGCT", 2),  # Two mismatches
    ]
    assert result == expected


def test_get_distances_with_mixed_candidates():
    read_barcode = "ACGT"
    candidates = ["ACGT", "AGGT", "TTTT", "ACGTA"]
    result = py_get_distances(read_barcode, candidates, constants.LEVENSHTEIN, 1)
    expected = [
        ("ACGT", 0),  # Perfect match
        ("AGGT", 1),  # One substitution
        ("ACGTA", 1),  # One deletion
    ]
    assert result == expected


def test_get_distances_all_exceed():
    read_barcode = "ACGT"
    candidates = ["AAAA", "TTTT"]
    result = py_get_distances(read_barcode, candidates, constants.HAMMING, 1)
    assert result == []


def test_get_distances_empty_candidates():
    read_barcode = "ACGT"
    result = py_get_distances(read_barcode, [], constants.HAMMING, 2)
    assert result == []


def test_get_distances_invalid_metric_choice():
    read_barcode = "ACGT"
    candidates = ["AAAA", "TTTT"]
    with pytest.raises(ValueError):
        py_get_distances(read_barcode, candidates, 999, 2)


def test_get_distances_edge_cases():
    read_barcode = "ACGT"
    candidates = ["AAAA", "TTTT"]
    result = py_get_distances(read_barcode, candidates, constants.LEVENSHTEIN, 3)
    assert result == [("AAAA", 3), ("TTTT", 3)]


def test_get_top_hits_basic():
    qual_hits = [("barcode1", 1), ("barcode2", 2), ("barcode3", 3)]
    result = py_get_top_hits(qual_hits, 1.5)
    assert result == [("barcode1", 1), ("barcode2", 2)]


def test_get_top_hits_empty():
    result = py_get_top_hits([], 1.5)
    assert result is None


def test_get_top_hits_ambiguity_factor_one():
    qual_hits = [("barcode1", 1), ("barcode2", 2), ("barcode3", 3)]
    result = py_get_top_hits(qual_hits, 1.0)
    assert result == [("barcode1", 1)]


def test_get_top_hits_single_hit():
    qual_hits = [("barcode1", 1)]
    result = py_get_top_hits(qual_hits, 2.0)
    assert result == [("barcode1", 1)]


def test_get_top_hits_multi_hit():
    qual_hits = [("barcode1", 1), ("barcode2", 1)]
    result = py_get_top_hits(qual_hits, 2.0)
    assert result == [("barcode1", 1), ("barcode2", 1)]


def test_get_top_hits_high_ambiguity_factor():
    qual_hits = [("barcode1", 1), ("barcode2", 5), ("barcode3", 10)]
    result = py_get_top_hits(qual_hits, 10.0)
    assert result == qual_hits
