import pytest
from unittest.mock import MagicMock
from taggd.core.demultiplex_sub_functions import demultiplex_record, py_trim_helpers
from taggd.misc.kmer_utils import get_kmers_dicts
import taggd.constants as constants


def test_trim_helpers():
    seq = "AAACCCGGGTTT"
    trim_sequences = [(3, 6), (9, 12)]
    result = py_trim_helpers(seq, trim_sequences)
    assert result == "AAAGGG"


def test_trim_helpers_empty():
    seq = "AAACCCGGGTTT"
    result = py_trim_helpers(seq, [])
    assert result == seq


def test_demultiplex_record_perfect_hit():
    rec = MagicMock()
    rec.annotation = "test_annotation"
    rec.sequence = "A" * 20
    rec.attributes = {"tags": [("BC", "AAAAAA")]}  # Mock attributes
    k = 3
    true_barcodes = {"AAAAAA"}
    kmer2seq = get_kmers_dicts(list(true_barcodes), k, False)
    result = demultiplex_record(
        rec,
        barcode_tag="BC",
        true_barcodes=true_barcodes,
        homopolymers=[],
        start_position=0,
        barcode_length=6,
        trim_sequences=None,
        pre_overhang=0,
        post_overhang=0,
        k=k,
        kmer2seq=kmer2seq,
        metric_choice=constants.HAMMING,
        max_edit_distance=1,
        ambiguity_factor=0.5,
        slider_increment=1,
        no_offset_speedup=False,
        multiple_hits_keep_one=False,
    )
    assert len(result) == 1
    assert result[0].match_type == constants.MATCHED_PERFECTLY


def test_demultiplex_record_unmatched():
    rec = MagicMock()
    rec.annotation = "test_annotation"
    rec.sequence = "A" * 20
    rec.attributes = {"tags": [("BC", "CCCCCC")]}  # Mock attributes
    k = 3
    true_barcodes = {"AAAAAA"}
    kmer2seq = get_kmers_dicts(list(true_barcodes), k, False)
    result = demultiplex_record(
        rec,
        barcode_tag="BC",
        true_barcodes=true_barcodes,
        homopolymers=[],
        start_position=0,
        barcode_length=6,
        trim_sequences=None,
        pre_overhang=0,
        post_overhang=0,
        k=k,
        kmer2seq=kmer2seq,
        metric_choice=constants.HAMMING,
        max_edit_distance=1,
        ambiguity_factor=0.5,
        slider_increment=1,
        no_offset_speedup=False,
        multiple_hits_keep_one=False,
    )
    assert len(result) == 1
    assert result[0].match_type == constants.UNMATCHED


@pytest.mark.parametrize(
    "barcode, true_barcodes, pre_overhang, post_overhang, expected_match_type",
    [
        ("AAGGTT", {"AAGGCC", "AAGGTT"}, 2, 0, constants.MATCHED_UNAMBIGUOUSLY),
        ("AAGGTT", {"AAGGCC"}, 2, 0, constants.MATCHED_UNAMBIGUOUSLY),
        ("AAGGTT", {"AAGGTC", "AAGGTA"}, 0, 0, constants.MATCHED_AMBIGUOUSLY),
        ("AAGGTT", {"AAGGTC"}, 0, 0, constants.MATCHED_UNAMBIGUOUSLY),
    ],
)
def test_demultiplex_record_matches(
    barcode, true_barcodes, pre_overhang, post_overhang, expected_match_type
):
    rec = MagicMock()
    rec.annotation = "test_annotation"
    rec.sequence = "A" * 20
    rec.attributes = {"tags": [("BC", barcode)]}  # Mock attributes
    k = 3
    kmer2seq = get_kmers_dicts(list(true_barcodes), k, False)
    result = demultiplex_record(
        rec,
        barcode_tag="BC",
        true_barcodes=true_barcodes,
        homopolymers=[],
        start_position=2 if pre_overhang > 0 else 0,
        barcode_length=6,
        trim_sequences=None,
        pre_overhang=pre_overhang,
        post_overhang=post_overhang,
        k=k,
        kmer2seq=kmer2seq,
        metric_choice=constants.HAMMING,
        max_edit_distance=2,
        ambiguity_factor=1,
        slider_increment=1,
        no_offset_speedup=False,
        multiple_hits_keep_one=False,
    )
    assert len(result) >= 1
    assert all(r.match_type == expected_match_type for r in result)


@pytest.mark.parametrize(
    "barcode, true_barcodes, pre_overhang, post_overhang, expected_match_type",
    [
        ("AAGGTT", {"AAGGCC", "AAGGTT"}, 2, 0, constants.MATCHED_UNAMBIGUOUSLY),
        ("AAGGTT", {"AAGGCC"}, 2, 0, constants.MATCHED_UNAMBIGUOUSLY),
        ("AAGGTT", {"AAGGCC", "AAGGCA"}, 0, 0, constants.MATCHED_AMBIGUOUSLY),
        ("AAGGTT", {"AAGGCC"}, 0, 0, constants.MATCHED_UNAMBIGUOUSLY),
    ],
)
def test_demultiplex_record_k_equals_barcode_length_ambiguous_and_partial_matches(
    barcode, true_barcodes, pre_overhang, post_overhang, expected_match_type
):
    rec = MagicMock()
    rec.annotation = "test_annotation"
    rec.sequence = "A" * 20
    rec.attributes = {"tags": [("BC", barcode)]}  # Mock attributes
    k = 3
    kmer2seq = get_kmers_dicts(list(true_barcodes), k, False, k)
    print(kmer2seq)
    result = demultiplex_record(
        rec,
        barcode_tag="BC",
        true_barcodes=true_barcodes,
        homopolymers=[],
        start_position=2 if pre_overhang > 0 else 0,
        barcode_length=6,
        trim_sequences=None,
        pre_overhang=pre_overhang,
        post_overhang=post_overhang,
        k=k,
        kmer2seq=kmer2seq,
        metric_choice=1,
        max_edit_distance=2,
        ambiguity_factor=1,
        slider_increment=k,
        no_offset_speedup=False,
        multiple_hits_keep_one=False,
    )
    assert len(result) >= 1
    assert all(r.match_type == expected_match_type for r in result)
