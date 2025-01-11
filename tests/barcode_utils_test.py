import pytest
from taggd.io.barcode_utils import read_barcode_file, estimate_min_edit_distance


class Barcode:
    def __init__(self, seq, attributes):
        self.sequence = seq
        self.attributes = attributes


def test_read_barcode_file_valid(tmp_path):
    file_content = "ACGT attribute1 attribute2\nTGCA attribute3 attribute4\n"
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text(file_content)

    result = read_barcode_file(str(barcode_file))
    assert len(result) == 2
    assert "ACGT" in result
    assert result["ACGT"].sequence == "ACGT"
    assert result["ACGT"].attributes == ["attribute1", "attribute2"]


def test_read_barcode_file_varying_lengths(tmp_path):
    file_content = "ACGT attribute1 attribute2\nTGCAA attribute3 attribute4\n"
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text(file_content)

    with pytest.raises(ValueError, match="varying lengths barcodes"):
        read_barcode_file(str(barcode_file))


def test_read_barcode_file_duplicate_barcodes(tmp_path):
    file_content = "ACGT attribute1 attribute2\nACGT attribute3 attribute4\n"
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text(file_content)

    with pytest.raises(ValueError, match="duplicate barcode: ACGT"):
        read_barcode_file(str(barcode_file))


def test_read_barcode_file_empty(tmp_path):
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text("")

    result = read_barcode_file(str(barcode_file))
    assert result == {}


def test_estimate_min_edit_distance_basic():
    barcodes = {
        "ACGT": Barcode("ACGT", []),
        "TGCA": Barcode("TGCA", []),
        "AAAA": Barcode("AAAA", []),
    }

    result = estimate_min_edit_distance(barcodes, max_iters=10)
    assert result > 0  # There should be a positive minimum distance


def test_estimate_min_edit_distance_max_iters():
    barcodes = {
        "ACGT": Barcode("ACGT", []),
        "TGCA": Barcode("TGCA", []),
        "AAAA": Barcode("AAAA", []),
    }

    result = estimate_min_edit_distance(barcodes, max_iters=2)
    assert result == 3  # Stops early due to max_iters


def test_estimate_min_edit_distance_empty():
    barcodes = {}
    result = estimate_min_edit_distance(barcodes, max_iters=10)
    assert result == 1000000000  # No comparisons, return default minimum


def test_estimate_min_edit_distance_single_barcode():
    barcodes = {"ACGT": Barcode("ACGT", [])}
    result = estimate_min_edit_distance(barcodes, max_iters=10)
    assert result == 1000000000  # No comparisons possible
