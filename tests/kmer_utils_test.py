from taggd.misc.kmer_utils import get_kmers_dicts, py_get_kmers


def test_get_kmers_dicts_basic():
    seqs = ["ACGTACGT", "TTTTACGT"]
    k = 3
    slider_increment = 1
    result = get_kmers_dicts(seqs, k, slider_increment=slider_increment)

    # Verify k-mers for first sequence
    assert "ACG" in result
    assert result["ACG"]["ACGTACGT"] == [0, 4]

    # Verify k-mers for second sequence
    assert "TTT" in result
    assert result["TTT"]["TTTTACGT"] == [0, 1]


def test_get_kmers_dicts_round_robin():
    seqs = ["ACGTAC"]
    k = 3
    slider_increment = 1
    result = get_kmers_dicts(
        seqs, k, round_robin=True, slider_increment=slider_increment
    )

    # Verify round-robin k-mers
    assert "ACG" in result
    assert result["ACG"]["ACGTAC"] == [0]
    assert "TAC" in result
    assert result["TAC"]["ACGTAC"] == [3]


def test_get_kmers_dicts_with_slider_increment():
    seqs = ["ACGTACGT"]
    k = 3
    slider_increment = 2
    result = get_kmers_dicts(seqs, k, slider_increment=slider_increment)

    # Verify k-mers with slider increment
    assert "ACG" in result
    assert result["ACG"]["ACGTACGT"] == [0, 4]
    assert "GTA" in result
    assert result["GTA"]["ACGTACGT"] == [2]


def test_get_kmers_dicts_empty_sequence():
    seqs = [""]
    k = 3
    result = get_kmers_dicts(seqs, k)
    assert result == {}


def test_get_kmers_dicts_k_larger_than_sequence():
    seqs = ["AC"]
    k = 3
    result = get_kmers_dicts(seqs, k)
    assert result == {}


def test_get_kmers_basic():
    seq = "ACGTACGT"
    k = 3
    slider_increment = 1
    result = py_get_kmers(seq, k, slider_increment=slider_increment)

    # Verify k-mers
    assert result == [
        ("ACG", 0),
        ("CGT", 1),
        ("GTA", 2),
        ("TAC", 3),
        ("ACG", 4),
        ("CGT", 5),
    ]


def test_get_kmers_round_robin():
    seq = "ACGTAC"
    k = 3
    slider_increment = 1
    result = py_get_kmers(seq, k, round_robin=True, slider_increment=slider_increment)

    # Verify round-robin k-mers
    assert result == [
        ("ACG", 0),
        ("CGT", 1),
        ("GTA", 2),
        ("TAC", 3),
        ("ACA", 4),
        ("CAC", 5),
    ]


def test_get_kmers_with_slider_increment():
    seq = "ACGTACGT"
    k = 3
    slider_increment = 2
    result = py_get_kmers(seq, k, slider_increment=slider_increment)

    # Verify k-mers with slider increment
    assert result == [("ACG", 0), ("GTA", 2), ("ACG", 4)]


def test_get_kmers_empty_sequence():
    seq = ""
    k = 3
    slider_increment = 1
    result = py_get_kmers(seq, k, slider_increment=slider_increment)
    assert result == []


def test_get_kmers_k_larger_than_sequence():
    seq = "AC"
    k = 3
    slider_increment = 1
    result = py_get_kmers(seq, k, slider_increment=slider_increment)
    assert result == []
