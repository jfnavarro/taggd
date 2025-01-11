import os
import tempfile
import filecmp
import time
import subprocess
import pytest


@pytest.fixture(scope="module")
def test_data():
    testdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    data = {
        "testdir": testdir,
        "inbarcodes": os.path.join(testdir, "data", "testbarcodes.tsv"),
        "insam": os.path.join(testdir, "data", "testset.sam"),
        "inbam": os.path.join(testdir, "data", "testset.bam"),
        "infq": os.path.join(testdir, "data", "testset.fq"),
        "infa": os.path.join(testdir, "data", "testset.fa"),
    }
    for key, filepath in data.items():
        if key != "testdir":
            assert os.path.exists(filepath), f"{filepath} does not exist"
    return data


def compare_to_expected_results(suffix, filename, file_description, outdir, testdir):
    filepath_from_test = os.path.join(outdir, filename)
    expected_result_dir = os.path.join(testdir, "expected_results", suffix)
    expected_result_filepath = os.path.join(expected_result_dir, filename)

    assert os.path.exists(
        expected_result_filepath
    ), f"{expected_result_filepath} does not exist"
    assert os.path.exists(filepath_from_test), f"{file_description} does not exist"

    comp = filecmp.cmp(filepath_from_test, expected_result_filepath, shallow=False)
    assert comp, f"{file_description} does not match the expected contents"


def validate_output_data(expName, suffix, outdir, testdir):
    print(f"# Validating test {expName}")

    # Verify existence of output files and temp files
    assert os.listdir(outdir), "Output folder is empty"

    time.sleep(2)
    compare_to_expected_results(
        suffix, f"outfile_matched.{suffix}", "Matched file", outdir, testdir
    )
    compare_to_expected_results(
        suffix, f"outfile_unmatched.{suffix}", "Unmatched file", outdir, testdir
    )
    compare_to_expected_results(
        suffix, f"outfile_ambiguous.{suffix}", "Ambiguous file", outdir, testdir
    )
    compare_to_expected_results(
        suffix, "outfile_results.tsv", "Results file", outdir, testdir
    )


def run_taggd_test(test_data, suffix, input_file, k, max_edit_distance, overhang):
    outdir = tempfile.mkdtemp(prefix=f"taggd_demultiplex_test_out_{suffix}_")
    args = [
        "taggd_demultiplex",
        "--k",
        str(k),
        "--max-edit-distance",
        str(max_edit_distance),
        "--overhang",
        str(overhang),
        "--subprocesses",
        "3",
        "--seed",
        "dsfiogwhgfsaeadsgfADSgsagaagd",
        test_data["inbarcodes"],
        input_file,
        os.path.join(outdir, "outfile"),
    ]

    # Run the demultiplexer
    try:
        print(f"\n# Running {suffix.upper()} test with parameters: {' '.join(args)}")
        subprocess.check_call(args)
        validate_output_data(
            f"Normal {suffix.upper()} test", suffix, outdir, test_data["testdir"]
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Running {suffix.upper()} test failed: {e}")


def test_normal_sam_run(test_data):
    run_taggd_test(
        test_data, "sam", test_data["insam"], k=7, max_edit_distance=7, overhang=2
    )


def test_normal_bam_run(test_data):
    run_taggd_test(
        test_data, "bam", test_data["inbam"], k=6, max_edit_distance=5, overhang=0
    )


def test_normal_fq_run(test_data):
    run_taggd_test(
        test_data, "fq", test_data["infq"], k=4, max_edit_distance=8, overhang=3
    )


def test_normal_fa_run(test_data):
    run_taggd_test(
        test_data, "fa", test_data["infa"], k=4, max_edit_distance=8, overhang=3
    )
