#!/usr/bin/env python
"""
Runs taggd, a tool to demultiplex (link molecular barcodes back to) a file of genetic reads,
typically obtained by sequencing. For matched reads, the barcode and its related properties
are added to the read. Unmatched reads, ambiguously matched reads, and stats are by default
produced as output files as well.

The input ID file should be tab-delimited with the following format:
<barcode>   <prop1>   <prop2>   ...
<barcode>   <prop1>   <prop2>   ...

The input files can be in FASTA, FASTQ, SAM or BAM format. Matched files will be appended
with the barcode and properties like this:

B0:Z:<barcode> B1:Z:<prop1> B2:Z:<prop3> ...

Source:          https://github.com/SpatialTranscriptomicsResearch/taggd
Python package:  https://pypi.python.org/pypi/taggd
Contact:         joel.sjostrand@gmail.com;jc.fernandez.navarro@gmail.com
"""

import os
import time
import multiprocessing as mp
import argparse
from taggd.io.barcode_utils import read_barcode_file, estimate_min_edit_distance  # type: ignore
from taggd.core.demultiplex import DemultipleReads


def parse_arguments(argv=None):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Demultiplex reads using barcodes.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument(
        "barcodes_infile",
        help="Path to the barcode file containing true barcodes and their properties.",
        metavar="BARCODES_INFILE",
    )
    parser.add_argument(
        "reads_infile",
        help="Path to the reads file (FASTA, FASTQ, SAM, or BAM format).",
        metavar="READS_INFILE",
    )
    parser.add_argument(
        "outfile_prefix",
        help="Prefix for output files.",
        metavar="OUTFILE_PREFIX",
    )

    # Optional arguments
    parser.add_argument(
        "--no-matched-output",
        help="Do not output matched reads.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--no-ambiguous-output",
        help="Do not output ambiguous reads.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--no-unmatched-output",
        help="Do not output unmatched reads.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--no-results-output",
        help="Do not output a tab-separated results file with stats on the reads.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--start-position",
        type=int,
        help="Start position for barcodes in reads (default: %(default)d).",
        default=0,
        metavar="INT",
    )
    parser.add_argument(
        "--k",
        type=int,
        help="K-mer length (default: %(default)d).",
        default=6,
        metavar="INT",
    )
    parser.add_argument(
        "--max-edit-distance",
        type=int,
        help="Maximum edit distance for allowing hits (default: %(default)d).",
        default=2,
        metavar="INT",
    )
    parser.add_argument(
        "--metric",
        help="Distance metric: Subglobal, Levenshtein, or Hamming (default: %(default)s).",
        default="Subglobal",
        metavar="STRING",
    )
    parser.add_argument(
        "--ambiguity-factor",
        type=float,
        help="Factor for determining ambiguous hits (default: %(default).1f).",
        default=1.0,
        metavar="FLOAT",
    )
    parser.add_argument(
        "--slider-increment",
        type=int,
        help="Space between k-mer searches, "
        "0 yields k-mer length (default: %(default)d).",
        default=0,
        metavar="INT",
    )
    parser.add_argument(
        "--overhang",
        type=int,
        help="Additional flanking bases around read barcodes "
        "to allow for insertions when matching (default: %(default)d).",
        default=2,
        metavar="INT",
    )
    parser.add_argument(
        "--seed",
        help="Random number generator seed for shuffling ambiguous hits (default: %(default)s).",
        default=None,
        metavar="STRING",
    )
    parser.add_argument(
        "--homopolymer-filter",
        type=int,
        help="Excludes reads where the barcode contains a homopolymer of the given length. "
        "0 means no filter (default: %(default)d).",
        default=0,
        metavar="INT",
    )
    parser.add_argument(
        "--subprocesses",
        type=int,
        help="Number of subprocesses to start (default: 0, yielding number of machine cores - 1).",
        default=0,
        metavar="INT",
    )
    parser.add_argument(
        "--estimate-min-edit-distance",
        type=int,
        help="If set, estimates the minimum edit distance among true barcodes "
        "by comparing the specified number of pairs. 0 means no estimation (default: %(default)d).",
        default=0,
        metavar="INT",
    )
    parser.add_argument(
        "--no-offset-speedup",
        help="Disables offset speedup routine, increasing runtime but potentially yielding more hits.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--multiple-hits-keep-one",
        help="When multiple k-mer hits are found for a record, "
        "keep one as unambiguous and the rest as ambiguous.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--trim-sequences",
        nargs="+",
        type=int,
        help="Trims specified ranges from the barcodes. Provide ranges as START END START END ... "
        "where START and END are integer positions (0-based).",
        default=None,
        metavar="INT",
    )
    parser.add_argument(
        "--barcode-tag",
        type=str,
        help="Uses the sequence in the specified tag for barcode demultiplexing. "
        "The tag must be a two-letter string and is applicable only for SAM/BAM input files.",
        default=None,
        metavar="STRING",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        help="The chunk size (number of reads) for parallel processing (default: %(default)d).",
        default=10000,
        metavar="INT",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.4.0",
    )

    return parser.parse_args(argv)


def validate_arguments(options):
    """
    Validate command-line arguments.

    Args:
        options (argparse.Namespace): Parsed command-line arguments.

    Raises:
        ValueError: If any argument fails validation.
    """
    # Validate barcode input file
    if not os.path.isfile(options.barcodes_infile):
        raise ValueError(f"Invalid barcodes input file path: {options.barcodes_infile}")

    # Validate reads input file
    if not os.path.isfile(options.reads_infile):
        raise ValueError(f"Invalid reads input file path: {options.reads_infile}")

    # Validate reads input file format
    valid_formats = [".fastq", ".fq", ".fasta", ".fa", ".sam", ".bam"]
    if not any(options.reads_infile.lower().endswith(ext) for ext in valid_formats):
        raise ValueError(
            "Invalid reads input file format. Must be one of: "
            f"{', '.join(valid_formats)}."
        )

    # Validate output file prefix
    if not options.outfile_prefix or not options.outfile_prefix.strip():
        raise ValueError("Output file prefix cannot be empty.")

    # Validate k-mer length
    if options.k <= 0:
        raise ValueError("K-mer length (--k) must be greater than 0.")

    # Validate max edit distance
    if options.max_edit_distance < 0:
        raise ValueError(
            "Max edit distance (--max-edit-distance) must be 0 or greater."
        )

    # Validate metric
    valid_metrics = ["Subglobal", "Levenshtein", "Hamming"]
    if options.metric not in valid_metrics:
        raise ValueError(
            f"Invalid metric (--metric). Must be one of: {', '.join(valid_metrics)}."
        )

    # Validate ambiguity factor
    if options.ambiguity_factor < 1.0:
        raise ValueError(
            "Ambiguity factor (--ambiguity-factor) must be 1.0 or greater."
        )

    # Validate slider increment
    if options.slider_increment < 0:
        raise ValueError("Slider increment (--slider-increment) must be 0 or greater.")

    # Set slider increment to k if it is 0
    if options.slider_increment == 0:
        options.slider_increment = options.k

    # Validate start position
    if options.start_position < 0:
        raise ValueError("Start position (--start-position) must be 0 or greater.")

    # Validate overhang
    if options.overhang < 0:
        raise ValueError("Overhang (--overhang) must be 0 or greater.")

    # Ensure overhang is 0 for Hamming metric
    if options.metric == "Hamming" and options.overhang > 0:
        raise ValueError(
            "Overhang (--overhang) must be 0 when using the Hamming metric."
        )

    # Validate subprocesses
    if options.subprocesses < 0:
        raise ValueError(
            "Number of subprocesses (--subprocesses) must be 0 or greater."
        )

    # Validate homopolymer filter
    if options.homopolymer_filter < 0:
        raise ValueError(
            "Homopolymer filter (--homopolymer-filter) must be 0 or greater."
        )

    # Validate trim sequences
    if options.trim_sequences:
        if len(options.trim_sequences) % 2 != 0:
            raise ValueError(
                "Invalid trim sequences (--trim-sequences). The number of positions "
                "must be even, specifying start and end pairs."
            )
        if min(options.trim_sequences) < 0:
            raise ValueError(
                "Invalid trim sequences (--trim-sequences). Positions must be non-negative."
            )

    # Validate barcode tag
    if options.barcode_tag:
        if len(options.barcode_tag) != 2:
            raise ValueError(
                f"Invalid barcode tag (--barcode-tag). Must be a two-letter string, "
                f"but got '{options.barcode_tag}'."
            )
        if not (
            options.reads_infile.lower().endswith(".sam")
            or options.reads_infile.lower().endswith(".bam")
        ):
            raise ValueError(
                "Barcode tag (--barcode-tag) is only valid for SAM or BAM formatted input files."
            )

    # Validate estimate min edit distance
    if options.estimate_min_edit_distance < 0:
        raise ValueError(
            "Estimate min edit distance (--estimate-min-edit-distance) must be 0 or greater."
        )

    # If no output files are selected, warn the user
    if (
        options.no_matched_output
        and options.no_ambiguous_output
        and options.no_unmatched_output
        and options.no_results_output
    ):
        raise ValueError(
            "All output files are disabled (--no-matched-output, --no-ambiguous-output, "
            "--no-unmatched-output, --no-results-output). At least one output must be enabled."
        )

    if options.chunk_size < 100:
        raise ValueError(
            "Chunk size (--chunk-size) must be 100 or greater to avoid excessive overhead."
        )


def main(argv=None):
    """
    Main application.
    Starts a timer, create parameter parsers, parsers parameters
    and run all the steps for the demultiplexing.
    """
    start_time = time.time()

    # Parse arguments and validate
    options = parse_arguments(argv)
    validate_arguments(options)

    # Read barcodes file
    true_barcodes = read_barcode_file(options.barcodes_infile)

    # Paths
    frmt = options.reads_infile.split(".")[-1]
    fn_bc = os.path.abspath(options.barcodes_infile)
    fn_reads = os.path.abspath(options.reads_infile)
    fn_prefix = os.path.abspath(options.outfile_prefix)
    fn_matched = None if options.no_matched_output else fn_prefix + "_matched." + frmt
    fn_ambig = None if options.no_ambiguous_output else fn_prefix + "_ambiguous." + frmt
    fn_unmatched = (
        None if options.no_unmatched_output else fn_prefix + "_unmatched." + frmt
    )
    fn_results = None if options.no_results_output else fn_prefix + "_results.tsv"

    # Subprocesses
    if options.subprocesses == 0:
        options.subprocesses = mp.cpu_count() - 1

    print(f"# Options: {str(options).split('Namespace')[-1]}")
    print(f"# Barcodes input file: {fn_bc}")
    print(f"# Reads input file: {fn_reads}")
    print(f"# Matched output file: {fn_matched}")
    print(f"# Ambiguous output file: {fn_ambig}")
    print(f"# Unmatched output file: {fn_unmatched}")
    print(f"# Results output file: {fn_results}")
    print(f"# Number of barcodes in input: {len(true_barcodes)}")
    lngth = len(list(true_barcodes.keys())[0])
    print(f"# Barcode length: {lngth}")
    print(
        f"# Barcode length when overhang added: "
        f"{lngth + min(options.start_position, options.overhang) + options.overhang}"
    )

    # Check barcodes file
    if options.estimate_min_edit_distance > 0:
        min_dist = estimate_min_edit_distance(
            true_barcodes, options.estimate_min_edit_distance
        )
        if min_dist <= options.max_edit_distance:
            raise ValueError(
                "Invalid max edit distance: exceeds or equal "
                "to estimated minimum edit distance among true barcodes."
            )
        print(
            f"# Estimate of minimum edit distance between true barcodes (may be less): {min_dist}"
        )
    else:
        print(
            "# Estimate of minimum edit distance between true barcodes (may be less): Not estimated"
        )

    # Make the input trim coordinates a list of tuples
    trim_sequences = None
    if options.trim_sequences is not None:
        trim_sequences = []
        for i in range(len(options.trim_sequences) - 1):
            if i % 2 == 0:
                trim_sequences.append(
                    (options.trim_sequences[i], options.trim_sequences[i + 1])
                )

    # Demultiplex
    print("# Starting demultiplexing...")
    demux = DemultipleReads(
        fn_reads,
        true_barcodes,
        options.k,
        options.metric,
        options.slider_increment,
        options.start_position,
        options.overhang,
        options.overhang,
        options.max_edit_distance,
        options.homopolymer_filter,
        options.ambiguity_factor,
        options.no_offset_speedup,
        options.seed,
        options.multiple_hits_keep_one,
        trim_sequences,
        options.barcode_tag,
        options.subprocesses,
        fn_matched,
        fn_ambig,
        fn_unmatched,
        fn_results,
        options.chunk_size,
    )
    demux.run()
    print("# ...finished demultiplexing")
    print("# Wall time in secs: " + str(time.time() - start_time))
    print(str(demux.stats))


if __name__ == "__main__":
    main()
