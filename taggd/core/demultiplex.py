"""
"""
import os
from typing import Tuple, Any, List, Dict, Optional
import time
from concurrent.futures import ThreadPoolExecutor
import queue
import threading
from taggd.core.statistics import Statistics
from taggd.misc.kmer_utils import get_kmers_dicts
from taggd.core.demultiplex_sub_functions import demultiplex_record

# Metric constants
KILL = -1
SUBGLOBAL = 0
LEVENSHTEIN = 1
HAMMING = 2


class DemultipleReads:
    """
    A class to process files (FASTA, FASTQ, SAM, BAM) in chunks using a Cython function
    and write the results to output files using a thread-safe queue.
    """

    def __init__(
        self,
        filename: str,
        true_barcodes: Dict[str, Tuple[str, str]],
        k: int,
        metric: str,
        slider_increment: int,
        start_position: int,
        pre_overhang: int,
        post_overhang: int,
        max_edit_distance: int,
        homopolymer_filter: int,
        ambiguity_factor: float,
        no_offset_speedup: bool,
        seed: int,
        multiple_hits_keep_one: bool,
        trim_sequences: List[Tuple[str, str]],
        barcode_tag: str,
        subprocesses: int,
        output_reads: str,
        output_matched: Optional[str] = None,
        output_ambiguous: Optional[str] = None,
        output_unmatched: Optional[str] = None,
        output_results: Optional[str] = None,
        chunk_size: int = 100,
    ):
        """
        Initializes the FileProcessor.
        """
        self.filename = filename
        self.true_barcodes = true_barcodes
        self.k = k
        self.metric = metric
        self.slider_increment = slider_increment
        self.start_position = start_position
        self.pre_overhang = pre_overhang
        self.post_overhang = post_overhang
        self.max_edit_distance = max_edit_distance
        self.homopolymer_filter = homopolymer_filter
        self.ambiguity_factor = ambiguity_factor
        self.no_offset_speedup = no_offset_speedup
        self.seed = seed
        self.multiple_hits_keep_one = multiple_hits_keep_one
        self.trim_sequences = trim_sequences
        self.barcode_tag = barcode_tag
        self.subprocesses = subprocesses
        self.output_reads = output_reads
        self.output_matched = output_matched
        self.output_ambiguous = output_ambiguous
        self.output_unmatched = output_unmatched
        self.output_results = output_results
        self.chunk_size = chunk_size
        self.write_queue = queue.Queue()  # Thread-safe queue for writing
        self.stop_signal = None  # Signal to stop the writer thread
        self.stats = Statistics()
        # Adjust the barcode length and the overhang if
        # we want to trim away helpers from the barcode
        self.barcode_length = len(list(true_barcodes.keys())[0])
        if trim_sequences is not None:
            for start, end in trim_sequences:
                self.barcode_length += end - start
        # Create k-mer mappings with ALL kmers
        self.kmer2seq = get_kmers_dicts(list(self.true_barcodes_.keys()), self.k, False)
        # define metric choice
        if self.metric == "Subglobal":
            self.metric_choice = SUBGLOBAL
        elif self.metric == "Levenshtein":
            self.metric_choice = LEVENSHTEIN
        else:
            self.metric_choice = HAMMING

    def _process_chunk(self, chunk: List[Any]) -> None:
        """
        Processes a single chunk using the Cython function and adds results to the queue.

        Args:
            chunk: A chunk of file records.
        """
        for record in chunk:
            matches = demultiplex_record(record, **self.options)
            # Iterate over all matches (only more than one if ambiguous)
            for match in matches:
                self.write_queue.put((match, record))

    def _writer_thread(self) -> None:
        """
        A thread that writes results from the queue to appropriate files.
        """
        while True:
            item = self.write_queue.get()
            if item == self.stop_signal:  # Stop signal
                break
            match, record = item
            self.stats.total_reads += 1
            # Write to results file.
            if self.output_results is not None:
                self.writers["output_results"].write(f"{match}\n")
            # No match.
            if match.match_type == match.UNMATCHED:
                if self.output_unmatched is not None:
                    self.writers["output_unmatched"].write(record)
                    self.stats.total_reads_wr += 1
                self.stats.unmatched += 1
                continue
            # Append record with properties. B0:Z:Barcode, B1:Z:Prop1, B2:Z:prop3 ...
            bc = self.true_barcodes[match.barcode]
            tags = list()
            # To avoid duplicated B0 tag when input is BAM/SAM we set instead of add
            tags.append(("B0:Z", match.barcode))
            for j in range(len(bc.attributes)):
                tags.append((f"B{j+1}:Z", bc.attributes[j]))
            record.add_tags(tags)
            # Write to output file.
            if match.match_type == match.MATCHED_PERFECTLY:
                if self.output_reads is not None:
                    self.writers["output_reads"].write(record)
                    self.stats.total_reads_wr += 1
                self.stats.perfect_matches += 1
                self.stats.edit_distance_counts[0] += 1
            elif match.match_type == match.MATCHED_UNAMBIGUOUSLY:
                if self.output_matched is not None:
                    self.writers["output_matched"].write(record)
                self.stats.imperfect_unambiguous_matches += 1
                self.stats.edit_distance_counts[match.edit_distance] += 1
            elif match.match_type == match.MATCHED_AMBIGUOUSLY:
                if self.output_ambiguous is not None:
                    self.writers["output_ambiguous"].write(record)
                self.stats.imperfect_ambiguous_matches += 1
            else:
                raise ValueError(f"Invalid match type {match.match_type}")

    def run(self) -> None:
        """
        Processes the input file in chunks and writes results to different files based on the output category.

        Args:
            output_dir (str): Directory where output files will be written.
        """
        start_time = time.time()

        # Prepare writers for each output category
        writers = {
            "category1": dnaio.open(
                os.path.join(output_dir, "category1.fasta"), mode="w"
            ),
            "category2": dnaio.open(
                os.path.join(output_dir, "category2.fasta"), mode="w"
            ),
            # Add more categories as needed
        }

        writer_thread = threading.Thread(target=self._writer_thread, args=(writers,))
        writer_thread.start()

        try:
            with ThreadPoolExecutor(max_workers=self.num_processes) as executor:
                futures = []
                current_chunk = []
                reader = self._get_reader()

                # Read the file and divide it into chunks
                for record in reader:
                    current_chunk.append(record)
                    if len(current_chunk) == self.chunk_size:
                        futures.append(
                            executor.submit(self._process_chunk, current_chunk)
                        )
                        current_chunk = []

                # Process any remaining records
                if current_chunk:
                    futures.append(executor.submit(self._process_chunk, current_chunk))

                # Ensure all tasks complete
                for future in futures:
                    future.result()

            # Send stop signal to writer thread
            self.write_queue.put(self.stop_signal)
            writer_thread.join()
        finally:
            # Close all writers
            for writer in writers.values():
                writer.close()
        self.stats.time = time.time() - start_time
