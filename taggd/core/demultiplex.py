"""
Main demultiplexing module to demultiplex reads in parallel.
"""
from typing import Tuple, Any, List, Dict, Optional
import time
import asyncio
import aiofiles
import random
from tqdm import tqdm
from taggd.core.statistics import Statistics
from taggd.misc.kmer_utils import get_kmers_dicts  # type: ignore
from taggd.core.demultiplex_sub_functions import demultiplex_record  # type: ignore
from taggd.io.reads_reader_writer import ReadsReaderWriter
from taggd.core.match import get_match_header
import taggd.constants as constants


class DemultipleReads:
    """
    A class to process files (FASTA, FASTQ, SAM, BAM) in chunks using a Cython function
    and write the results to output files using a thread-safe queue.
    """

    def __init__(
        self,
        filename: str,
        true_barcodes: Dict[str, Any],
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
        trim_sequences: List[Tuple[int, int]],
        barcode_tag: str,
        subprocesses: int,
        output_matched: Optional[str] = None,
        output_ambiguous: Optional[str] = None,
        output_unmatched: Optional[str] = None,
        output_results: Optional[str] = None,
        chunk_size: int = 10000,
    ):
        """
        Initializes the DemultipleReads class.
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
        self.homopolymers = []
        if self.homopolymer_filter > 0:
            for c in "ACGT":
                self.homopolymers.append(c * homopolymer_filter)
        self.ambiguity_factor = ambiguity_factor
        self.no_offset_speedup = no_offset_speedup
        self.seed = seed
        self.multiple_hits_keep_one = multiple_hits_keep_one
        self.trim_sequences = trim_sequences
        self.barcode_tag = barcode_tag
        self.subprocesses = subprocesses
        self.output_matched = output_matched
        self.output_ambiguous = output_ambiguous
        self.output_unmatched = output_unmatched
        self.output_results = output_results
        self.chunk_size = chunk_size
        self.stop_signal = None  # Unique stop signal
        self.stats = Statistics(self.max_edit_distance)
        # Adjust the barcode length and the overhang if
        # we want to trim away helpers from the barcode
        self.barcode_length = len(list(true_barcodes.keys())[0])
        if trim_sequences is not None:
            for start, end in trim_sequences:
                self.barcode_length += end - start
        # Create k-mer mappings with ALL kmers
        self.kmer2seq = get_kmers_dicts(
            list(self.true_barcodes.keys()),
            self.k,
            round_robin=False,
            slider_increment=1,
        )
        print(f"Obtained {len(self.kmer2seq)} k-mers from the barcodes.")
        # define metric choice
        if self.metric == "Subglobal":
            self.metric_choice = constants.SUBGLOBAL
        elif self.metric == "Levenshtein":
            self.metric_choice = constants.LEVENSHTEIN
        else:
            self.metric_choice = constants.HAMMING
        # Create the reader/writer
        self.reader_writter = ReadsReaderWriter(self.filename)

    def _process_chunk(self, chunk: List[Any]) -> List[Tuple[Any, Any]]:
        """
        Processes a single chunk of records and returns the matches.

        Args:
            chunk: A chunk of file records.

        Returns:
            A list of tuples (match, record) to be added to the async queue.
        """
        results = []  # type: ignore
        for record in chunk:
            try:
                matches = demultiplex_record(
                    record,
                    self.barcode_tag,
                    self.true_barcodes,
                    self.homopolymers,
                    self.start_position,
                    self.barcode_length,
                    self.trim_sequences,
                    self.pre_overhang,
                    self.post_overhang,
                    self.k,
                    self.kmer2seq,
                    self.metric_choice,
                    self.max_edit_distance,
                    self.ambiguity_factor,
                    self.slider_increment,
                    self.no_offset_speedup,
                    self.multiple_hits_keep_one,
                )
                # Collect results for the queue
                results.extend((match, record) for match in matches)
            except Exception as e:
                print(f"Error processing record {record.annotation}: {e}")
                raise e
        return results

    async def _async_writer_results(self, queue: asyncio.Queue, file_path: str) -> None:
        """
        Asynchronously writes data (results) from a queue to a file (TSV).

        Args:
            queue: An asyncio.Queue containing data to write.
            file_path: Path to the output file.
        """
        async with aiofiles.open(file_path, mode="w") as writer:
            while True:
                item = await queue.get()
                if item is self.stop_signal:
                    break
                await writer.write(item)

    async def _async_writer_records(self, queue: asyncio.Queue, file_path: str) -> None:
        """
        Asynchronously writes data (records) from a queue to a file (BAM/SAM/FASTQ/FASTA).

        Args:
            queue: An asyncio.Queue containing data to write.
            file_path: Path to the output file.
        """
        with self.reader_writter.get_writer(file_path) as writer:
            while True:
                item = await queue.get()
                if item is self.stop_signal:
                    break
                writer.write(item)

    async def _async_writer_manager(self) -> None:
        """
        Creates and manages multiple async writer tasks for each output file.
        """
        # Prepare writers for each output category
        queues: Dict[str, asyncio.Queue] = {}
        tasks = []
        if self.output_results:
            queues["output_results"] = asyncio.Queue(maxsize=1000)
            tasks.append(
                asyncio.create_task(
                    self._async_writer_results(
                        queues["output_results"], self.output_results
                    )
                )
            )
            # Add a header to the output results
            await queues["output_results"].put(get_match_header() + "\n")
        if self.output_unmatched:
            queues["output_unmatched"] = asyncio.Queue(maxsize=1000)
            tasks.append(
                asyncio.create_task(
                    self._async_writer_records(
                        queues["output_unmatched"], self.output_unmatched
                    )
                )
            )
        if self.output_matched:
            queues["output_matched"] = asyncio.Queue(maxsize=1000)
            tasks.append(
                asyncio.create_task(
                    self._async_writer_records(
                        queues["output_matched"], self.output_matched
                    )
                )
            )
        if self.output_ambiguous:
            queues["output_ambiguous"] = asyncio.Queue(maxsize=1000)
            tasks.append(
                asyncio.create_task(
                    self._async_writer_records(
                        queues["output_ambiguous"], self.output_ambiguous
                    )
                )
            )

        # Writer loop
        try:
            while True:
                item = await self.async_write_queue.get()
                if item is self.stop_signal:
                    # for queue in queues.values():
                    #   await queue.put(self.stop_signal)
                    break

                match, record = item
                self.stats.total_reads += 1

                # Write match data to results file
                if "output_results" in queues:
                    await queues["output_results"].put(f"{str(match)}\n")

                # Write unmatched records
                if match.match_type == constants.UNMATCHED:
                    self.stats.unmatched += 1
                    if "output_unmatched" in queues:
                        await queues["output_unmatched"].put(record.unwrap())
                    continue

                # Append record with properties. B0:Z:Barcode, B1:Z:Prop1, B2:Z:prop3 ...
                tags = []
                bc = self.true_barcodes[match.barcode]
                # To avoid duplicated B0 tag when input is BAM/SAM we set instead of add
                tags.append(("B0:Z", match.barcode))
                for j in range(len(bc.attributes)):  # type: ignore
                    tags.append((f"B{j+1}:Z", bc.attributes[j]))  # type: ignore
                record.add_tags(tags)

                # Write matched records
                if match.match_type == constants.MATCHED_PERFECTLY:
                    self.stats.perfect_matches += 1
                    self.stats.edit_distance_counts[0] += 1
                    if "output_matched" in queues:
                        await queues["output_matched"].put(record.unwrap())

                if match.match_type == constants.MATCHED_UNAMBIGUOUSLY:
                    self.stats.imperfect_unambiguous_matches += 1
                    self.stats.edit_distance_counts[match.edit_distance] += 1
                    if "output_matched" in queues:
                        await queues["output_matched"].put(record.unwrap())

                # Write ambiguous records
                if match.match_type == constants.MATCHED_AMBIGUOUSLY:
                    self.stats.imperfect_ambiguous_matches += 1
                    if "output_ambiguous" in queues:
                        await queues["output_ambiguous"].put(record.unwrap())
        finally:
            # Ensure all queues receive the stop signal
            for queue in queues.values():
                await queue.put(self.stop_signal)
            await asyncio.gather(*tasks)

    def run(self) -> None:
        """
        Processes the input file in chunks and writes results asynchronously to output files.
        """
        start_time = time.time()
        random.seed(self.seed)
        print(f"Using {self.subprocesses} threads and {self.chunk_size} as chunk size.")

        # Queue for async writing
        self.async_write_queue: asyncio.Queue = asyncio.Queue(
            maxsize=self.chunk_size * self.subprocesses
        )

        async def main():
            writer_task = asyncio.create_task(self._async_writer_manager())
            processed_records = 0
            pbar = tqdm(desc="Processing records", unit="record", dynamic_ncols=True)

            try:
                current_chunk = []
                # Read the input file and divide it into chunks
                for record in self.reader_writter.reader_open():
                    current_chunk.append(record)
                    if len(current_chunk) == self.chunk_size:
                        # Process the chunk asynchronously in a thread
                        results = await asyncio.to_thread(
                            self._process_chunk, current_chunk
                        )
                        # Put the matches in the queue
                        for result in results:
                            await self.async_write_queue.put(result)
                        processed_records += len(current_chunk)
                        pbar.update(len(current_chunk))
                        current_chunk = []

                # Process any remaining records
                if current_chunk:
                    results = await asyncio.to_thread(
                        self._process_chunk, current_chunk
                    )
                    for result in results:
                        await self.async_write_queue.put(result)
                    processed_records += len(current_chunk)
                    pbar.update(len(current_chunk))

                pbar.close()
                print(f"Processed {processed_records} records.")

                # Send stop signal to queue
                await self.async_write_queue.put(self.stop_signal)
                await writer_task
            finally:
                # Close reader
                self.reader_writter.reader_close()

        asyncio.run(main())
        self.stats.total_reads_wr = (
            self.stats.perfect_matches
            + self.stats.imperfect_unambiguous_matches
            + self.stats.imperfect_ambiguous_matches
        )
        self.stats.time = time.time() - start_time
