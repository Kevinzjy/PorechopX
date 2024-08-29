import os
import sys
import gzip
import logging
import subprocess
import multiprocessing
import numpy as np
from tqdm import tqdm
from pathlib import Path
from itertools import cycle
from typing import NamedTuple
from traceback import format_exc
from collections import defaultdict
from porechopx.porechop import NanoporeRead
from porechopx.misc import int_to_str, print_table
logger = logging.getLogger("porechopx")

# Nanopore read data class
Read = NamedTuple("Read", [('id', str), ('seq', str), ('qual', str)])


class FastqReaderPy(object):
    """
    A class to iterate over records in a FASTQ file.

    Each record is returned as a NamedTuple: (id, seq, qual)
    """
    def __init__(self, filepath, chunk_size=1_000_000, verbose=True) -> None:
        """
        Initialize the FastqReader with the path to the FASTQ file.
        """
        self.input = filepath
        self.chunk_size = chunk_size

        # Open handles
        self.file = open(self.input, 'rb')
        self.fh = gzip.open(self.file, 'rt')

        # Progress bar
        self.cursor = 0
        if verbose:
            self.pbar = tqdm(
                total=os.path.getsize(self.input),
                ncols=75, unit='b', unit_scale=True, unit_divisor=1024
            )
        else:
            self.pbar = None

    def __iter__(self):
        """
        Returns the iterator object itself.
        """
        return self

    def __next__(self):
        """
        Return the next record from the FASTQ file.
        """
        chunk = []
        for line in self.fh:
            # Load next read
            header = line.rstrip().lstrip("@")
            sequence = self.fh.readline().rstrip()
            separator = self.fh.readline().rstrip()
            quality = self.fh.readline().rstrip()

            # Push into chunk
            chunk.append(Read(header, sequence, quality))
            if len(chunk) >= self.chunk_size:
                break

        # Update progress
        processed = self.file.tell() - self.cursor
        if self.pbar is not None:
            self.pbar.update(processed)
        self.cursor += processed

        # End of file
        if not chunk:
            self.close()
            raise StopIteration

        return chunk

    def close(self):
        """
        Close the FASTQ file if it's still open.
        """
        if not self.fh.closed:
            self.fh.close()

        if not self.file.closed:
            self.file.close()

        if self.pbar is not None:
            self.pbar.update(self.pbar.total - self.pbar.n)
            self.pbar.close()
            sys.stderr.flush()


class FastqReaderRS(object):
    """
    A class to iterate over records in a FASTQ file. but using needletail

    Each record is returned as an needtail instance
        object (_type_): _description_
    """
    def __init__(self, filepath, chunk_size=1_000_000, verbose=True) -> None:
        """
        Initialize the FastqReaderX with the path to the FASTQ file.
        """
        from needletail import parse_fastx_file
        self.input = filepath
        self.chunk_size = chunk_size
        self.file = parse_fastx_file(filepath)

        # Progressbar
        if verbose:
            self.pbar = tqdm(
                total=os.path.getsize(self.input),
                ncols=75, unit='b', unit_scale=True, unit_divisor=1024
            )
        else:
            self.pbar = None

    def __iter__(self):
        """
        Returns the iterator object itself.
        """
        return self

    def __next__(self):
        """
        Return the next record from the FASTQ file.
        """
        try:
            chunk, processed = [], 0
            for record in self.file:
                chunk.append(Read(record.id, record.seq, record.qual))
                processed += len(record.id) + len(record.seq)*2 + 6 # Byte to store a read
                if len(chunk) >= self.chunk_size:
                    break

            # Update progress
            # A rough estimation for guppy basecalled reads = bytes / 2
            if self.pbar is not None:
                self.pbar.update(min(processed/2, self.pbar.total - self.pbar.n))

            # End of file
            if not chunk:
                self.close()
                raise StopIteration

        except KeyboardInterrupt:
            # Catch KeyboardInterrupt for terminating needletail iterator
            # https://stackoverflow.com/questions/21120947/catching-keyboardinterrupt-in-python-during-program-shutdown
            logger.error('Interrupted')
            try:
                sys.exit(130)
            except SystemExit:
                os._exit(130)

        return chunk

    def close(self):
        """
        Close the FASTQ file if it's still open.
        """
        if self.pbar is not None:
            self.pbar.update(self.pbar.total - self.pbar.n)
            self.pbar.close()
            sys.stderr.flush()


# Only use FastqReaderRS when needletail is corrected installed
try:
    import needletail
    FastqReader = FastqReaderRS
except ImportError:
    FastqReader = FastqReaderPy


class FastqChoper(object):
    """
    Class for trmming nanopore reads
    """
    def __init__(self, filepath, matching_sets, verbosity, end_size, extra_trim_size,
                 end_threshold, scoring_scheme_vals, min_trim_size,
                 threads, check_barcodes, barcode_threshold, barcode_diff,
                 require_two_barcodes, forward_or_reverse_barcodes,
                 middle_threshold, extra_middle_trim_good_side, extra_middle_trim_bad_side,
                 discard_middle, discard_unassigned, no_split,
                 out_format, output, barcode_stats_csv,
                 min_split_read_size,
                 barcode_dir, barcode_labels, extended_labels, untrimmed,
                 chunk_size=100_000, cache_size=10) -> None:

        # Global parameters
        self.input = filepath
        self.matching_sets = matching_sets
        self.verbosity = verbosity

        # Parameters for end trimming
        self.end_size = end_size
        self.extra_trim_size = extra_trim_size
        self.end_threshold = end_threshold
        self.scoring_scheme_vals = scoring_scheme_vals
        self.min_trim_size = min_trim_size
        self.threads = threads
        self.check_barcodes = check_barcodes
        self.barcode_threshold = barcode_threshold
        self.barcode_diff = barcode_diff
        self.require_two_barcodes = require_two_barcodes
        self.forward_or_reverse_barcodes = forward_or_reverse_barcodes

        # Parameters for middle splitting
        self.middle_threshold = middle_threshold
        self.extra_middle_trim_good_side = extra_middle_trim_good_side
        self.extra_middle_trim_bad_side = extra_middle_trim_bad_side
        self.discard_middle = discard_middle
        self.discard_unassigned = discard_unassigned
        self.no_split = no_split

        # Parameters for output
        self.out_format = out_format
        self.output = Path(output) if output is not None else None
        self.barcode_stats_csv = barcode_stats_csv
        self.min_split_size = min_split_read_size
        self.barcode_dir = Path(barcode_dir) if barcode_dir is not None else None
        self.barcode_labels = barcode_labels
        self.extended_labels = extended_labels
        self.untrimmed = untrimmed

        self.chunk_size = chunk_size

        # Saving output files
        trimmed_or_untrimmed = 'untrimmed' if self.untrimmed else 'trimmed'
        if self.barcode_dir is not None:
            verb = 'Saving'
            destination = 'barcode-specific files'
        elif self.output is None:
            verb = 'Outputting'
            destination = 'stdout'
        else:
            verb = 'Saving'
            destination = 'file'
        logger.info(f'{verb} {trimmed_or_untrimmed} reads to {destination}')

        if subprocess.getstatusoutput('which pigz')[0] == 0:
            logger.info('pigz found - using it to compress instead of gzip')
        else:
            logger.info('pigz not found - using gzip to compress')

        # Init multiprocessing
        manager = multiprocessing.Manager()
        self.exit = multiprocessing.Event()
        self.errors = manager.list()
        self.pool = []

        # Init queues for dandling input and output data
        self.queues_in = [multiprocessing.Queue(maxsize=cache_size) for _ in range(self.threads)]
        self.queues_out = [multiprocessing.Queue(maxsize=cache_size) for _ in range(self.threads)]

        # Init Workers
        self.pool = []
        for q_in, q_out in zip(self.queues_in, self.queues_out):
            p = multiprocessing.Process(target=self.worker, args=(q_in, q_out))
            self.pool.append(p)
            p.start()

        # Start loading fastq
        self.process_read = multiprocessing.Process(target=self.loader)
        self.process_read.start()
        self.process_write = multiprocessing.Process(target=self.writer)
        self.process_write.start()

        # Wait
        self.join()

    def loader(self):
        """
        Process for reading fastq into chunk, and pass into queues_in
        """
        try:
            # Iterator
            # TODO: infer data format
            fq = FastqReader(self.input, self.chunk_size)

            # Pass data into queues
            for chunk, q in zip(fq, cycle(self.queues_in)):
                # Exit if error
                if self.exit.is_set():
                    os._exit(1)

                q.put(chunk)

            # End of file
            for q in self.queues_in:
                q.put(None)

        except Exception as e:
            for q in self.queues_in:
                q.put(None)

            # Exit process if exception is caught
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)

    def worker(self, queue_in, queue_out):
        try:
            adapters = []
            start_sequence_names = set()
            end_sequence_names = set()
            for matching_set in self.matching_sets:
                adapters.append(matching_set.start_sequence)
                start_sequence_names.add(matching_set.start_sequence[0])
                if matching_set.end_sequence:
                    end_sequence_names.add(matching_set.end_sequence[0])
                    if matching_set.end_sequence[1] == matching_set.start_sequence[1]:
                        continue
                    adapters.append(matching_set.end_sequence)

            while True:
                # Exit if error
                if self.exit.is_set():
                    os._exit(1)

                # Process data
                chunk = queue_in.get()
                if chunk is None:
                    break

                results = []
                for _read in chunk:
                    read = NanoporeRead(_read.id, _read.seq, _read.qual)

                    # find_adapters_at_read_ends
                    read.find_start_trim(
                        self.matching_sets, self.end_size, self.extra_trim_size, self.end_threshold,
                        self.scoring_scheme_vals, self.min_trim_size, self.check_barcodes,
                        self.forward_or_reverse_barcodes)
                    read.find_end_trim(
                        self.matching_sets, self.end_size, self.extra_trim_size, self.end_threshold,
                        self.scoring_scheme_vals, self.min_trim_size, self.check_barcodes,
                        self.forward_or_reverse_barcodes)

                    if self.check_barcodes:
                        read.determine_barcode(
                            self.barcode_threshold, self.barcode_diff, self.require_two_barcodes)

                    # elif verbosity == 2:
                    #     sys.stderr.write(read.formatted_start_and_end_seq(end_size, extra_trim_size, check_barcodes),
                    #         file=print_dest)
                    # elif verbosity > 2:
                    #     print(read.full_start_end_output(end_size, extra_trim_size, check_barcodes),
                    #         file=print_dest)

                    # find_adapters_in_read_middles
                    if not self.no_split:
                        if self.discard_unassigned and read.barcode_call == "none":
                            # if we are discarding unassigned barcodes then there is no point in checking for middle adapters
                            continue

                        read.find_middle_adapters(
                            adapters, self.middle_threshold, self.extra_middle_trim_good_side,
                            self.extra_middle_trim_bad_side, self.scoring_scheme_vals,
                            start_sequence_names, end_sequence_names)

                        # if read.middle_adapter_positions and verbosity > 1:
                            # print(read.middle_adapter_results(
                                # verbosity), file=print_dest, flush=True)

                    results.append(read)

                queue_out.put(results)

            # Finish
            queue_out.put(None)

        except Exception as e:
            # Exit process if exception is caught
            queue_out.put(None)
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)

    def writer(self):
        try:
            # Barcode summary
            if self.barcode_stats_csv is not None:
                custom_output = open(self.barcode_stats_csv, "w")
                custom_output.write("name,start_time,barcode_call,start_name,start_id,end_name,end_id,middle_name,middle_id\n")

            # TODO: auto set output format
            self.out_format = 'fastq.gz'
            # if out_format == 'auto':
            #     if output is None:
            #         out_format = read_type.lower()
            #         if barcode_dir is not None and input_filename.lower().endswith('.gz'):
            #             out_format += '.gz'
            #     elif '.fasta.gz' in output.lower():
            #         out_format = 'fasta.gz'
            #     elif '.fastq.gz' in output.lower():
            #         out_format = 'fastq.gz'
            #     elif '.fasta' in output.lower():
            #         out_format = 'fasta'
            #     elif '.fastq' in output.lower():
            #         out_format = 'fastq'
            #     else:
            #         out_format = read_type.lower()

            gzipped_out = False
            gzip_command = 'gzip'
            if self.out_format.endswith('.gz') and \
                (self.barcode_dir is not None or self.output is not None):
                gzipped_out = True
                self.out_format = self.out_format[:-3]
                if subprocess.getstatusoutput('which pigz')[0] == 0:
                    gzip_command = f'pigz -p {self.threads}'

            # Output reads to barcode bins
            if self.barcode_dir is not None:
                if not self.barcode_dir.is_dir():
                    self.barcode_dir.mkdir(parents=True, exist_ok=True)
                barcode_files = {}
                barcode_read_counts, barcode_base_counts = defaultdict(int), defaultdict(int)

            # Output all reads to stdout
            elif self.output is None:
                pass

            # Output all reads to file
            else:
                if gzipped_out:
                    out_filename = self.output.with_name(f'TEMP_{os.getpid()}.{self.output.stem}')
                else:
                    out_filename = self.output
                out_file = open(out_filename, 'wt')

            # Iter through results
            total_reads = 0
            start_trim_total, start_trim_count = 0, 0
            end_trim_total, end_trim_count = 0, 0
            middle_trim_count = 0

            is_finish = False
            while True:
                # Exit if exception
                if self.exit.is_set():
                    os._exit(1)

                # Exit if finish
                if is_finish:
                    break

                # Get output
                for q in cycle(self.queues_out):
                    chunk = q.get()

                    # If finish
                    if chunk is None:
                        is_finish = True
                        break

                    # Summarize trimmed bases
                    total_reads += len(chunk)

                    start_trim_amount = np.array([read.start_trim_amount for read in chunk])
                    start_trim_total += start_trim_amount.sum()
                    start_trim_count += (start_trim_amount > 0).sum()

                    end_trim_amount = np.array([read.end_trim_amount for read in chunk])
                    end_trim_total += end_trim_amount.sum()
                    end_trim_count += (end_trim_amount > 0).sum()

                    middle_trim_count += sum([1 if read.middle_adapter_positions else 0 for read in chunk])

                    # Save to files
                    for read in chunk:
                        # Barcode information
                        if self.barcode_stats_csv is not None:
                            read.write_to_stats_csv(custom_output)

                        # Drop unassigned reads if specified
                        if self.discard_unassigned and read.barcode_call == 'none':
                            continue

                        # Prepare output fasta/fastq str
                        if self.out_format == 'fasta':
                            read_str = read.get_fasta(
                                self.min_split_size, self.discard_middle, self.untrimmed,
                                self.barcode_labels, self.extended_labels)
                        else:
                            read_str = read.get_fastq(
                                self.min_split_size, self.discard_middle, self.untrimmed,
                                self.barcode_labels, self.extended_labels)
                        if not read_str:
                            continue

                        # Output reads to barcode bin
                        if self.barcode_dir is not None:
                            barcode_name = read.barcode_call
                            if barcode_name not in barcode_files:
                                barcode_files[barcode_name] = open(
                                    self.barcode_dir / f'{barcode_name}.{self.out_format}', 'wt'
                                )
                            barcode_files[barcode_name].write(read_str)
                            barcode_read_counts[barcode_name] += 1
                            if self.untrimmed:
                                seq_length = len(read.seq)
                            else:
                                seq_length = read.seq_length_with_start_end_adapters_trimmed()
                            barcode_base_counts[barcode_name] += seq_length

                        # Output all reads to stdin.
                        elif self.output is None:
                            print(read_str, end='')

                        # Output all reads to file.
                        else:
                            out_file.write(read_str)

            # Barcoding summary
            if self.barcode_stats_csv is not None:
                custom_output.close()
            if self.barcode_dir is not None:
                if gzipped_out:
                    logger.info("Compressing output")

                table = [['Barcode', 'Reads', 'Bases', 'File']]
                for barcode_name in sorted(barcode_files.keys()):
                    barcode_files[barcode_name].close()
                    bin_filename = self.barcode_dir / f'{barcode_name}.{self.out_format}'

                    if gzipped_out:
                        if not bin_filename.is_file():
                            continue
                        bin_filename_gz = bin_filename.with_name(bin_filename.name + ".gz")
                        if bin_filename_gz.is_file():
                            bin_filename_gz.unlink()

                        try:
                            subprocess.check_output(f'{gzip_command} {bin_filename}',
                                                    stderr=subprocess.STDOUT, shell=True)
                        except subprocess.CalledProcessError:
                            pass
                        bin_filename = bin_filename_gz

                    table_row = [barcode_name, int_to_str(barcode_read_counts[barcode_name]),
                                int_to_str(barcode_base_counts[barcode_name]), str(bin_filename)]
                    table.append(table_row)

                sys.stderr.write('\n')
                print_table(table, sys.stderr, alignments='LRRL', max_col_width=60, col_separation=2)
                sys.stderr.write('\n')

            elif self.output is None:
                pass

            else:
                out_file.close()
                if gzipped_out:
                    logger.info("Compressing output")
                    subprocess.check_output(f'{gzip_command} -c {out_filename} > {self.output}',
                                            stderr=subprocess.STDOUT, shell=True)
                    out_filename.unlink()
                logger.info(f'Saved result to {self.output}')

            logger.info("Finished trimming")

            # Trimming and splitting summary
            sys.stderr.write('\n' + int_to_str(start_trim_count).rjust(len(int_to_str(total_reads))) +
                ' / ' + int_to_str(total_reads) + ' reads had adapters trimmed from their start (' +
                int_to_str(start_trim_total) + ' bp removed)\n')
            sys.stderr.write(int_to_str(end_trim_count).rjust(len(int_to_str(total_reads))) + ' / ' +
                int_to_str(total_reads) + ' reads had adapters trimmed from their end (' +
                int_to_str(end_trim_total) + ' bp removed)\n')
            verb = 'discarded' if self.discard_middle else 'split'
            sys.stderr.write(int_to_str(middle_trim_count).rjust(len(int_to_str(total_reads))) +
                             ' / ' + int_to_str(total_reads) + ' reads were ' +
                             verb + ' based on middle adapters\n\n')

        except Exception as e:
            # Exit process if exception is caught
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)

    def join(self):
        """
        Wait for program finish
        """
        self.process_read.join()
        for p in self.pool:
            p.join()
        self.process_write.join()
        self.raise_exec()

    def raise_exec(self):
        """
        Raise exception if error
        """
        if self.exit.is_set():
            for pid, e in self.errors:
                logging.error(f"worker {pid}:\n {e}")
            logging.error("Terminated")
            sys.exit(1)
