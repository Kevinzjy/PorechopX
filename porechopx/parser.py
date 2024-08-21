import os
import sys
import gzip
import logging
from tqdm import tqdm
import multiprocessing
from itertools import cycle
from typing import NamedTuple
from traceback import format_exc
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
        self.filepath = filepath
        self.chunk_size = chunk_size

        # Open handles
        self.file = open(self.filepath, 'rb')
        self.fh = gzip.open(self.file, 'rt')

        # Progress bar
        self.cursor = 0
        if verbose:
            self.pbar = tqdm(
                total=os.path.getsize(self.filepath),
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
        self.filepath = filepath
        self.chunk_size = chunk_size
        self.file = parse_fastx_file(filepath)

        # Progressbar
        if verbose:
            self.pbar = tqdm(
                total=os.path.getsize(self.filepath),
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
    def __init__(self, filepath, threads, chunk_size=100_000, cache_size=10) -> None:
        # Global parameters
        self.input = filepath
        self.threads = threads
        self.chunk_size = chunk_size

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
            # Exit process if exception is caught
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)

    def worker(self, queue_in, queue_out):
        try:
            while True:
                # Exit if error
                if self.exit.is_set():
                    os._exit(1)

                # Process data
                chunk = queue_in.get()
                if chunk is None:
                    break

                queue_out.put(chunk)

            # Finish
            queue_out.put(None)

        except Exception as e:
            # Exit process if exception is caught
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)

    def writer(self):
        try:
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
                    res = q.get()

                    # If finish
                    if res is None:
                        is_finish = True
                        break

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


class TrimFastQ(object):
    """
    Class for trimming nanopore reads
    """
    def __init__(self, infile) -> None:
        pass

    def data_loader(self):
        """Iter through query fastq.gz"""
        try:
            chunk_size = self.chunk_size
            exit_event = self.exit

            with open(self.in_file, 'rb') as f, gzip.open(f, 'rt') as g:
                chunk, cursor = [], 0
                for line in g:
                    if exit_event.is_set():
                        os._exit(1)

                    # Read fastq unit
                    header = line.rstrip().lstrip('@')
                    read_id = header.split()[0]
                    seq = g.readline().rstrip()
                    sep = g.readline().rstrip()
                    qual = g.readline().rstrip()
                    chunk.append([read_id, seq, qual])

                    # Send to queue
                    if len(chunk) >= chunk_size:
                        pos = f.tell()
                        self.q_fastx.put((pos-cursor, chunk))
                        chunk, cursor = [], pos

                # Send last chunk to queue
                if len(chunk) > 0:
                    self.q_fastx.put((f.tell()-cursor, chunk))

        except Exception as e:
            # Exit process if exception is caught
            pid = os.getpid()
            self.errors.append((pid, format_exc()))
            self.exit.set()
            os._exit(1)