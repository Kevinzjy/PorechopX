import os
import sys
import gzip
import logging
from tqdm import tqdm
from traceback import format_exc
from typing import NamedTuple
logger = logging.getLogger("porechopx")

# Nanopore read data class
Read = NamedTuple("Read", [('id', str), ('seq', str), ('qual', str)])


class FastqReaderPy(object):
    """
    A class to iterate over records in a FASTQ file.

    Each record is returned as a NamedTuple: (id, seq, qual)
    """
    def __init__(self, filepath, chunk_size=4_000, n_reads=None) -> None:
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
        self.pbar = tqdm(total=os.path.getsize(self.filepath),
                         ncols=75, unit='b', unit_scale=True, unit_divisor=1024)

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

        self.pbar.update(self.pbar.total - self.pbar.n)
        self.pbar.close()
        sys.stderr.flush()


class FastqReaderRS(object):
    """
    A class to iterate over records in a FASTQ file. but using needletail

    Each record is returned as an needtail instance
        object (_type_): _description_
    """
    def __init__(self, filepath, chunk_size=4_000, n_reads=None) -> None:
        """
        Initialize the FastqReaderX with the path to the FASTQ file.
        """
        from needletail import parse_fastx_file
        self.filepath = filepath
        self.chunk_size = chunk_size
        self.file = parse_fastx_file(filepath)

        # Progressbar
        self.pbar = tqdm(total=os.path.getsize(self.filepath),
                         ncols=75, unit='b', unit_scale=True, unit_divisor=1024)
    def __iter__(self):
        """
        Returns the iterator object itself.
        """
        return self

    def __next__(self):
        """
        Return the next record from the FASTQ file.
        """
        chunk, processed = [], 0
        for record in self.file:
            chunk.append(record)
            processed += len(record.id) + len(record.seq)*2 + 6 # Byte to store a read
            if len(chunk) >= self.chunk_size:
                break

        # Update progress
        # A rough estimation for guppy basecalled reads = bytes / 2
        self.pbar.update(min(processed/2, self.pbar.total - self.pbar.n))

        # End of file
        if not chunk:
            self.close()
            raise StopIteration

        return chunk

    def close(self):
        """
        Close the FASTQ file if it's still open.
        """
        self.pbar.update(self.pbar.total - self.pbar.n)
        self.pbar.close()
        sys.stderr.flush()


# Only use FastqReaderRS when needletail is corrected installed
try:
    import needletail
    FastqReader = FastqReaderRS
except ImportError:
    FastqReader = FastqReaderPy


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