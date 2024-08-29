# PorechopX

PorechopX is a customized and enhanced version of the [ARTICnetwork's fork of Porechop](https://github.com/artic-network/Porechop), a tool originally developed for finding and trimming adapters from Oxford
Nanopore reads. PorechopX introduces several key improvements to improve performance:

## Key Features and Modifications

- Replaced the argparse module with click for command-line parsing.
- Reconstructed the multiprocessing architecture to enable real-time writing of results to the output. This replaces the original behavior where results were only written after all reads were processed, improving efficiency and reducing memory usage.
- Switched from SeqAn to parasail for semi-global adapter alignment.
- Pure python implementation without manual compile of SeqAn library.

## Requirements

* Linux
* [Python](https://www.python.org/) >=3.10, <3.12

## Installation

Simple install from PyPI:

`pip install porechopx`

## Quick usage examples

__Basic adapter trimming:__<br>
`porechopx -i input_reads.fastq.gz -o output_reads.fastq.gz`

__Trimmed reads to stdout, if you prefer:__<br>
`porechopx -i input_reads.fastq.gz > output_reads.fastq`

__Demultiplex barcoded reads:__<br>
`porechopx -i input_reads.fastq.gz -b output_dir`

__Demultiplex barcoded reads, straight from Albacore output directory:__<br>
`porechopx -i albacore_dir -b output_dir`

__Also works with FASTA:__<br>
`porechopx -i input_reads.fasta -o output_reads.fasta`

__More verbose output:__<br>
`porechopx -i input_reads.fastq.gz -o output_reads.fastq.gz --verbosity 2`

__Got a big server?__<br>
`porechopx -i input_reads.fastq.gz -o output_reads.fastq.gz --threads 40`

## Customize adapters

- The ARTIC's version of Porechop allows user specific additional adapters

| Barcode name | Direction | 5' start barcode | 3' end barcode |
| - | - | - |
	使用 csv 格式指定自己的 barcode，格式为
第一列：Barcode name （必须包含“ Barcode “，前后带空格）
第二列：1
第三列：5’start barcode
第二列：3’end barcode

以截图为例：

自定义barcode文件为：
Custom Barcode 01,1,ACTTGTACTTCGTTCAGTTGCGTATTGCTTTAACGGTAGAGTTTGATCCTGGCTCAG,AAGTCGTAACAAGGTAACCGTAGTAACGTAAGCAATGCGTAA


## Usage

PorechopX provides the same command-line interface (CLI) as porechop. Just replace `porechop` with
`porechopx` for better performance!

```txt
Usage: porechopx [OPTIONS]

  PorechopX: a tool for finding adapters in Oxford Nanopore reads, trimming
  them from the ends and splitting reads with internal adapters

Main options:
  -i, --input TEXT                FASTA/FASTQ of input reads or a directory
                                  which will be recursively searched for FASTQ
                                  files  [required]
  -o, --output TEXT               Filename for FASTA or FASTQ of trimmed reads
                                  (if not set, trimmed reads will be printed
                                  to stdout)
  --barcode_stats_csv TEXT        Path to a csv file with start/ end/ middle
                                  barcode names and percentage identities for
                                  each given read ( if not set, no information
                                  will be printed)
  --format [auto|fasta|fastq|fasta.gz|fastq.gz]
                                  Output format for the reads - if auto, the
                                  format will be chosen based on the output
                                  filename or the input read format
  -v, --verbosity INTEGER         Level of progress information: 0 = none, 1 =
                                  some, 2 = lots, 3 = full - output will go to
                                  stdout if reads are saved to a file and
                                  stderr if reads are printed to stdout
  -t, --threads INTEGER           Number of threads to use for adapter
                                  alignment

Barcode binning settings:
  Control the binning of reads based on barcodes (i.e. barcode demultiplexing)

  -b, --barcode_dir TEXT          Reads will be binned based on their barcode
                                  and saved to separate files in this
                                  directory (incompatible with --output)
  --barcode_labels                Reads will have a label added to their
                                  header with their barcode
  --extended_labels               Reads will have an extended label added to
                                  their header with the barcode_call (if any),
                                  the best start/ end barcode hit and their
                                  identities, and whether a barcode is found
                                  in middle of read. (Dependent on
                                  --barcode_labels).
  --native_barcodes               Only attempts to match the 24 native
                                  barcodes
  --pcr_barcodes                  Only attempts to match the 96 PCR barcodes
  --rapid_barcodes                Only attempts to match the 12 rapid barcodes
  --limit_barcodes_to TEXT        Specify a list of barcodes to look for
                                  (numbers refer to native, PCR or rapid)
  --custom_barcodes TEXT          CSV file containing custom barcode sequences
  --barcode_threshold FLOAT       A read must have at least this percent
                                  identity to a barcode to be binned
  --barcode_diff FLOAT            If the difference between a read's best
                                  barcode identity and its second-best barcode
                                  identity is less than this value, it will
                                  not be put in a barcode bin (to exclude
                                  cases which are too close to call)
  --require_two_barcodes          Reads will only be put in barcode bins if
                                  they have a strong match for the barcode on
                                  both their start and end (default: a read
                                  can be binned with a match at its start or
                                  end)
  --untrimmed                     Bin reads but do not trim them (default:
                                  trim the reads)
  --discard_unassigned            Discard unassigned reads (instead of
                                  creating a "none" bin)

Adapter search settings:
  Control how the program determines which adapter sets are present

  --adapter_threshold FLOAT       An adapter set has to have at least this
                                  percent identity to be labelled as present
                                  and trimmed off (0 to 100)
  --check_reads INTEGER           This many reads will be aligned to all
                                  possible adapters to determine which adapter
                                  sets are present
  --scoring_scheme TEXT           Comma-delimited string of alignment scores:
                                  match, mismatch, gap open, gap extend

End adapter settings:
  Control the trimming of adapters from read ends

  --end_size INTEGER              The number of base pairs at each end of the
                                  read which will be searched for adapter
                                  sequences
  --min_trim_size INTEGER         Adapter alignments smaller than this will be
                                  ignored
  --extra_end_trim INTEGER        This many additional bases will be removed
                                  next to adapters found at the ends of reads
  --end_threshold FLOAT           Adapters at the ends of reads must have at
                                  least this percent identity to be removed (0
                                  to 100)

Middle adapter settings:
  Control the splitting of read from middle adapters

  --no_split                      Skip splitting reads based on middle
                                  adapters (default: split reads when an
                                  adapter is found in the middle)
  --discard_middle                Reads with middle adapters will be discarded
                                  (default: reads with middle adapters are
                                  split)
  --middle_threshold FLOAT        Adapters in the middle of reads must have at
                                  least this percent identity to be found (0
                                  to 100)
  --extra_middle_trim_good_side INTEGER
                                  This many additional bases will be removed
                                  next to middle adapters on their "good" side
  --extra_middle_trim_bad_side INTEGER
                                  This many additional bases will be removed
                                  next to middle adapters on their "bad" side
  --min_split_read_size INTEGER   Post-split read pieces smaller than this
                                  many base pairs will not be outputted

Help:
  --help                          Show this message and exit.
  --version                       Show the version and exit.
```

## Credits

PorechopX is based on the orginal version of [Porechop](https://github.com/rrwick/Porechop.git) and
the [modified version of Porechop](https://github.com/artic-network/Porechop) by ARTIC Network.
Many thanks for developing a convenient software for processing nanopore data.

## Documentation

For detailed description of the adapter trimming strategy, please refer to [Porechop Documentation](https://github.com/rrwick/Porechop/blob/master/README.md)

##  License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)