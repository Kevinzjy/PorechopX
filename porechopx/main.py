import os
import sys
import click
import importlib.metadata
from datetime import datetime
from porechopx import porechop, adapters, misc
from porechopx.parser import FastqChoper, FastqReader
from porechopx.logger import get_logger
__version__ = importlib.metadata.version('porechopx')


@click.command(help='PorechopX: a tool for finding adapters in Oxford Nanopore reads, trimming '
                    'them from the ends and splitting reads with internal adapters',
               context_settings={'show_default': True})
@click.version_option(__version__)
# Main options
@click.option('-i', '--input', 'input_filepath', required=True,
              help='FASTA/FASTQ of input reads or a directory which will be '
                   'recursively searched for FASTQ files')
@click.option('-o', '--output', 'output_filepath',
              help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                   'trimmed reads will be printed to stdout)')
@click.option('--barcode_stats_csv', default=None,
              help='Path to a csv file with start/ end/ middle barcode names '
                   'and percentage identities for each given read ( if not set, '
                   'no information will be printed)')
@click.option('--format', 'out_format', type=click.Choice(['auto', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz']),
              default='auto',
              help='Output format for the reads - if auto, the '
                   'format will be chosen based on the output filename or the input '
                   'read format')
@click.option('-v', '--verbosity', type=int, default=1,
              help='Level of progress information: 0 = none, 1 = some, 2 = lots, '
                   '3 = full - output will go to stdout if reads are saved to '
                   'a file and stderr if reads are printed to stdout')
@click.option('-t', '--threads', type=int, default=lambda: min(os.cpu_count(), 16),
              help='Number of threads to use for adapter alignment')
@click.option('-c', '--chunk_size', type=int, default=10_000,
              help='Number of reads per chunk')
# Barcode binning settings:
#   Control the binning of reads based on barcodes (i.e. barcode demultiplexing)
@click.option('-b', '--barcode_dir',
              help='Reads will be binned based on their barcode and saved to '
                   'separate files in this directory (incompatible with --output)')
@click.option('--barcode_labels', is_flag=True,
              help='Reads will have a label added to their header with their barcode')
@click.option('--extended_labels', is_flag=True,
              help='Reads will have an extended label added to their header with the '
                   'barcode_call (if any), the best start/ end barcode hit and their '
                   'identities, and whether a barcode is found in middle of read. '
                   '(Dependent on --barcode_labels).')
@click.option('--native_barcodes', is_flag=True,
              help='Only attempts to match the 24 native barcodes')
@click.option('--pcr_barcodes', is_flag=True,
              help='Only attempts to match the 96 PCR barcodes')
@click.option('--rapid_barcodes', is_flag=True,
              help='Only attempts to match the 12 rapid barcodes')
@click.option('--limit_barcodes_to', multiple=True,
              help='Specify a list of barcodes to look for (numbers refer to native, PCR or rapid)')
@click.option('--custom_barcodes',
              help='CSV file containing custom barcode sequences')
@click.option('--barcode_threshold', type=float, default=75.0,
              help='A read must have at least this percent identity to a barcode to be binned')
@click.option('--barcode_diff', type=float, default=5.0,
              help="If the difference between a read's best barcode identity and its second-best "
                   "barcode identity is less than this value, it will not be put in a barcode bin "
                   "(to exclude cases which are too close to call)")
@click.option('--require_two_barcodes', is_flag=True,
              help='Reads will only be put in barcode bins if they have a strong match for the '
                   'barcode on both their start and end (default: a read can be binned with '
                   'a match at its start or end)')
@click.option('--untrimmed', is_flag=True,
              help='Bin reads but do not trim them (default: trim the reads)')
@click.option('--discard_unassigned', is_flag=True,
              help='Discard unassigned reads (instead of creating a "none" bin)')
# Adapter search settings:
#   Control how the program determines which adapter sets are present
@click.option('--adapter_threshold', type=float, default=90.0,
              help='An adapter set has to have at least this percent identity to be labelled as '
                   'present and trimmed off (0 to 100)')
@click.option('--check_reads', 'check_read_count', type=int, default=10000,
              help='This many reads will be aligned to all possible adapters to determine which '
                   'adapter sets are present')
@click.option('--scoring_scheme', type=str, default='3,-6,5,2',
              help='Comma-delimited string of alignment scores: '
                   'match, mismatch, gap open, gap extend')
# End adapter settings:
#   Control the trimming of adapters from read ends
@click.option('--end_size', type=int, default=150,
              help='The number of base pairs at each end of the read which will be searched for '
                   'adapter sequences')
@click.option('--min_trim_size', type=int, default=10,
              help='Adapter alignments smaller than this will be ignored')
@click.option('--extra_end_trim', type=int, default=2,
              help='This many additional bases will be removed next to adapters found at '
                   'the ends of reads')
@click.option('--end_threshold', type=float, default=75.0,
              help='Adapters at the ends of reads must have at least this percent identity '
                   'to be removed (0 to 100)')
# Middle adapter settings:
#   Control the splitting of read from middle adapters
@click.option('--no_split', is_flag=True,
              help='Skip splitting reads based on middle adapters (default: split reads when '
                   'an adapter is found in the middle)')
@click.option('--discard_middle', is_flag=True,
              help='Reads with middle adapters will be discarded (default: reads with '
                   'middle adapters are split)')
@click.option('--middle_threshold', type=float, default=90.0,
              help='Adapters in the middle of reads must have at least this percent identity '
                   'to be found (0 to 100)')
@click.option('--extra_middle_trim_good_side', type=int, default=10,
              help='This many additional bases will be removed next to middle adapters '
                   'on their "good" side')
@click.option('--extra_middle_trim_bad_side', type=int, default=100,
              help='This many additional bases will be removed next to middle adapters '
                   'on their "bad" side')
@click.option('--min_split_read_size', type=int, default=1000,
              help='Post-split read pieces smaller than this many base pairs '
                   'will not be outputted')
def main(input_filepath, output_filepath, barcode_stats_csv, out_format, verbosity, threads,
         chunk_size, barcode_dir, barcode_labels, extended_labels,
         native_barcodes, pcr_barcodes, rapid_barcodes, limit_barcodes_to, custom_barcodes,
         barcode_threshold, barcode_diff, require_two_barcodes, untrimmed, discard_unassigned,
         adapter_threshold, check_read_count, scoring_scheme,
         end_size, min_trim_size, extra_end_trim, end_threshold, no_split,
         discard_middle, middle_threshold, extra_middle_trim_good_side, extra_middle_trim_bad_side,
         min_split_read_size):

     # Checking arguments
     try:
          scoring_scheme_vals = [int(x) for x in scoring_scheme.split(',')]
     except ValueError:
          sys.exit('Error: incorrectly formatted scoring scheme')
     if len(scoring_scheme_vals) != 4:
          sys.exit('Error: incorrectly formatted scoring scheme')

     if barcode_dir and output_filepath:
          sys.exit('Error: only one of the following options may be used: --output, --barcode_dir')

     if untrimmed and not barcode_dir:
          sys.exit('Error: --untrimmed can only be used with --barcode_dir')

     if barcode_dir:
          discard_middle = True

     if threads < 1:
          sys.exit('Error: at least one thread required')

     # Program starts
     logger = get_logger("porechopx", debug=verbosity > 1)

     start_time = datetime.now()

     search_adapters = None
     matching_sets = None
     forward_or_reverse_barcodes = None

     fq = FastqReader(input_filepath, check_read_count, False)
     check_reads = fq.__next__()
     fq.close()

     # Check adapters
     if native_barcodes or pcr_barcodes or rapid_barcodes:
          if native_barcodes:
               logger.info('Using native barcodes')
               # construct a smaller set of search adapters with only the 24 barcodes to speed up
               # the initial step
               barcodes_set = adapters.NATIVE_BARCODES
               forward_or_reverse_barcodes = 'reverse'
               # barcode_adapters = []
               # for i in range(1, 25):
               #     barcode_adapters.append(make_full_native_barcode_adapter(i))
          elif pcr_barcodes:
               logger.info('Using PCR barcodes')
               barcodes_set = adapters.PCR_BARCODES
               forward_or_reverse_barcodes = 'forward'
          else:
               logger.info('Using rapid barcodes')
               barcodes_set = adapters.RAPID_BARCODES
               forward_or_reverse_barcodes = 'forward'

          if limit_barcodes_to:
               logger.info(f'Limiting barcodes to : {limit_barcodes_to}')

          logger.info(f'{len(barcodes_set)} barcodes in search set')

          matching_sets = []
          if limit_barcodes_to:
               for barcode_number in limit_barcodes_to:
                    if barcode_number < 1 or barcode_number > len(barcodes_set):
                         logger.error('Barcode number out of range of chosen set (1-24 for native, '
                                      '1-12 for rapid, 1-96 for PCR)')
                         sys.exit(1)
                    matching_sets.append(barcodes_set[barcode_number - 1])
          else:
               matching_sets.extend(barcodes_set)

          check_barcodes = True

     else:
          if limit_barcodes_to:
               logger.error('To limit search to specific barcodes, specify whether using native, '
                            'PCR, or rapid barcodes')
               sys.exit(1)

          if custom_barcodes:
               matching_sets = adapters.load_custom_barcodes(custom_barcodes)
               logger.info('Using custom barcodes')
               logger.info(f'{len(matching_sets)} barcodes in search set')

          if not matching_sets:
               if not search_adapters:
                    # just add all of the adapters
                    search_adapters = adapters.NATIVE_BARCODES
                    search_adapters.extend(adapters.PCR_BARCODES)
                    search_adapters.extend(adapters.RAPID_BARCODES)
                    search_adapters.extend(adapters.OTHER_ADAPTERS)

               matching_sets = porechop.find_matching_adapter_sets(
                    check_reads, verbosity, end_size,
                    scoring_scheme_vals, adapter_threshold,
                    threads, search_adapters
               )
               matching_sets = porechop.exclude_end_adapters_for_rapid(matching_sets)
               matching_sets = porechop.fix_up_1d2_sets(matching_sets)
               matching_sets = porechop.add_full_barcode_adapter_sets(matching_sets)
               porechop.display_adapter_set_results(search_adapters, matching_sets)

               if barcode_dir or barcode_labels:
                    forward_or_reverse_barcodes = porechop.choose_barcoding_kit(matching_sets)

          if matching_sets:
               check_barcodes = (barcode_dir is not None or barcode_labels is not False)

     # Perform adapter & barcode trimming
     if matching_sets:
          if no_split:
               logger.info('Trimming adapters from read ends\n')
          else:
               verb = 'discarding' if discard_middle else 'splitting'
               logger.info(f'Trimming adapters from read ends, and {verb} reads containing middle adapters\n')

          name_len = max(max(len(x.start_sequence[0]) for x in matching_sets),
                         max(len(x.end_sequence[0]) if x.end_sequence else 0 for x in matching_sets))
          for matching_set in matching_sets:
               sys.stderr.write('  ' + matching_set.start_sequence[0].rjust(name_len) + ': ' +
                    misc.red(matching_set.start_sequence[1]) + '\n')
               if matching_set.end_sequence:
                    sys.stderr.write('  ' + matching_set.end_sequence[0].rjust(name_len) + ': ' +
                         misc.red(matching_set.end_sequence[1]) + '\n')
          sys.stderr.write('\nthreads: ' + str(threads) + '\n\n')

          # Multiprocess trimming
          FastqChoper(
               input_filepath, matching_sets, verbosity, end_size,
               extra_end_trim, end_threshold,
               scoring_scheme_vals, min_trim_size,
               threads, check_barcodes, barcode_threshold,
               barcode_diff, require_two_barcodes,
               forward_or_reverse_barcodes, middle_threshold,
               extra_middle_trim_good_side, extra_middle_trim_bad_side,
               discard_middle, discard_unassigned, no_split,
               out_format, output_filepath, barcode_stats_csv, min_split_read_size,
               barcode_dir, barcode_labels, extended_labels, untrimmed,
               chunk_size, 10,
          )

     else:
          logger.info('No adapters found - output reads are unchanged from input reads\n')

     # Finish
     time = datetime.now() - start_time
     logger.info(f'Time taken:  {time}')


if __name__ == '__main__':
    main()
