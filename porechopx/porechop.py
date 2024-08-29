"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class for a Nanopore read.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""
import sys
import edlib
import logging
import parasail
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
from porechopx.adapters import make_full_native_barcode_adapter, make_full_rapid_barcode_adapter
from porechopx.misc import yellow, red, add_line_breaks_to_sequence, END_FORMATTING, RED, YELLOW
logger = logging.getLogger("porechopx")


class NanoporeRead(object):

    def __init__(self, name, seq, quals):
        self.name = name

        self.seq = seq
        self.quals = quals
        if len(quals) < len(seq):
            self.quals += '+' * (len(seq) - len(quals))

        self.start_trim_amount = 0
        self.end_trim_amount = 0
        self.start_adapter_alignments = []
        self.end_adapter_alignments = []

        self.middle_adapter_positions = set()
        self.middle_trim_positions = set()
        self.middle_hit_str = ''
        self.middle_hit_id = []

        self.start_barcode_scores = {}
        self.end_barcode_scores = {}

        self.best_start_barcode = ('none', 0.0)
        self.best_end_barcode = ('none', 0.0)
        self.second_best_start_barcode = ('none', 0.0)
        self.second_best_end_barcode = ('none', 0.0)
        self.barcode_call = 'none'
        self.albacore_barcode_call = None

    def get_name(self, barcode_labels=False, extended_labels=False):
        if barcode_labels:
            if extended_labels:
                start_name, start_id, end_name, end_id, middle_name, middle_id = self.get_barcode_hits()
                bcode_info = " barcode={} start_name={} start_id={} end_name={} end_id={}".format(self.barcode_call, start_name, start_id, end_name, end_id)
                if middle_name:
                    return self.name + bcode_info + " middle_name=" + middle_name + " middle_id=" + middle_id
                else:
                    return self.name + bcode_info + " middle_name=none middle_id=none"
            else:
                return self.name + " barcode=" + self.barcode_call

        return self.name

    def get_barcode_hits(self):
        start_name, start_id = self.best_start_barcode
        end_name, end_id = self.best_end_barcode

        if self.middle_hit_id != []:
            best_middle_hit = sorted(self.middle_hit_id, key = lambda x : float(x[1]), reverse=True)[0]

            return start_name, start_id, end_name, end_id, best_middle_hit[0], best_middle_hit[1]

        return start_name, start_id, end_name, end_id, None, None

    def write_to_stats_csv(self, custom_output):
        """
        Write barcoding status into csv:
            start_time has been removed to support both dorado and guppy fastq
        """
        # TODO: unify read information for guppy and dorado basecalled fastq

        name = self.name.split()[0]
        start_name, start_id, end_name, end_id, middle_name, middle_id = self.get_barcode_hits()
        if not middle_name:
            middle_name, middle_id = ("none","none")
        custom_output.write(f"{name},{self.barcode_call},{start_name},{start_id},{end_name},{end_id},{middle_name},{middle_id}\n")

    def get_seq_with_start_end_adapters_trimmed(self):
        if not self.start_trim_amount and not self.end_trim_amount:
            return self.seq
        start_pos = self.start_trim_amount
        end_pos = len(self.seq) - self.end_trim_amount
        trimmed_seq = self.seq[start_pos:end_pos]
        return trimmed_seq

    def seq_length_with_start_end_adapters_trimmed(self):
        return len(self.get_seq_with_start_end_adapters_trimmed())

    def get_quals_with_start_end_adapters_trimmed(self):
        if not self.start_trim_amount and not self.end_trim_amount:
            return self.quals
        start_pos = self.start_trim_amount
        end_pos = len(self.quals) - self.end_trim_amount
        trimmed_quals = self.quals[start_pos:end_pos]
        return trimmed_quals

    def get_split_read_parts(self, min_split_read_size):
        """
        Returns the read split into parts as determined by the middle_trim_positions set.
        """
        trimmed_seq = self.get_seq_with_start_end_adapters_trimmed()
        trimmed_quals = self.get_quals_with_start_end_adapters_trimmed()
        split_read_parts = []
        part_seq, part_quals = [], []
        for i in range(len(trimmed_seq)):
            if i in self.middle_trim_positions:
                if part_seq:
                    split_read_parts.append((''.join(part_seq), ''.join(part_quals)))
                    part_seq, part_quals = [], []
            else:
                part_seq.append(trimmed_seq[i])
                part_quals.append(trimmed_quals[i])
        if part_seq:
            split_read_parts.append((''.join(part_seq), ''.join(part_quals)))
        split_read_parts = [x for x in split_read_parts if len(x[0]) >= min_split_read_size]
        return split_read_parts

    def get_fasta(self, min_split_read_size, discard_middle, untrimmed=False, barcode_labels=False, extended_labels=False):
        if not self.middle_trim_positions:
            if untrimmed:
                seq = self.seq
            else:
                seq = self.get_seq_with_start_end_adapters_trimmed()
            if not seq:  # Don't return empty sequences
                return ''
            return ''.join(['>', self.get_name(barcode_labels, extended_labels), '\n', add_line_breaks_to_sequence(seq, 70)])
        elif discard_middle:
            return ''
        else:
            fasta_str = ''
            for i, split_read_part in enumerate(self.get_split_read_parts(min_split_read_size)):
                read_name = add_number_to_read_name(self.get_name(), i + 1)
                if not split_read_part[0]:  # Don't return empty sequences
                    return ''
                seq = add_line_breaks_to_sequence(split_read_part[0], 70)
                fasta_str += ''.join(['>', read_name, '\n', seq])
            return fasta_str

    def get_fastq(self, min_split_read_size, discard_middle, untrimmed=False, barcode_labels=False, extended_labels=False):
        if not self.middle_trim_positions:
            if untrimmed:
                seq = self.seq
                quals = self.quals
            else:
                seq = self.get_seq_with_start_end_adapters_trimmed()
                quals = self.get_quals_with_start_end_adapters_trimmed()
            if not seq:  # Don't return empty sequences
                return ''
            return ''.join(['@', self.get_name(barcode_labels,extended_labels), '\n', seq, '\n+\n', quals, '\n'])
        elif discard_middle:
            return ''
        else:
            fastq_str = ''
            for i, split_read_part in enumerate(self.get_split_read_parts(min_split_read_size)):
                read_name = add_number_to_read_name(self.get_name(), i + 1)
                if not split_read_part[0]:  # Don't return empty sequences
                    return ''
                fastq_str += ''.join(['@', read_name, '\n', split_read_part[0], '\n+\n',
                                      split_read_part[1], '\n'])
            return fastq_str

    def find_start_trim(self, adapters, end_size, extra_trim_size, end_threshold,
                        scoring_scheme_vals, min_trim_size, check_barcodes, forward_or_reverse):
        """
        Aligns one or more adapter sequences and possibly adjusts the read's start trim amount based
        on the result.
        """
        read_seq_start = self.seq[:end_size]
        for adapter in adapters:
            full_score, partial_score, read_start, read_end = \
                align_adapter(read_seq_start, adapter.start_sequence[1], *scoring_scheme_vals)
            if partial_score > end_threshold and read_end != end_size and \
                    read_end - read_start >= min_trim_size:
                trim_amount = read_end + extra_trim_size
                self.start_trim_amount = max(self.start_trim_amount, trim_amount)
                self.start_adapter_alignments.append((adapter, full_score, partial_score,
                                                      read_start, read_end))
            if check_barcodes and adapter.is_barcode() and \
                    adapter.barcode_direction() == forward_or_reverse:
                self.start_barcode_scores[adapter.get_barcode_name()] = full_score

    def find_end_trim(self, adapters, end_size, extra_trim_size, end_threshold,
                      scoring_scheme_vals, min_trim_size, check_barcodes, forward_or_reverse):
        """
        Aligns one or more adapter sequences and possibly adjusts the read's end trim amount based
        on the result.
        """
        read_seq_end = self.seq[-end_size:]
        for adapter in adapters:
            if not adapter.end_sequence:
                continue
            full_score, partial_score, read_start, read_end = \
                align_adapter(read_seq_end, adapter.end_sequence[1], *scoring_scheme_vals)
            if partial_score > end_threshold and read_start != 0 and \
                    read_end - read_start >= min_trim_size:
                trim_amount = (end_size - read_start) + extra_trim_size
                self.end_trim_amount = max(self.end_trim_amount, trim_amount)
                self.end_adapter_alignments.append((adapter, full_score, partial_score,
                                                    read_start, read_end))
            if check_barcodes and adapter.is_barcode() and \
                    adapter.barcode_direction() == forward_or_reverse:
                self.end_barcode_scores[adapter.get_barcode_name()] = full_score

    def find_middle_adapters(self, adapters, middle_threshold, extra_middle_trim_good_side,
                             extra_middle_trim_bad_side, scoring_scheme_vals,
                             start_sequence_names, end_sequence_names):
        """
        Aligns an adapter sequence to the whole read to find places where the read should be split.
        """
        masked_seq = self.get_seq_with_start_end_adapters_trimmed()
        valid_base = len(masked_seq)
        for adapter_name, adapter_seq in adapters:

            # We keep aligning adapters as long we get strong hits, so we can find multiple
            # occurrences in a single read.
            while True:
                # Check if all bases are masked
                if valid_base == 0:
                    break

                full_score, _, read_start, read_end = align_adapter(
                    masked_seq, adapter_seq, *scoring_scheme_vals)
                if full_score >= middle_threshold:
                    masked_seq = masked_seq[:read_start] + '-' * (read_end - read_start) + \
                                 masked_seq[read_end:]
                    valid_base -= read_end - read_start
                    self.middle_adapter_positions.update(range(read_start, read_end))

                    self.middle_hit_str += '  ' + adapter_name + ' (read coords: ' + \
                                           str(read_start) + '-' + str(read_end) + ', ' + \
                                           'identity: ' + '%.1f' % full_score + '%)\n'
                    self.middle_hit_id.append((adapter_name, full_score))
                    trim_start = read_start - extra_middle_trim_good_side
                    if adapter_name in start_sequence_names:
                        trim_start = read_start - extra_middle_trim_bad_side

                    trim_end = read_end + extra_middle_trim_good_side
                    if adapter_name in end_sequence_names:
                        trim_end = read_end + extra_middle_trim_bad_side

                    self.middle_trim_positions.update(range(trim_start, trim_end))
                else:
                    break

    def formatted_start_seq(self, end_size, extra_trim_size):
        """
        Returns the start of the read sequence, with any found adapters highlighted in red.
        """
        start_seq = self.seq[:end_size]
        if not self.start_trim_amount:
            return start_seq
        red_bases = self.start_trim_amount - extra_trim_size
        formatted_str = ''
        if red_bases:
            formatted_str = red(start_seq[:red_bases])
        formatted_str += yellow(start_seq[red_bases:red_bases+extra_trim_size])
        formatted_str += start_seq[red_bases+extra_trim_size:]
        return formatted_str

    def formatted_end_seq(self, end_size, extra_trim_size):
        """
        Returns the end of the read sequence, with any found adapters highlighted in red.
        """
        end_seq = self.seq[-end_size:]
        if not self.end_trim_amount:
            return end_seq
        red_bases = self.end_trim_amount - extra_trim_size
        formatted_str = ''
        if red_bases:
            formatted_str = red(end_seq[-red_bases:])
        formatted_str = yellow(end_seq[-(red_bases+extra_trim_size):-red_bases]) + formatted_str
        formatted_str = end_seq[:-(red_bases+extra_trim_size)] + formatted_str
        return formatted_str

    def formatted_whole_seq(self, extra_trim_size):
        """
        Returns the entire read sequence, with any found adapters highlighted in red.
        """
        if not self.start_trim_amount and not self.end_trim_amount:
            return self.seq

        red_start_bases, red_end_bases = 0, 0
        if self.start_trim_amount:
            red_start_bases = self.start_trim_amount - extra_trim_size
        if self.end_trim_amount:
            red_end_bases = self.end_trim_amount - extra_trim_size
        if red_start_bases + red_end_bases >= len(self.seq):
            return red(self.seq)

        formatted_start, formatted_end = '', ''
        if self.start_trim_amount:
            formatted_start = red(self.seq[:red_start_bases])
        if self.end_trim_amount:
            formatted_end = red(self.seq[-red_end_bases:])
        middle = self.seq[red_start_bases:len(self.seq)-red_end_bases]

        if len(middle) <= extra_trim_size * 2:
            middle = yellow(middle)
        else:
            if self.start_trim_amount:
                middle = yellow(middle[:extra_trim_size]) + middle[extra_trim_size:]
            if self.end_trim_amount:
                middle = middle[:-extra_trim_size] + yellow(middle[-extra_trim_size:])
        return formatted_start + middle + formatted_end

    def formatted_start_and_end_seq(self, end_size, extra_trim_size, check_barcodes):
        read_seq = ''
        if check_barcodes:
            start_name, start_id = self.best_start_barcode
            end_name, end_id = self.best_end_barcode
            read_seq += 'start: ' + start_name + ' (' + '%.1f' % start_id + '%), '
            read_seq += 'end: ' + end_name + ' (' + '%.1f' % end_id + '%), '
            read_seq += 'barcode call: ' + self.barcode_call + '   '
        if len(self.seq) <= 2 * end_size:
            read_seq += self.formatted_whole_seq(extra_trim_size)
        else:
            read_seq += (self.formatted_start_seq(end_size, extra_trim_size) + '...' +
                         self.formatted_end_seq(end_size, extra_trim_size))
        return read_seq

    def full_start_end_output(self, end_size, extra_trim_size, check_barcodes):
        def get_alignment_string(aln):
            return aln[0].name + ', full score=' + str(aln[1]) + ', partial score=' + \
                   str(aln[2]) + ', read position: ' + str(aln[3]) + '-' + str(aln[4])
        output = self.get_name() + '\n'
        output += '  start: ' + self.formatted_start_seq(end_size, extra_trim_size) + '...\n'
        if self.start_adapter_alignments:
            output += '    start alignments:\n'
            for a in self.start_adapter_alignments:
                output += '      ' + get_alignment_string(a) + '\n'
        output += '  end:   ...' + self.formatted_end_seq(end_size, extra_trim_size) + '\n'
        if self.end_adapter_alignments:
            output += '    end alignments:\n'
            for a in self.end_adapter_alignments:
                output += '      ' + get_alignment_string(a) + '\n'
        if check_barcodes:
            start_name, start_id = self.best_start_barcode
            end_name, end_id = self.best_end_barcode
            output += '  Barcodes:\n'
            all_start_barcodes_str = ', '.join([b[0] + ' (' + '%.1f' % b[1] + '%)'
                                                for b in self.start_barcode_scores.items()])
            all_end_barcodes_str = ', '.join([b[0] + ' (' + '%.1f' % b[1] + '%)'
                                              for b in self.end_barcode_scores.items()])
            output += '    start barcodes:        ' + all_start_barcodes_str + '\n'
            output += '    end barcodes:          ' + all_end_barcodes_str + '\n'
            output += '    best start barcode:    ' + start_name + ' (' + '%.1f' % start_id + '%)\n'
            output += '    best end barcode:      ' + end_name + ' (' + '%.1f' % end_id + '%)\n'
            if self.albacore_barcode_call is not None:
                output += '    albacore barcode call: ' + self.albacore_barcode_call + '\n'
            output += '    final barcode call:    ' + self.barcode_call + '\n'
        return output

    def formatted_middle_seq(self):
        """
        If a middle adapter was found, this returns the relevant part of the read sequence, with
        the adapter highlighted in red.
        """
        if not self.middle_adapter_positions:
            return

        trimmed_seq = self.get_seq_with_start_end_adapters_trimmed()

        range_start = max(0, min(self.middle_trim_positions) - 100)
        range_end = min(len(trimmed_seq),
                        max(self.middle_trim_positions) + 100)
        formatted_str = '' if range_start == 0 else '(' + str(range_start) + ' bp)...'

        last_colour = None
        for i in range(range_start, range_end):
            char_colour = None
            if i in self.middle_trim_positions:
                char_colour = 'yellow'
            if i in self.middle_adapter_positions:
                char_colour = 'red'
            if char_colour != last_colour:
                formatted_str += END_FORMATTING
                if char_colour == 'yellow':
                    formatted_str += YELLOW
                if char_colour == 'red':
                    formatted_str += RED

            formatted_str += trimmed_seq[i]
            last_colour = char_colour
        if last_colour is not None:
            formatted_str += END_FORMATTING

        formatted_str += '' if range_end == len(trimmed_seq) \
            else '...(' + str(len(trimmed_seq) - range_end) + ' bp)'
        return formatted_str

    def middle_adapter_results(self, verbosity):
        if not self.middle_adapter_positions:
            return ''
        results = self.get_name() + '\n' + self.middle_hit_str
        if verbosity > 1:
            results += self.formatted_middle_seq() + '\n'
        return results

    def determine_barcode(self, barcode_threshold, barcode_diff, require_two_barcodes):
        """
        This function works through the logic of choosing a barcode for the read based on the
        settings and the read's barcode alignments. It stores its result in self.barcode_call.
        """
        start_barcode_scores = sorted(self.start_barcode_scores.items(), reverse=True,
                                      key=lambda x: x[1])
        end_barcode_scores = sorted(self.end_barcode_scores.items(), reverse=True,
                                    key=lambda x: x[1])

        if len(start_barcode_scores) >= 1:
            self.best_start_barcode = start_barcode_scores[0]
        if len(start_barcode_scores) >= 2:
            self.second_best_start_barcode = start_barcode_scores[1]
        if len(end_barcode_scores) >= 1:
            self.best_end_barcode = end_barcode_scores[0]
        if len(end_barcode_scores) >= 2:
            self.second_best_end_barcode = end_barcode_scores[1]

        try:
            # If the user set --require_two_barcodes, then the criteria are much more stringent.
            # Both the start and end barcodes need to be over the threshold, they both need to be
            # sufficiently better than their second-best barcode hit, and they need to match.
            if require_two_barcodes:
                start_over_threshold = (self.best_start_barcode[1] >= barcode_threshold)
                end_over_threshold = (self.best_end_barcode[1] >= barcode_threshold)
                start_good_diff = (self.best_start_barcode[1] >=
                                   self.second_best_start_barcode[1] + barcode_diff)
                end_good_diff = (self.best_end_barcode[1] >=
                                 self.second_best_end_barcode[1] + barcode_diff)
                start_end_match = (self.best_start_barcode[0] == self.best_end_barcode[0])

                assert (start_over_threshold and end_over_threshold and
                        start_good_diff and end_good_diff and start_end_match)
                self.barcode_call = self.best_start_barcode[0]

            # If the user didn't set --require_two_barcodes, then the criteria aren't so strict.
            # The start/end barcodes are analysed all together.
            else:
                # Combine the start and end barcodes into a single list (i.e. we no longer care
                # whether the hit was at the start or end of the read), only keeping the best score
                # for each barcode.
                all_barcode_scores = []
                included_barcodes = set()
                for name, score in sorted(start_barcode_scores + end_barcode_scores, reverse=True,
                                          key=lambda x: x[1]):
                    if name not in included_barcodes:
                        all_barcode_scores.append((name, score))
                        included_barcodes.add(name)

                if len(all_barcode_scores) >= 1:
                    best_overall_barcode = all_barcode_scores[0]
                else:
                    best_overall_barcode = ('none', 0.0)
                if len(all_barcode_scores) >= 2:
                    second_best_overall_barcode = all_barcode_scores[1]
                else:
                    second_best_overall_barcode = ('none', 0.0)

                over_threshold = (best_overall_barcode[1] >= barcode_threshold)
                good_diff = (best_overall_barcode[1] >=
                             second_best_overall_barcode[1] + barcode_diff)
                assert over_threshold
                assert good_diff

                self.barcode_call = best_overall_barcode[0]

        except AssertionError:
            self.barcode_call = 'none'

        # If the read has been binned by Albacore, then Porechop and Albacore must agree on the
        # barcode. If they don't, the read is unclassified.
        if self.albacore_barcode_call is not None and \
                self.barcode_call != self.albacore_barcode_call:
            self.barcode_call = 'none'


def align_adapter(read, adapter, match=3, mismatch=-6, gap_open=5, gap_extend=2):
    """
    SSW local alignment of (read, adapter)
    """
    score_size = 1
    matrix = parasail.matrix_create("ACGT", match, mismatch)
    profile = parasail.ssw_init(read, matrix, score_size)
    hit = parasail.ssw_profile(profile, adapter, gap_open, gap_extend)

    if hit is None or hit.read_end1 == hit.read_begin1: # If no hit
        read_start = -1
        read_end = 0
        aligned_region_percent_identity = 0.0
        full_adapter_percent_identity = 0.0
    else:
        read_start = hit.read_begin1
        read_end = hit.read_end1+1

        _read = read[hit.read_begin1:hit.read_end1+1]
        _adapter = adapter[hit.ref_begin1:hit.ref_end1+1]

        region_distance = edlib.align(_read, _adapter)['editDistance']
        aligned_region_percent_identity = 100 * (1 - region_distance / len(_adapter))

        # adapter_distance = edlib.align(_read, adapter)['editDistance']
        # full_adapter_percent_identity = 100 * (1 - adapter_distance / len(adapter))

        # Re-alignment of candidate adapter
        _read = read[hit.read_begin1-hit.ref_begin1:hit.read_end1+len(adapter)-hit.ref_end1]
        adapter_distance = edlib.align(_read, adapter)['editDistance']
        full_adapter_percent_identity = 100 * (1 - adapter_distance / len(adapter))

    return full_adapter_percent_identity, aligned_region_percent_identity, read_start, read_end


def align_adapter_global(read, adapter, match=3, mismatch=-6, gap_open=5, gap_extend=2):
    """
    Global alignment of (read, adapter)
    """
    matrix = parasail.matrix_create("ACGT", match, mismatch)
    hit = parasail.sg_qx_trace(read, adapter, gap_open, gap_extend, matrix)
    # hit = parasail.sg_trace(read, adapter, gap_open, gap_extend, matrix)
    if hit is None: # If no hit
        read_start = -1
        read_end = 0
        aligned_region_percent_identity = 0.0
        full_adapter_percent_identity = 0.0
    else:
        adapter_str = hit.traceback.ref
        read_start = len(adapter_str) - len(adapter_str.lstrip("-"))
        read_end = len(adapter_str) - len(adapter_str.rstrip("-"))

        # _read = read[hit.read_begin1:hit.read_end1+1]
        # _adapter = adapter[hit.ref_begin1:hit.ref_end1+1]

        n_match = sum([i == k for i, j, k in zip(
            hit.traceback.query, hit.traceback.comp, hit.traceback.ref)])
        aligned_region_percent_identity = 100 * n_match / len(adapter_str.strip("-"))
        full_adapter_percent_identity = 100 * n_match / len(adapter_str.strip("-"))

    return full_adapter_percent_identity, aligned_region_percent_identity, read_start, read_end


def align_adapter_set(read, adapter_set, end_size, scoring_scheme_vals):
    """
    This function aligns the adapter to the reads and updates the best score for the adapter.
    This is not to determine where to trim the reads, but rather to figure out which adapter
    sets are present in the data.
    """
    read_seq_start = read.seq[:end_size]
    score, _, _, _ = align_adapter(read_seq_start, adapter_set.start_sequence[1],
                                      *scoring_scheme_vals)
    adapter_set.best_start_score = max(adapter_set.best_start_score, score)
    if adapter_set.end_sequence:
        read_seq_end = read.seq[-end_size:]
        score, _, _, _ = align_adapter(read_seq_end, adapter_set.end_sequence[1],
                                          *scoring_scheme_vals)
        adapter_set.best_end_score = max(adapter_set.best_end_score, score)


def find_matching_adapter_sets(check_reads, verbosity, end_size, scoring_scheme_vals,
                               adapter_threshold, threads, search_adapters=None):
    """
    Aligns all of the adapter sets to the start/end of reads to see which (if any) matches best.
    """
    read_count = len(check_reads)
    logger.info('Looking for known adapter sets')

    search_adapter_count = len(search_adapters)

    # multi-thread is not useful here due to the GIL
    for read in tqdm(check_reads, ncols=75):
        for adapter_set in search_adapters:
            align_adapter_set(read, adapter_set, end_size, scoring_scheme_vals)
    sys.stderr.write('\n')

    return [x for x in search_adapters if x.best_start_or_end_score() >= adapter_threshold]


def exclude_end_adapters_for_rapid(matching_sets):
    """
    Rapid reads shouldn't have end adapters, so we don't want to look for them if this seems to be
    a rapid read set.
    """
    if 'Rapid adapter' in [x.name for x in matching_sets]:
        for s in matching_sets:
            s.end_sequence = None
    return matching_sets


def fix_up_1d2_sets(matching_sets):
    """
    The 1D^2 adapters share a lot of common sequence with the old SQK-MAP006_short adapters, so if
    there are strong 1D^2 hits, we can exclude the redundant SQK-MAP006_short adapters.
    """
    matching_set_names = [x.name for x in matching_sets]
    if '1D^2 part 1' in matching_set_names and '1D^2 part 2' in matching_set_names and \
            'SQK-MAP006 Short' in matching_set_names:
        part_1_score = [x for x in matching_sets
                        if x.name == '1D^2 part 1'][0].best_start_or_end_score()
        part_2_score = [x for x in matching_sets
                        if x.name == '1D^2 part 2'][0].best_start_or_end_score()
        sqk_score = [x for x in matching_sets
                     if x.name == 'SQK-MAP006 Short'][0].best_start_or_end_score()
        if part_1_score >= sqk_score and part_2_score >= sqk_score:
            matching_sets = [
                x for x in matching_sets if x.name != 'SQK-MAP006 Short']
    return matching_sets


def add_full_barcode_adapter_sets(matching_sets):
    """
    This function adds some new 'full' adapter sequences based on what was already found. For
    example, if the ligation adapters and the reverse barcode adapters are found, it assumes we are
    looking at a native barcoding run and so it adds the complete native barcoding adapter
    sequences (with the barcode's upstream and downstream context included).
    """
    matching_set_names = [x.name for x in matching_sets]

    for i in range(1, 97):

        # Native barcode full sequences
        if all(x in matching_set_names
               for x in ['SQK-NSK007', 'Barcode ' + str(i) + ' (reverse)']):
            matching_sets.append(make_full_native_barcode_adapter(i))

        # Rapid barcode full sequences
        if all(x in matching_set_names
               for x in ['SQK-NSK007', 'Rapid', 'Barcode ' + str(i) + ' (forward)']):
            matching_sets.append(make_full_rapid_barcode_adapter(i))

    return matching_sets


def display_adapter_set_results(search_adapters, matching_sets):
    from porechopx.misc import print_table
    table = [['Set', 'Best read start %ID', 'Best read end %ID']]
    row_colours = {}
    matching_set_names = [x.name for x in matching_sets]
    for adapter_set in search_adapters:
        start_score = '%.1f' % adapter_set.best_start_score
        end_score = '%.1f' % adapter_set.best_end_score
        table.append([adapter_set.name, start_score, end_score])
        if adapter_set.name in matching_set_names:
            row_colours[len(table) - 1] = 'green'
    print_table(table, sys.stderr, alignments='LRR', row_colour=row_colours,
                fixed_col_widths=[35, 8, 8])
    sys.stderr.write('\n')


def choose_barcoding_kit(adapter_sets):
    """
    If the user is sorting reads by barcode bin, choose one barcode configuration (rev comp
    barcodes at the start of the read or at the end of the read) and ignore the other.
    """
    forward_barcodes = 0
    reverse_barcodes = 0
    for adapter_set in adapter_sets:
        score = adapter_set.best_start_or_end_score()
        if 'Native' in adapter_set.name:
            reverse_barcodes += score
        elif 'Barcode' in adapter_set.name and '(forward)' in adapter_set.name:
            forward_barcodes += score
        elif 'Barcode' in adapter_set.name and '(reverse)' in adapter_set.name:
            reverse_barcodes += score

    if forward_barcodes > reverse_barcodes:
        logger.info('Barcodes determined to be in forward orientation')
        return 'forward'
    elif reverse_barcodes > forward_barcodes:
        logger.info('Barcodes determined to be in reverse orientation',)
        return 'reverse'
    else:
        return None


def add_number_to_read_name(read_name, number):
    if ' ' not in read_name:
        return read_name + '_' + str(number)
    else:
        return read_name.replace(' ', '_' + str(number) + ' ', 1)
