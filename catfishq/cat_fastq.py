from __future__ import print_function

import argparse
import logging
import math
import os
import re
from typing import List, Set, Optional

import pysam
import sys
from pathlib import Path
from datetime import datetime, timedelta

from . import __version__



LOOKUP = []

for q in range(100):
    LOOKUP.append(pow(10, -0.1 * q))


channel_regex_pattern = re.compile("(^|\s)ch=(?P<channel>\d+)")
"""
simple pattern for integer `channel` id in fastq comment section
"""
channel_range_regex_pattern = re.compile("^(?P<c1>\d+)(-(?P<c2>\d+))?$")
"""
patter for integer entry (`c1`) or range (`c1-c2`) 
"""


def _compute_mean_qscore(scores):
    """Returns the phred score corresponding to the mean of the probabilities
    associated with the phred scores provided.

    :param scores: Iterable of phred scores.

    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if not scores:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += LOOKUP[val]
    mean_prob = sum_prob / len(scores)
    return -10.0 * math.log10(mean_prob)


def _parse_channels_input(channels_input: str) -> List[int]:
    match = channel_range_regex_pattern.search(channels_input.strip())
    if match is None:
        raise ValueError(f"Channels input '{channels_input}' does not specify a[-b] single[range] integer pattern")
    le = int(match.group('c1'))
    he = le + 1
    if match.group('c2') is not None:
        he = int(match.group('c2')) + 1
    if he <= le:
        raise ValueError(f"Channels input '{channels_input}' clopen range has higher end '{he}' <= than lower end '{le}'")
    return list(range(le, he))


def get_channels_set(channels_input_list: List[str]) -> Set[int]:
    result: Set[int] = set()
    for entry in channels_input_list:
        for channel in _parse_channels_input(channels_input=entry):
            result.add(channel)
    return result


def get_channel_from_comment(comment: str) -> Optional[int]:
    if not comment:
        return None
    match = channel_regex_pattern.search(comment)
    if match is not None:
        return int(match.group('channel'))
    return None


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Cat long lists of FASTQ files"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="OUT",
        type=str,
        default="/dev/stdout",
        help="Output file. (default: stdout)",
    )
    parser.add_argument(
        "-l", "--min-length", dest="MIN_LEN", type=int, default=0, help="Minimum read length"
    )
    parser.add_argument(
        "-q", "--min-qscore", dest="MIN_QSCORE", type=int, default=0, help="Minimum q-score"
    )
    parser.add_argument(
        "--max-sequencing-time", dest="MAX_SEQ_TIME", type=int, default=None, help="Only output reads that where sequenced at or up to the given time (minutes)."
    )
    parser.add_argument(
        "--min-sequencing-time", dest="MIN_SEQ_TIME", type=int, default=None, help="Only output reads that where sequenced at or after the given time (minutes)."
    )
    parser.add_argument(
        "--start-time", dest="START_TIME", type=str, default=None, help="Starttime of the run as guppy time stamp (only required with --sequencing-time). If 'min' is given as argument the minimal time is detected automatically."
    )

    parser.add_argument(
        "--filter-id", dest="FILTER_ID", type=str, default=None, help="Only print reads with IDs present in file."
    )

    parser.add_argument(
        "--channels", dest="channels_input", type=str, nargs="*",
        help="List of individual `a`/ranges `a-b` of integer channel ids to print reads from. "
             "Ranges are inclusive on both sides. First `ch=\\d+` entry in header is considered.",
    )

    parser.add_argument(
        "--filter-as", dest="FILTER_AS", type=str, default=None, help="Adaptive sampling CSV file created by guppy."
    )

    parser.add_argument(
        "--filter-as-state", dest="FILTER_AS_STATE", default=None, action='append', choices=['unblock', 'unblock_hit_outside_bed', 'stop_receiving', 'no_decision'], help="Fiter reads by adaptiver sampling decision. (requires --filter-as)"
    )

    parser.add_argument(
        "--print-start-time", dest="PRINT_START_TIME", action="store_true", help="Print the minimal start_time of all fastq files"
    )

    parser.add_argument(
        "-n", "--max_n", dest="MAX_N", type=int, default=0, help="Stop after <max_n> reads"
    )
    parser.add_argument(
        "-b", "--max_mbp", dest="MAX_BP", type=int, default=0, help="Stop after <max_bp> mega base pairs"
    )
    parser.add_argument(
        "-r",
        "--recursive",
        dest="RECURSIVE",
        action="store_true",
        help="Search folders recursively",
    )

    parser.add_argument(
        "-d",
        "--dedup",
        dest="DEDUP",
        action="store_true",
        help="Remove duplicated reads.",
    )

    parser.add_argument(
        "--comments",
        choices=['forward', 'skip', 'wrap'],
        default='wrap',
        help="How to treat FASTQ header comments. "
             "`forward` them 'as is', `wrap` them into 'CO:Z:xxx' tag, `skip` them in the output.",
    )

    parser.add_argument(
        "FASTQ",
        nargs="+",
        type=str,
        help="FASTQ files or folders containing FASTQ files",
    )

    parser.add_argument(
        "-v",
        '--version',
        action='version',
        version='catfish ' + __version__,
        help="Print version",
    )


    args = parser.parse_args(argv)

    return args


def find_file_in_folder(
    folder,
    recursive=True,
    patterns=[
        "*.fastq",
        "*.fastq.gz",
        "*.fasta",
        "*.fasta.gz",
        "*.fa",
        "*.fa.gz",
        "*.fq",
        "*.fq.gz",
    ],
):
    if os.path.isfile(folder):
        return folder
    files = []

    glob = Path(folder).glob
    if recursive:
        glob = Path(folder).rglob

    for pattern in patterns:
        for file in glob(pattern):
            files.append(file)

    if len(files) == 0:
        logging.warning("Could not find {} files in {}".format(pattern, folder))
    return files


def parse_timestamp(time_str):
    try:
        time_obj = datetime.strptime(time_str,'%Y-%m-%dT%H:%M:%SZ')
    except ValueError:
        time_obj = datetime.strptime(time_str,'%Y-%m-%dT%H:%M:%S.%f%z').replace(tzinfo=None)
    return time_obj


def check_seq_time(comment, max_start_time,min_start_time):
    #This tests if the start time of the respective read is between
    #max_sequencing_time and min_sequencing_time
    #If one of the times is not given the condition is automatically considered true
    if (max_start_time == None and min_start_time == None):
        return True
    else:
        matchObj = re.search( r'start_time=([^ ]+)', comment, re.M|re.I)
        start_str = matchObj.group(1)
        start = parse_timestamp(start_str)

        bool_min=0
        bool_max=0

        if (max_start_time == None or start<=max_start_time):
            bool_max=1
        if (min_start_time == None or start>=min_start_time):
            bool_min=1

        if (bool_min == 1 and bool_max == 1):
            return True
        else:
            return False

def compare_start_time(comment,min_start_time):
    #Checks if a given min start time is smaller than the time of an entry
    #The smaller time is returned
    matchObj = re.search( r'start_time=([^ ]+)', comment, re.M|re.I)
    start_time_str = matchObj.group(1)
    start_time = parse_timestamp(start_time_str)

    if(min_start_time==0):
        return start_time
    elif(min_start_time<=start_time):
        return min_start_time
    else:
        return start_time


def parse_fastqs(filename, min_len=0, min_qscore=0, max_start_time=None, min_start_time=None, comments='wrap',
                 channels: Optional[Set[int]] = None):
    with pysam.FastxFile(filename) as fh:
        for entry in fh:
            if min_len and len(entry.sequence) < min_len:
                continue
            if (
                min_qscore
                and _compute_mean_qscore(entry.get_quality_array()) < min_qscore
            ):
                continue
            if not check_seq_time(entry.comment, max_start_time, min_start_time):
                continue
            if entry.comment and comments == 'wrap':
                entry.comment = "CO:Z:{}".format(entry.comment)
            elif comments == 'skip':
                entry.comment = None
            if channels is not None and get_channel_from_comment(entry.comment) not in channels:
                continue
            yield entry


def get_file_names(path, recursive):
    filenames = []
    if os.path.exists(path):
        filenames = [path]
        if os.path.isdir(path):
            logging.debug("Searching {} for FASTQ files".format(path))
            filenames = find_file_in_folder(path, recursive=recursive)
    else:
        logging.warning("Could not find {}".format(path))

    return filenames

def get_start_time(paths,recursive=False):
    """
    Only print the start time.
    This function automatically detects the minmal start_time of
    all the given fastq files


    :param paths: Input FASTQ files or folders containing FASTQ files
    :return: min_start_time
    """
    start = None
    max_start_time = None
    min_start_time=0
    for path in paths:
        filenames = get_file_names(path,recursive)
        for filename in filenames:
            with pysam.FastxFile(filename) as fh:
                for entry in fh:
                    min_start_time=compare_start_time(entry.comment,min_start_time)
    return min_start_time


def format_fq(paths, out_filename, min_len=0, min_qscore=0, max_n=0, max_bp=0, recursive=False, dedup=False,
              max_seq_time=0, min_seq_time=0, start_time=0, filter_read_ids_file=None, comments='wrap',
              filter_read_as_file=None, filter_read_as_decision=None, channels: Optional[Set[int]]=None):
    """
    Concatenate FASTQ files

    :param paths: Input FASTQ files or folders containing FASTQ files
    :param out_filename: Output FASTQ file
    :return: None
    """
    start = None
    max_start_time = None
    min_start_time = None

    keep_ids = None
    if filter_read_ids_file:
        keep_ids = set()
        with open(filter_read_ids_file, "r") as fh:
            for line in fh:
                read_id = line.strip()
                keep_ids.add(read_id)
            logging.info("Found {} read ids.".format(len(keep_ids)))

    if filter_read_as_file and filter_read_as_decision:
        if not keep_ids:
            keep_ids = set()
        with open(filter_read_as_file, "r") as fh:
            for line in fh:
                cols = line.strip().split(',')
                read_id = cols[4]
                read_decision = cols[6]
                if read_decision in filter_read_as_decision:
                    keep_ids.add(read_id)
            logging.info("Found {} reads with {}.".format(len(keep_ids), filter_read_as_decision))

    if start_time:
        if not start_time=="min":
            start = parse_timestamp(start_time)

            if(max_seq_time):
                max_start_time = start + timedelta(minutes=max_seq_time)
            if(min_seq_time):
                min_start_time = start + timedelta(minutes=min_seq_time)
        else:
            #This option allows to automatically use the minmal start_time of
            #all the given fastq files as input for --start-time
            start=get_start_time(paths,recursive)

            if(max_seq_time):
                max_start_time = start + timedelta(minutes=max_seq_time)
            if(min_seq_time):
                min_start_time = start + timedelta(minutes=min_seq_time)
    read_ids = set()

    n = 0
    n_bp = 0
    with open(out_filename, mode="w") as fout:
        for path in paths:
            filenames = get_file_names(path, recursive)
            logging.debug("Found {} files".format(len(filenames)))
            for filename in filenames:
                for entry in parse_fastqs(
                    filename, min_len=min_len, min_qscore=min_qscore, max_start_time=max_start_time,
                    min_start_time=min_start_time, comments=comments, channels=channels,
                ):
                    if dedup and entry.name in read_ids:
                        continue

                    if keep_ids is not None and entry.name not in keep_ids:
                        continue

                    fout.write(str(entry) + "\n")
                    if dedup:
                        read_ids.add(entry.name)
                    n += 1
                    n_bp += len(entry.sequence)
                    if max_n and n >= max_n or max_bp and n_bp > max_bp:
                        return


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    if(args.PRINT_START_TIME):
        min_start_time=get_start_time(args.FASTQ,args.RECURSIVE)
        print(min_start_time.strftime('%Y-%m-%dT%H:%M:%SZ'))
    else:

        if args.FILTER_AS and not args.FILTER_AS_STATE or not args.FILTER_AS and args.FILTER_AS_STATE:
            logging.error("--filter-as-state and --filter-as must either be both specified or both skipped.")
            return

        channels_set: Optional[Set[int]] = get_channels_set(args.channels_input) if args.channels_input else None

        format_fq(
            args.FASTQ,
            args.OUT,
            min_len=args.MIN_LEN,
            min_qscore=args.MIN_QSCORE,
            max_n=args.MAX_N,
            max_bp=args.MAX_BP * 1000 * 1000,
            recursive=args.RECURSIVE,
            dedup=args.DEDUP,
            max_seq_time=args.MAX_SEQ_TIME,
            min_seq_time=args.MIN_SEQ_TIME,
            start_time=args.START_TIME,
            filter_read_ids_file=args.FILTER_ID,
            filter_read_as_file=args.FILTER_AS,
            filter_read_as_decision=args.FILTER_AS_STATE,
            comments=args.comments,
            channels=channels_set,
        )


if __name__ == "__main__":
    main()
