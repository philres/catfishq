# catfishq
Takes paths to an arbritary number of zipped and unzipped FASTQ files and/or folders containing zipped or unzipped FASTQ files, concatenates them and prints them to standard out (default) or an unzipped output file.

Supported file extensions are: `'*.fastq', '*.fastq.gz', '*.fasta', '*.fasta.gz', '*.fa', '*.fa.gz', '*.fq', '*.fq.gz'`

May also be used to filter FQ reads by read ID, read length, q-score, and min/max sequencing time.


# Install
``` bash
pip install catfishq
```

# Examples

Check full command list:

```bash
$ catfishq --help;
```

Merge all FQ files within a target directory:

```bash
$ catfishq test/ > test.fastq;
```

Merge all FQ files within a target directory and its sub-directories (recurive):

```bash
$ catfishq -r test/ > test.fastq;
```

Merge the first 1000 reads:

```bash
$ catfishq -n 1000 test/ > test_1st_1000.fastq;
```

Merge reads with a length >=50bp and a q-score >=10:
```bash
$ catfishq -l 50 -q 10 -l test/ > test_filt.fastq;
```

Merge reads collected <60mins from sequencing start:

```bash
catfishq --min-sequencing-time 0 --max-sequencing-time 60 test/ > test_60_min.fastq;  #merge reads
```

Note that when looping catfishq over multiple folders from the same run, it is quicker to grab the start time via `--print-start-time` and providing it to catfishq via `--start-time "$timestamp"`.

```bash
$ t0="$(catfishq --print-start-time test1)";
catfishq --max-sequencing-time 60 --start-time "$t0" test1/ > test1_60_min.fastq;
catfishq --max-sequencing-time 60 --start-time "$t0" test2/ > test2_60_min.fastq;
```
