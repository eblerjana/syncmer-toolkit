# Toolkit to analyze syncmer graphs

This code is experimental and currently only used to play around and analyze syncmere graphs built with syng. Use at your own risk.

## Installation

``` bat
git clone --recurse-submodules https://github.com/eblerjana/syncmer-toolkit.git
cd syncmer-toolkit
mkdir build && cd build
cmake ..
make
```


## Usage

``` bat
syncmer-toolkit histogram <.1path file> <.1khash file> <outname>
syncmer-toolkit histogram <.1gbwt file> <.1khash file> <outname>
syncmer-toolkit distances <.1path file> <.1khash file> <syncmerset> <outname>
syncmer-toolkit query <.1path file> <.1khash file> <unique-syncmers.tsv> <fasta filename> <seqid:start-end> <outname>

```
## Commands

### histogram

Computes a histogram over all "unique" syncmeres that occur at most once per input file provided to syng to build the graph. In other words, each such syncmere is allowed to occur at most once across all paths that are contained in the same input file. The histogram then counts in how many files such syncmeres are seen.

Arguments:
- `<.1path file>` or `<.1gbwt file>`: ``.1path`` or ``.1gwt`` files produced by syng.
- `<.1khash file>`: ``.1khash`` file produced by syng.
- `<outname>`: name of the output file TSV file, with columns `syncmer_ID` and `total_count`.


Currently, running this command with the `.1gbwt` file seems very slow, as the traversal of paths in the GBWT takes a lot of time. If available, use the `.1path` file instead.


### distances

For the subset of unique syncmers, this command parses the ``.1path`` file again to determine the distances between subsequent syncmeres along the paths. For each source file, it then outputs a histogram of distances.

Arguments:
- `<.1path file>`: ``.1path`` file produced by syng.
- `<.1khash file>`: ``.1khash`` file produced by syng.
- `<unique-syncmers.tsv>`:  the `<outname>_syncmers.tsv` file produced by `histogram`.
- `<outname>`: name of the output file TSV file, with columns `source_id`, `distance`, `count`.

### query

Given a genomic interval on one of the original input FASTA files, this command reports the unique syncmers (the ones `histogram` would report, passed in via `<unique-syncmers.tsv>`) that fall into that interval.

Arguments:
- `<.1path file>`, `<.1khash file>`: the same graph files used by the other commands.
- `<unique-syncmers.tsv>`: the `<outname>_syncmers.tsv` file produced by `histogram`.
- `<fasta filename>`: the name of one of the FASTA files originally given to `syng` to build the graph, matched against the source filenames recorded in the `.1path` file.
- `<seqid:start-end>`: the region to query. `start` and `end` are 0-based and half-open, i.e. a syncmer at offset `p` is included iff `start <= p < end`.
- `<outname>`: path of the output TSV, with columns `syncmer_ID`, `position`, `count` and `sequence`.

`seqid` identifies the sequence within the FASTA file, either by name or by its 1-based position in the file:
- If the graph was built with sequence names recorded (the default; disabled in `syng` via `-noNames`), `seqid` must match that name.
- If the graph has no stored sequence names for that FASTA file, `seqid` must instead be given as the sequence's 1-based position within the file (e.g. `2:1000-2000` for the second sequence in that FASTA file).

If a syncmer is traversed in reverse orientation by the path, `sequence` is reported as the reverse complement.
