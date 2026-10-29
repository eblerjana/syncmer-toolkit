# query

## Why

Right now there is no way to find out which unique syncmers (the ones `histogram` reports) fall inside a specific genomic interval without manually scanning raw `.1path` offsets. A `query` command is needed to answer "what unique syncmers exist in this region of this sequence" directly.

## What

A new `query` subcommand for the `syncmer-toolkit` CLI:

```
syncmer-toolkit query <.1path file> <.1khash file> <unique-syncmers.tsv> <fasta filename> <region> <outname>
```

Given the `.1path` file, the `.1khash` file, the unique-syncmer list produced by `histogram` (`<outname>_syncmers.tsv`), the name of one of the original input FASTA files, and a region `seqid:start-end` (0-based, half-open — same convention as the raw `o`-line offsets used by `distances`), it writes `<outname>` as a TSV with header `syncmer_ID\tposition\tcount\tsequence` listing every unique syncmer on the matching sequence whose offset falls in `[start, end)`. `sequence` is the syncmer's DNA sequence as traversed by the path: reverse-complemented when the path traverses it in reverse orientation (negative node id on the `z` line), forward otherwise — see "Orientation" below.

### Orientation (reverse-traversed syncmers)

A path's `z` line stores signed node ids: a negative id means the path traverses that syncmer in reverse orientation. `external/syng/kmerhash.h:51` already provides exactly this: `kmerHashSeq(KmerHash *kh, I64 i, char *buf)` returns the stored sequence reverse-complemented when `i < 0`. The core function must therefore look up each syncmer's sequence using its **signed** id straight from the `z` line (not `abs(sync_id)`, which is only used for the subset/count lookup), via a `KmerHash` loaded from the `.1khash` file with `kmerHashReadOneFile` (`external/syng/kmerhash.h:57`).

`seqid` identifies the sequence within the FASTA file either by name or by its 1-based position in the file — see "Sequence identification" below, since `I` lines (sequence names) are not always present in a `.1path` file.

### Sequence identification (`I` lines are optional)

The `D I` line carrying a path's sequence name is only written when the graph was built without `-noNames`; some `.1path` files have no `I` lines at all (e.g. the existing `tests/data/small.1path` fixture). The `O P` line's third field ("sequence number in file", 1-based) is always present regardless, so it is used as a fallback identifier.

Matching rule, applied per candidate path (i.e. per path whose `source_id` matches the resolved fasta filename):
- If the path has an `I` line, match requires `seqid == <I-line identifier>` (string equality).
- If the path has no `I` line, match requires `seqid` to parse as a positive integer equal to the path's sequence-number-in-file (third field of its `P` line).

This means a `.1path` file without any names is queried by passing the 1-based sequence index instead of a name (e.g. `2:0-100` for the second sequence in that FASTA file). If `seqid` doesn't parse as an integer and the candidate paths have no `I` lines, fail with a clear error telling the user the graph has no stored sequence names and an index must be given instead.

Verifiable by: running it against `tests/data` fixtures with a hand-computed expected region and comparing output.

## Constraints

### Must

- Follow the existing CLI dispatch pattern (`src/commands.cpp`, `src/main.cpp`): argc validation, usage string, `check_open_file` checks before processing.
- Follow the existing ONElib/syng parsing idioms in `src/syncmerstats.cpp` (P/I/z/o line handling), modeled closely on `compute_syncmer_distances_from_paths` (`src/syncmerstats.cpp:367`).
- Match `<fasta filename>` against the `<` reference lines stored in the `.1path` file (exact string match against `ipath->reference[i].filename`; `.count` gives the source_id).
- Match the region's `seqid` against candidate paths (those with `source_id` equal to the one found above) using the per-path fallback rule in "Sequence identification" above: by `I`-line name if present, otherwise by 1-based sequence-number-in-file.
- Treat the region as 0-based, half-open: a syncmer at offset `p` is included iff `start <= p < end`.
- Require the unique-syncmer subset file (histogram's `_syncmers.tsv`) as an explicit CLI argument, same pattern as `distances`'s `syncmerset` argument — do not recompute uniqueness from scratch.
- Load the actual syncmer sequences from the `.1khash` file via `kmerHashReadOneFile` and report each one via `kmerHashSeq` using the signed node id from the `z` line, so reverse-traversed syncmers come out reverse-complemented (see "Orientation").
- Add unit tests (core function in `test_syncmerstats.cpp`, CLI wrapper in `test_commands.cpp`) and document the command in `README.md`.

### Must Not

- No new external dependencies.
- Do not change the behavior or output format of `histogram` or `distances`.
- Do not add `.1gbwt` support for this command (path-only input, like `distances`).

### Out of Scope

- Querying multiple regions in one invocation.
- Any interactive/REPL-style query mode.
- Resolving ambiguity between a sequence literally *named* a number and index-based lookup — the fallback only engages when a path has no `I` line at all (see "Sequence identification").

## Current State

- `src/main.cpp` — usage banner + `if/else if` dispatch on `argv[1]` (`"histogram"`, `"distances"`); needs a `"query"` branch.
- `src/commands.hpp` / `src/commands.cpp` — per-command CLI wrappers (`create_histogram`, `compute_distances`): argc checks, usage string, `ends_with`/`check_open_file` helpers.
- `src/syncmerstats.hpp` / `src/syncmerstats.cpp` — core ONE-file logic. `compute_syncmer_distances_from_paths` (`src/syncmerstats.cpp:367`) is the closest analog: it opens the `.1path` file, reads a syncmer subset TSV into a lookup structure, and scans `P`/`z`/`o` lines per path — the same shape `query` needs, plus matching on the `<` reference and `I` lines.
- `external/syng/syng.h:62` (`syngSchemaText`) — schema reference: `O P` = path length/source_id/seq-index-in-file, `D I` = sequence identifier, `D z`/`D o` = syncmer node ids and their offsets (1:1), `<` lines (via `oneAddReference`/`OneReference`) = source filename → source_id mapping.
- `external/syng/kmerhash.h` — `KmerHash` holds the actual syncmer sequences; `kmerHashReadOneFile(OneFile*)` loads one from an opened `.1khash` file, `kmerHashSeq(KmerHash*, I64 i, char* buf)` returns the sequence for id `i`, reverse-complemented automatically when `i < 0`. `src/syncmerstats.cpp`'s current `get_max()` (line 43) only reads the `t` line for `nSyncmers`/bounds-checking and does not load actual sequences — `query` needs the fuller `kmerHashReadOneFile` load instead/in addition.
- `tests/data/small.1path` has **no** `I` lines currently (verified with `ONEview`) — this is already a ready-made fixture for the index-fallback matching path. A new fixture with `I` lines is still needed to test name-based matching. `tests/data/README.md` documents how existing hand-edited fixtures were produced (dump with `ONEview`, edit text, reconvert with `ONEview -b`).
- `tests/test_helpers.hpp` provides `Args`, `data_path`, `output_path`, `read_file`, `file_exists` used by all CLI-level tests.

## Tasks

### T1: Core query function

What: Implement `compute_syncmer_query_from_paths(pathfile_path, khashfile_path, syncmerset_path, fasta_filename, region, outfile_path)`. Parse `region` (`seqid:start-end`, 0-based half-open; error on malformed input or `start >= end`). Open the `.1path` file and resolve `fasta_filename` to a `source_id` via its reference list (error listing available filenames if no match). Load the unique-syncmer subset TSV into a `sync_id -> count` lookup (same TSV parsing as `compute_syncmer_distances_from_paths`). Scan `P`/`I`/`z`/`o` lines; for each path with matching `source_id`, apply the per-path fallback rule from "Sequence identification" (match `seqid` against its `I` line if present, else against its `P` line's sequence-number-in-file) to find the one matching path (error if none match, with a message noting whether the file has sequence names at all). Load a `KmerHash` from the `.1khash` file via `kmerHashReadOneFile`. For that path's `z`/`o` pairs, write a row (`syncmer_ID`, `position`, `count`, `sequence`) for each syncmer that is both in the subset and has offset in `[start, end)`, in increasing offset order; `sequence` comes from `kmerHashSeq(kh, signed_sync_id, buf)` using the **signed** id straight from the `z` line so reverse-traversed syncmers are reverse-complemented. Destroy the `KmerHash` with `kmerHashDestroy` before returning.
Files: `src/syncmerstats.hpp`, `src/syncmerstats.cpp`
Verify: covered by T4 unit tests.

### T2: CLI wiring

What: Add `compute_query(int argc, char* argv[])` to `src/commands.hpp`/`src/commands.cpp`: validate `argc == 8`, print a usage string on mismatch, `check_open_file` the `.1path`/`.1khash` files, call `compute_syncmer_query_from_paths`. Add a `"query"` branch to the dispatch in `src/main.cpp` and extend its usage banner.
Files: `src/commands.hpp`, `src/commands.cpp`, `src/main.cpp`
Verify: `./build/src/syncmer-toolkit query` with missing/extra args prints usage and exits 1 (manual).

### T3: Test fixtures with sequence names

What: Create a `.1path` fixture with `I` lines, following the hand-edit workflow in `tests/data/README.md` (dump `small.1path` with `ONEview`, insert an `I` line after each `P` line with a sequence name, reconvert to binary with `ONEview -b -t path`). Document the new fixture and the names/regions it supports in `tests/data/README.md`.
Files: `tests/data/README.md`, new `tests/data/small_with_names.1path`
Verify: `./ONEview tests/data/small_with_names.1path` shows the expected `I` lines.

### T4: Unit tests for the core function

What: Add `TEST_CASE("Test syncmerstats::compute_syncmer_query_from_paths")` covering: happy path with known expected syncmer IDs/positions/sequences using the named fixture; happy path using `small.1path` (no `I` lines) addressed by 1-based sequence index; a forward-traversed syncmer (positive `z` id) reporting the sequence as stored; a reverse-traversed syncmer (negative `z` id, e.g. the `-4`/`-5` entries in `small.1path`'s path 1) reporting the reverse complement; region before/after all syncmers on a path (header-only output); unknown fasta filename; unknown sequence name/index under a known fasta filename; non-numeric `seqid` against a file with no `I` lines; malformed region string; a syncmer inside the region but excluded because it's absent from the subset file.
Files: `tests/test_syncmerstats.cpp`
Verify: `./build/tests/unit_tests` (or `ctest --test-dir build --output-on-failure`)

### T5: CLI-level tests

What: Add `TEST_CASE("Test commands::compute_query")` mirroring the `compute_distances` tests: valid input produces the expected output file; too few/too many arguments; path/khash file cannot be opened; unknown fasta filename.
Files: `tests/test_commands.cpp`
Verify: same as T4.

### T6: Documentation

What: Add a `### query` section under `## Commands` in `README.md` describing the arguments, the region convention (0-based, half-open), and the name-vs-index sequence identification fallback, and extend the `## Usage` block with the new command line.
Files: `README.md`
Verify: manual review.

## Validation

- `cmake --build build -j`
- `ctest --test-dir build --output-on-failure` (runs `unit_tests`, includes new T4/T5 cases)
- Manual check: `./build/src/syncmer-toolkit query tests/data/small_with_names.1path tests/data/small.1khash tests/data/small_subset.tsv seq.fa <seqname>:0-6 /tmp/query_out` and confirm `/tmp/query_out` matches hand-computed expectations, including that the `sequence` column is reverse-complemented for syncmers traversed with a negative `z`-line id.
