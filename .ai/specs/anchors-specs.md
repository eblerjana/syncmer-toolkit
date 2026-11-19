# anchors

## Why

`histogram` reports the unique syncmers (at most once per source file) together with how many source files contain them, but there is currently no way to get the subset of these syncmers that is shared by *all* haplotypes, laid out along one chosen reference haplotype. These shared syncmers act as anchors: they mark the regions where all haplotypes agree, i.e. everything outside the top-level bubbles of the pangenome. An `anchors` command is needed to produce this list, in reference order, so that downstream analyses can build on it.

## What

A new `anchors` subcommand for the `syncmer-toolkit` CLI:

```
syncmer-toolkit anchors <.1path file> <unique-syncmers.tsv> <reference fasta filename> <outname>
```

Given the `.1path` file, the unique-syncmer list produced by `histogram` (`<outname>_syncmers.tsv`) and the name of the reference haplotype, it writes `<outname>` as a TSV listing the anchor syncmers, one per line, in the order in which they occur along the reference haplotype.

### Haplotypes and the reference haplotype

Each source (FASTA) file given to `syng` is one haplotype; this is the same unit that `histogram`'s `total_count` counts over. The reference haplotype is identified by its FASTA filename, matched exactly against the `<` reference lines stored in the `.1path` file (`ipath->reference[i].filename`, with `.count` giving the `source_id`), the same way `query` resolves `<fasta filename>`. The total number of haplotypes is the number of source references in the `.1path` file (`oneReferenceCount(ipath)`).

### Anchor definition

Anchors are selected in two steps:

1. **Candidates:** a syncmer is a candidate iff it is listed in `<unique-syncmers.tsv>` (so it occurs at most once per haplotype) **and** its `total_count` equals the total number of haplotypes (so it occurs exactly once in every haplotype).
2. **Consistent order and orientation:** from the candidates, only the subset whose order and orientation are the same in all haplotypes is kept as anchors. Order and orientation are compared relative to the overall orientation of the contig the candidate lies on, so a contig stored reverse-complemented in a haplotype's FASTA file is not penalized. A candidate that occurs out of order or in a different orientation (relative to its contig) in at least one haplotype (e.g. due to a rearrangement or inversion) is not an anchor.

Anchors are therefore a subset of the input unique-syncmer list. The concrete algorithm for step 2 (including which consistent subset is selected) is described in "Anchor selection algorithm" below.

### Contigs

Each haplotype may consist of several contigs (sequences in the FASTA file, i.e. several paths with the same `source_id`), and contigs can be arranged arbitrarily within a FASTA file. Contig breaks in the non-reference haplotypes must therefore not be penalized: each contig of a non-reference haplotype is compared to the reference on its own, so the order of two candidates only has to match if they are located on the same contig of that haplotype. Across different contigs of a non-reference haplotype, the order does not need to be the same.

Contig breaks in the reference are not handled the same way: the reference contigs are concatenated in file order (see "Anchor selection algorithm"), so if a non-reference contig spans several reference contigs that are arranged differently, anchors can be lost (see "Out of Scope").

### Anchor selection algorithm

- Step 1:
    - parse the reference path(s) in the `.1path` input file. Process the reference paths sorted by their sequence number in file (third field of the `P` line), since paths are not necessarily stored in file order (cf. `small_unordered.1path`). Store a mapping of `syncmer_id -> index` for all candidates on the reference, where the index increases along the concatenated reference contigs.
    - additionally store the orientation of each candidate in the reference, as well as its contig and its position on the contig
    - initialize the anchors list with all candidates (see "Anchor definition")

- Step 2:
    - for each path in the input `.1path` file that is not a reference path (i.e. for each contig of each non-reference haplotype), perform the steps listed below
    - parse the sequence of syncmer IDs of the haplotype path from the input file
    - store the list of reference indices of the candidates on this path, in path order (called `H` in the following; looked up in the mapping from Step 1). `H` always contains **all** candidates on the path, including those already removed from the anchors list by a previously processed path, so that the result does not depend on the order in which paths are processed.
    - while parsing, also count the candidates on this path that are in the same and in reversed orientation compared to the reference
    - if strictly more candidates are in reversed orientation, the contig is considered reversed: reverse `H`. On a tie, the contig is considered forward.
    - remove every candidate whose orientation on this path, flipped if the contig is considered reversed, differs from its orientation in the reference from the anchors list
    - iterate `H` from the front and compute for each position `i` the running maximum `prefix_max[i] = max(H[0..i])`
    - iterate `H` from the back and compute for each position `i` the running minimum `suffix_min[i] = min(H[i..end])`
    - for each position `i`, unless `H[i] == prefix_max[i] && H[i] == suffix_min[i]`, remove the respective candidate from the anchors list

- Step 3:
    - output all remaining anchors in reference order, together with their contig information and position (see "Output order and format")

Note: the prefix-max/suffix-min rule deliberately removes *every* candidate involved in an order violation, not a minimal set. E.g. for `H = [1, 2, 10, 3, 4]`, only 1 and 2 are kept; 3 and 4 are removed as well because they come after 10. This conservative behavior is intended and must not be replaced by selecting a largest consistent subset (e.g. longest increasing subsequence).

### Output order and format

The output lists the anchors in the order in which they occur along the reference haplotype: ordered by the reference contig's sequence number in file (third field of the `P` line), then by increasing position within each contig. The output TSV has the header `syncmer_ID\tcontig\tposition`:
- `syncmer_ID`: the **signed** node id from the reference path's `z` line, so the orientation of the anchor on the reference is retained (consistent with `query`, see commit `cee6f1b`).
- `contig`: the reference contig the anchor is located on: the path's `I`-line name if present, otherwise its 1-based sequence number in file (same fallback as `query`).
- `position`: the syncmer's position on that reference contig (raw `o`-line value, 0-based).

Verifiable by: running it against a `tests/data` fixture in which some syncmers are shared by all source files, some of them in a different order or orientation in one haplotype, with hand-computed expected anchors, and comparing output.

## Constraints

### Must

- Follow the existing CLI dispatch pattern (`src/commands.cpp`, `src/main.cpp`): argc validation, usage string, `check_open_file` checks before processing.
- Follow the existing ONElib/syng parsing idioms in `src/syncmerstats.cpp` (P/I/z/o line handling), modeled closely on `compute_syncmer_distances_from_paths` (`src/syncmerstats.cpp:367`) and `compute_syncmer_query_from_paths` (`src/syncmerstats.cpp:566`).
- Resolve the reference FASTA filename to a `source_id` via the `.1path` file's reference list; on no match, fail with an error listing the available filenames (same as `query`).
- Require the unique-syncmer subset file (histogram's `_syncmers.tsv`) as an explicit CLI argument and read both its `syncmer_ID` and `total_count` columns — do not recompute uniqueness or counts from the paths.
- Do not take a `.1khash` file: unlike `distances` and `query`, `anchors` needs neither the number of syncmers nor their sequences. Store the candidate and reference-index lookups in hash-based containers keyed by syncmer ID (e.g. `unordered_map`), so that no upper bound on the IDs is needed. Consequently, syncmer IDs in the TSV are not bounds-checked: an ID that does not occur in the graph simply never matches a path and is ignored.
- Only output candidates whose order and orientation (relative to their contig) are the same in all haplotypes, as described in "Anchor definition" and "Anchor selection algorithm".
- Do not penalize contig breaks or reverse-complemented contigs in the non-reference haplotypes (see "Contigs").
- Keep the conservative removal rule of the "Anchor selection algorithm" (remove every candidate involved in an order violation).
- Write the anchors in reference order as described in "Output order and format".
- Add unit tests (core function in `test_syncmerstats.cpp`, CLI wrapper in `test_commands.cpp`) and document the command in `README.md`.

### Must Not

- No new external dependencies.
- Do not change the behavior or output format of `histogram`, `distances` or `query`.
- Do not add `.1gbwt` support for this command (path-only input, like `distances` and `query`).

### Out of Scope

- Computing or reporting the bubbles between consecutive anchors.
- If haplotype sequences are broken into contigs in reference, but joined in one contig in haplotype, some anchors might be missed (contig breaks are only handled for non-reference haplotypes, see "Contigs").
- Restricting the output to a region or a subset of the reference's sequences.
- Using a reference that is not one of the haplotypes the graph was built from.

## Notes for Downstream Use

These points do not change what `anchors` has to implement. They document properties of its output that a future command must take into account. That command would report a VCF-like record for every region between consecutive anchors in which the haplotype sequences differ.

### N1: Non-candidate syncmers can connect intervals in the graph

The order check only considers candidates. A syncmer that is unique but not present in all haplotypes, or that occurs more than once in a haplotype, is never checked and can lie in different anchor intervals in different haplotypes. Example with haplotypes A (reference), B and C. Syncmers 1–6 are candidates, X occurs only in A and B (`total_count` 2), so X is not a candidate:

```
A (ref): 1  X  2  3  4     5  6
B:       1     2  3  4  X  5  6
C:       1     2  3  4     5  6
```

All of 1–6 are reported as anchors. In the graph, however, the edges `1→X` (from A) and `X→5` (from B) create the walk `1 → X → 5`, which bypasses 2, 3 and 4. So in the strict graph sense, 2–4 are not cut points and the true top-level bubble spans 1 to 5.

Consequence: the anchors still split every haplotype consistently into non-overlapping intervals, so a sequence-based comparison between consecutive anchors is unaffected. In the example it yields two independent records (an insertion of X in A between 1 and 2, and in B between 4 and 5). What gets lost is the information that both insertions are the same sequence (e.g. a transposed segment). This only matters if such events are to be annotated or linked later. If true graph bubbles are needed, a possible extension is to record, for each non-candidate syncmer, the anchor interval it falls into per haplotype, and to merge intervals (i.e. remove the anchors between them) whenever a syncmer falls into different intervals.

### N2: Contig ends are not constrained

The order check runs per contig of each non-reference haplotype. The part of a contig before its first anchor (or after its last) is only required to lie before (or after) that anchor, and may belong to any earlier (or later) interval. A contig without any anchors is not constrained at all. The same applies to contig breaks in the reference (see "Out of Scope").

Consequence: for a downstream VCF-like output, a haplotype whose contig starts or ends inside an interval has only a partial allele for it, and an interval spread over several contigs has no contiguous allele in that haplotype. Such cases must be reported as missing or partial (e.g. `.` as genotype) rather than as a regular allele.

## Current State

- `src/main.cpp` — usage banner + `if/else if` dispatch on `argv[1]` (`"histogram"`, `"distances"`, `"query"`); needs an `"anchors"` branch and an updated error message listing valid program names.
- `src/commands.hpp` / `src/commands.cpp` — per-command CLI wrappers (`create_histogram`, `compute_distances` at line 93, `compute_query` at line 133): argc checks, usage string, `ends_with`/`check_open_file` helpers.
- `src/syncmerstats.hpp` / `src/syncmerstats.cpp` — core ONE-file logic. `compute_syncmer_distances_from_paths` (line 367) shows the subset-TSV parsing (currently only reads the `syncmer_ID` column; its `get_max` bounds check against the `.1khash` file is not needed for `anchors`); `compute_syncmer_query_from_paths` (line 566) shows resolving a FASTA filename to a `source_id` via `oneReferenceCount`/`ipath->reference[i]` and the `I`-line / sequence-index handling.
- `histogram` writes `<outname>_syncmers.tsv` with header `syncmer_ID\ttotal_count` (`src/syncmerstats.cpp:166`), where `total_count` is the number of source files containing the syncmer.
- `tests/data/small.1path` contains two source files (`seq.fa`, `seq2.fa`), but no syncmer occurs in both (`small_subset.tsv` has `total_count` 1 for all entries), so it yields no anchors. `small_repeats.1path` has syncmer 3 in multiple source files. A new fixture with several syncmers shared by both source files is needed. `tests/data/README.md` documents the hand-edit workflow for fixtures (dump with `ONEview`, edit text, reconvert with `ONEview -b`).
- `tests/test_helpers.hpp` provides `Args`, `data_path`, `output_path`, `read_file`, `file_exists` used by all CLI-level tests.

## Tasks

### T1: Core anchors function

What: Implement `compute_syncmer_anchors_from_paths(pathfile_path, syncmerset_path, reference_filename, outfile_path)`. Open the `.1path` file, resolve `reference_filename` to a `source_id` (error listing available filenames if no match) and determine the number of haplotypes via `oneReferenceCount`. Read the unique-syncmer TSV into a hash-based `sync_id -> total_count` lookup (error on malformed lines, no bounds check, see "Must"); mark as candidates the IDs whose `total_count` equals the number of haplotypes. Scan `P`/`I`/`z`/`o` lines and select the anchors from the candidates following "Anchor selection algorithm". Write the anchors in reference order to the TSV with header `syncmer_ID\tcontig\tposition`, using the signed ID from the reference path's `z` line and the contig name/index as described in "Output order and format". Print a summary (number of haplotypes, number of candidates, number of anchors) to `cerr`, like the other commands.
Files: `src/syncmerstats.hpp`, `src/syncmerstats.cpp`
Verify: covered by T4 unit tests.

### T2: CLI wiring

What: Add `compute_anchors(int argc, char* argv[])` to `src/commands.hpp`/`src/commands.cpp`: validate `argc == 6`, print a usage string on mismatch, `check_open_file` the `.1path` file, call `compute_syncmer_anchors_from_paths`. Add an `"anchors"` branch to the dispatch in `src/main.cpp`, extend its usage banner and the "program name must be specified as ..." error.
Files: `src/commands.hpp`, `src/commands.cpp`, `src/main.cpp`
Verify: `./build/src/syncmer-toolkit anchors` with missing/extra args prints usage and exits 1 (manual).

### T3: Test fixtures with shared syncmers

What: Create a `.1path` fixture (e.g. `small_shared.1path`) derived from `small.1path` in which several syncmers occur exactly once in both source files, including at least one traversed in reverse orientation on the reference, one shared syncmer that is not unique (occurs twice in one source file), at least one candidate that occurs in a different order in the non-reference haplotype, at least one candidate that occurs in a different orientation in the non-reference haplotype (both must be excluded), haplotypes with multiple contigs that are arranged in a different order in the two source files, and a contig that is stored reverse-complemented in the non-reference haplotype (both must not lead to exclusions). Additional fixtures (e.g. `small_shared_ref_break.1path`, `small_shared_tie.1path`) cover a reference contig split across two contigs of the non-reference haplotype, a non-reference contig spanning two reference contigs that are arranged differently (the known loss of anchors from "Out of Scope"), and a contig with a tie between forward and reversed candidates. Create a matching `_syncmers.tsv` fixture (e.g. `small_shared_syncmers.tsv`) with the correct `total_count` values, as `histogram` would produce it. Document both fixtures and the expected anchors per reference in `tests/data/README.md`.
Files: `tests/data/README.md`, new `tests/data/small_shared*.1path`, new `tests/data/small_shared*_syncmers.tsv`
Verify: `./ONEview tests/data/small_shared.1path` shows the expected paths; running `histogram` on it reproduces `small_shared_syncmers.tsv`.

### T4: Unit tests for the core function

What: Add `TEST_CASE("Test syncmerstats::compute_syncmer_anchors_from_paths")` covering: happy path with known expected anchors, using each of the two source files as reference (same anchor set, orientation reported relative to the respective reference); a reverse-traversed anchor reported with a negative ID; a candidate occurring in a different order in one haplotype being excluded; a candidate occurring in a different orientation in one haplotype being excluded; contigs arranged in a different order between haplotypes not causing exclusions; a reverse-complemented contig in the non-reference haplotype not causing exclusions; a reference contig split across two non-reference contigs not causing exclusions; a non-reference contig spanning two differently arranged reference contigs losing the expected anchors (documents the "Out of Scope" behavior); a tie between forward and reversed candidates treating the contig as forward; the conservative removal rule (e.g. reference indices `[1, 2, 10, 3, 4]` on a non-reference contig keep only 1 and 2); reference paths not stored in file order (cf. `small_unordered.1path`); output ordered by reference contig, then position, with `contig` reported by name (with `I` lines) and by 1-based index (without); a non-unique shared syncmer being excluded; no anchors (e.g. `small.1path` + `small_subset.tsv`) giving header-only output; unknown reference filename; a syncmer ID in the TSV that does not occur in the graph being ignored (e.g. `small_subset_invalid.tsv`); malformed subset file; empty subset file.
Files: `tests/test_syncmerstats.cpp`
Verify: `./build/tests/unit_tests` (or `ctest --test-dir build --output-on-failure`)

### T5: CLI-level tests

What: Add `TEST_CASE("Test commands::compute_anchors")` mirroring the `compute_query` tests: valid input produces the expected output file; too few/too many arguments; path file cannot be opened; unknown reference filename.
Files: `tests/test_commands.cpp`
Verify: same as T4.

### T6: Documentation

What: Add a `### anchors` section under `## Commands` in `README.md` describing the arguments, the anchor definition (exactly once in every haplotype, same order and orientation in all haplotypes, contig breaks not penalized), the reference ordering and the output columns, and extend the `## Usage` block with the new command line.
Files: `README.md`
Verify: manual review.

## Validation

- `cmake --build build -j`
- `ctest --test-dir build --output-on-failure` (runs `unit_tests`, includes new T4/T5 cases)
- Manual check: `./build/src/syncmer-toolkit anchors tests/data/small_shared.1path tests/data/small_shared_syncmers.tsv seq.fa /tmp/anchors_out` and confirm `/tmp/anchors_out` matches the hand-computed anchors in reference order, with negative IDs for anchors traversed in reverse on the reference and without the candidates that are out of order or differently oriented in the other haplotype.
