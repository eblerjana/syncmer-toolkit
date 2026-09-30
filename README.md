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

```
## Commands

### histogram

Computes a histogram over all "unique" syncmeres that occur at most once per input file provided to syng to build the graph. In other words, each such syncmere is allowed to occur at most once across all paths that are contained in the same input file. The histogram then counts in how many files such syncmeres are seen.

### distances

For the subset of unique syncmers, this command parses the ``.1path`` file again to determine the distances between subsequent syncmeres along the paths.  
