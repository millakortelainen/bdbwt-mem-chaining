# bdbwt-mem-chaining
This project contains the implementation side of an on-going Master's Thesis project. All code, excluding external libraries noted below, is written for, and as part of the thesis.

This algorithm takes two sequences as an input in the form of two fasta files. And computes all maximal exact matches of length higher than given threshold (see configuration details). And computes the edit distance between the two using chaining.

For further technical details see [TBA]
## Used external libraries:
All external libraries are packed into this repository

[BDBWT](https://github.com/algbio/bdbwt) and it's dependencies which are included in /include (BDBWT) and /sdsl-lite (SDSL) folders.
[SDSL-lite](https://github.com/simongog/sdsl-lite)

[Edlib](https://github.com/Martinsos/edlib) sequence alignment libary is used, and included under MIT licence. Folder "edlib" (and it's subfolders) are taken directly from said library and contain necessary source code and header file. Edlib works out of the box, and does not warrant additional  installation on this programs implementation.

[OpenMP](https://www.openmp.org//) is used for implementing parallelization for the algorithms.

## Installation and Usage
Similar to BDBWT, the project requires requires building SDSL to function which can be done as follows, starting from the root folder of the project:
```
cd sdsl-lite
sh ./install.sh
```
Afterwards the rest of the code can be compiled with `make all`, which creates the executable `.\main`. The algorithm itself is executed by giving it an configuration file as an parameter
```
./main config
```
The configuration file contains all the options available for running the algorithm (see example file: 'config')
### Configuration details
|Field | Description |
- | -	
|Text1| Path to the first fasta file.|
|Index1| Integer Index of the sequence contained in the fasta (defaults to the last found if index is higher than actual amount in the file).
|Text2| Path to the second fasta file.|
|Index2| Integer index analogously to index1.|
|Mode | Integer Defining the search mode of the algorithm, 0: BDBWT, 1: Minimizer, 2: Hybrid, 3: Edlib-only.|
|Depth | Integer denoting the minimum match size, doubles as k-mer size 'k'|
|Window | Integer denoting the minimizer window size, typically desireable to have Window = Depth|
|Mergers | Integer denoting iteration count for merging k-mers. Default: 0 |
|BWTThrd | Integer denoting whether to thread BWT MEM search or not, 0: off 1: on.|
|VerbCA | Integer denoting whether to print verbose representation of each chain and absent segments 0:off 1: on|
|RawChain| Integer denoting whether to print chains as raw tuples (i,j,d) 0:off, 1: on|
|strChain| Integer denoting whether to print chains as substrings from the texts 0:off, 1: on|


As of current commit (2020-08-27), the algorithm is able to efficiently (especially with minimizers) compute edit distance between two strings. Given matches appearing in relatively similar positions, optimal alignment can be gained faster than with Edlib alone.
