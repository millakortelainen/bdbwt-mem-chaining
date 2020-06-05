# bdbwt-mem-chaining

## Required libraries:
[BDBWT](https://github.com/algbio/bdbwt) and it's dependencies which are included in /include (BDBWT) and /sdsl-lite (SDSL) folders.

[Edlib](https://github.com/Martinsos/edlib) sequence alignment libary is used, and included under MIT licence. Folder "edlib" (and it's subfolders) are taken directly from said library and contain necessary source code and header file. Edlib works out of the box, and does not warrant addinational installation on this programs implementation.

Similar to BDBWT, requires building sdsl to be built:
```
cd sdsl-lite
sh ./install.sh
```
Afterwards the rest of the code can be compiled with `make all`, which creates the executable `.\main`.

As of current commit (2020-06-05), the algorithm program is able to quickly find all MEM's between two strings, and return an optimal chaining path between them.  

