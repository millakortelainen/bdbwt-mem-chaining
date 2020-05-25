# bdbwt-mem-chaining

Requires [BDBWT](https://github.com/algbio/bdbwt) and it's dependencies which are included in /include (BDBWT) and /sdsl-lite (SDSL) folders.

Similar to BDBWT, requires building sdsl to be built:
```
cd sdsl-lite
sh ./install.sh
```
Afterwards the rest of the code can be made with `make all` followed by executing the created executable `./mem` 

As of current commit (2020-05-25), the implementation, according to manual testing, finds all maximal exact matches between two strings using two bi-directional burrows-wheeler transforms. But relies on non-ideal tricks to yield said result. Namely by considering `a == c && b == d` instead of `a != c && b != d` in the subroutine `bwt_mem2_subroutine` contradictory to the pseudocode the implementation is based on.
