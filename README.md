# bdbwt-mem-chaining

Requires [BDBWT](https://github.com/algbio/bdbwt) and it's dependencies which are included in /include (BDBWT) and /sdsl-lite (SDSL) folders.

Similar to BDBWT, requires building sdsl to be built:
```
cd sdsl-lite
sh ./install.sh
```
Afterwards the rest of the code can be made with `make all` followed by executing the created executable `./main` 

As of current commit (2020-05-28), the algorithm should be working as intended on both naive (space-heavy) output and more efficient batched output style.
