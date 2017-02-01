# backwardMEM

Porting of *backwardMEM* by Ohlebusch *et al.* [1] to sdsl-lite [2]. The original source can be found [here](https://www.uni-ulm.de/in/theo/research/seqana/).

*backwardMEM* is an algorithm based on enhanced compressed suffix arrays for computing Maximal Exact Matches between two strings.

Tested with gcc version 5.4.0 and 5.4.1

Clone and test with:
```bash
git clone --recursive git@github.com:ldenti/backwardMEM.git
cd backwardMEM
make prerequisites
make
./bin/backwardMEM -l=3 ./example/T.fa ./example/P.fa
```

[1] Ohlebusch Enno, Simon Gog, and Adrian Kügel. "Computing matching statistics and maximal exact matches on compressed full-text indexes." *International Symposium on String Processing and Information Retrieval.* Springer Berlin Heidelberg, 2010.

[2] Gog, Simon, et al. "From theory to practice: Plug and play with succinct data structures." *International Symposium on Experimental Algorithms.* Springer International Publishing, 2014.
