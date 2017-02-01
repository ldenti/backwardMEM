# backwardMEM
## A tool for computing maximal exact matches
Porting of *backwardMEM* by Ohlebusch *et al.* to sdsl-lite (the original source can be found [here](https://www.uni-ulm.de/in/theo/research/seqana/)).

*backwardMEM* is a tool for computing Maximal Exact Matches between two strings.

Tested with gcc version 5.4.0 and 5.4.1

Clone and test with:
```bash
git clone git@github.com:ldenti/backwardMEM.git
cd backwardMEM
make prerequisites
make
./bin/backwardMEM -l=3 ./example/T.fa ./example/P.fa
```
