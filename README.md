# backwardMEM

Porting of backwardMEM by Ohlebusch *et al.* ([link](https://www.uni-ulm.de/in/theo/research/seqana/)) to sdsl-lite.
Tested with gcc version 5.4.0

Clone and test with:
```bash
git clone git@github.com:
cd backwardMEM
make prerequisites
make
./bin/backwardMEM -l=3 ./example/T.fa example/P.fa
```
