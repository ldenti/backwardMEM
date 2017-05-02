/* sdsl - succinct data structures library
   Copyright (C) 2009 Simon Gog

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

#include "../libs/sdsl-lite/include/sdsl/suffix_trees.hpp"
#include "../libs/sdsl-lite/include/sdsl/suffix_array_algorithm.hpp"
#include "../libs/sdsl-lite/include/sdsl/csa_wt.hpp"
#include "../src/testutils.hpp"
#include "../libs/sdsl-lite/include/sdsl/io.hpp"
#include "../libs/sdsl-lite/include/sdsl/construct.hpp"

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <iomanip>
#include <cctype>

#include <getopt.h>
#include <cctype> //std::tolower(), uppercase/lowercase conversion

#include "../sparseMEM/src/fasta.hpp"

using namespace sdsl;
using namespace std;

#define XSTR(s) STR(s)
#define STR(s)


#ifdef BWTK
typedef csa_wt<csa_wt<>::wavelet_tree_type, BWTK, 10000> tCSA_WT;
#else
typedef csa_wt<csa_wt<>::wavelet_tree_type, 4, 10000> tCSA_WT;
#endif

#ifdef LCPCOMPRESS
typedef cst_sct3<tCSA_WT, lcp_dac<> > tCST;
#else
typedef cst_sct3<tCSA_WT, lcp_dac<> > tCST;
#endif

typedef tCST::size_type size_type;
typedef vector<size_type> tVI;
typedef pair<size_type, size_type> tPII;
typedef vector<tPII> tVPII;

void usage(string prog);

enum mum_t { MUM, MAM, MEM };

tCST::size_type backward_search(const tCST::csa_type &csa,
                                const unsigned char c,
                                tCST::size_type &lb,
                                tCST::size_type &rb) {
    typedef tCST::size_type size_type;
    size_type cc = csa.char2comp[c];
    if(cc == 0) {
        lb=rb=0;
        return 0;
    }
    size_type c_before_l = csa.wavelet_tree.rank(lb, c);
    size_type c_before_r = csa.wavelet_tree.rank(rb, c);
    lb = csa.C[cc]+c_before_l;
    rb = csa.C[cc]+c_before_r;
    return rb-lb;
}

void print_maximal_exact_match(tCST::size_type p1,
                               tCST::size_type p2,
                               tCST::size_type len,
                               bool _4column,
                               const vector<string> &refdescr,
                               const vector<long> &startpos,
                               long maxdescrlen) {
    if(_4column == false) {
        cout << "- " << p1+1 << "," << p2+1 << "," << len << endl;
    }
    else {
        long refseq=0, refpos=0;
        // Use binary search to locate index of sequence and position
        // within sequence.
        vector<long>::const_iterator it = upper_bound(startpos.begin(), startpos.end(), (long)p1);
        refseq = distance(startpos.begin(), it) - 1;
        assert(refseq>=0);
        it--;
        refpos = p1 - *it;
        printf("  %s", refdescr[refseq].c_str());
        for(long j = 0; j < maxdescrlen - (long)refdescr[refseq].size() + 1; j++) {
            putchar(' ');
        }
        cout << refpos+1 << "," << p2+1 << "," << len << endl;
    }
}

struct path_item {
    tCST::size_type c, lb, rb, p;
    path_item(size_type c=0, size_type lb=0, size_type rb=0, size_type p=0):c(c),lb(lb),rb(rb),p(p) {
    }
};

void report_maximal_exact_match(const tCST &cst1,
                                const unsigned char *s2,
                                size_type c_,
                                size_type p2_,
                                size_type k,
                                bool _4column,
                                const vector<string> &refdescr,
                                const vector<long> &startpos,
                                long maxdescrlen,
                                tCST::size_type p2_offset=0) {
    size_type sa_k = cst1.csa[k];
    if(p2_ == 0) {
        print_maximal_exact_match(sa_k, p2_+p2_offset, c_, _4column, refdescr, startpos, maxdescrlen);
    } else {
        size_type lf_k = cst1.csa.lf[k];
        size_type cc2 = cst1.csa.char2comp[s2[p2_-1]];

        if( !(cst1.csa.C[cc2] <= lf_k and lf_k < cst1.csa.C[cc2+1]) ) {
            print_maximal_exact_match(sa_k, p2_+p2_offset, c_, _4column, refdescr, startpos, maxdescrlen);
        }
    }
}

//! Calculate all maximal exact matches of string s1 and s2 with minimal length l.
/*!
 *	\param cst1 	The (compressed) suffix tree of s1 supporting the parent operation.
 *      \param s2  		String s2 for which we calculate all maximal exact matches with string s1.
 *	\param n2		The length of s2.
 *	\param l		Threshold value for the length of the minimal maximal matches. \f$ l > 0 \f$
 *	\param _4column Force 4 column output
 *	\param refdescr	Description of the input sequences
 *	\param startpos Starting positions of the input sequences in the file
 *	\param maxdescrlen	Maximum length of a input sequence description
 */
void maximal_exact_matches(const tCST &cst1,
                           const unsigned char *s2,
                           tCST::size_type n2,
                           tCST::size_type l,
                           bool _4column,
                           const vector<string> &refdescr,
                           const vector<long> &startpos,
                           long maxdescrlen,
                           tCST::size_type p2_offset=0) {
    typedef tCST::node_type node_type;
    typedef tCST::size_type size_type;
    typedef list<path_item> path_type;

    assert(l > 0);
    size_type p2 = n2;
    size_type i = 0, j = cst1.size();
    size_type c = 0;

    while(p2 > 0) {
        path_type path;
        size_type lb = i, rb = j;
        backward_search(cst1.csa, s2[p2-1], lb, rb);
        while(lb != rb and p2 > 0) {
            ++c;
            if(c >= l) {
                path.push_back(path_item(c, lb, rb, p2-1) );
            }
            i = lb; j = rb;
            --p2;
            backward_search(cst1.csa, s2[p2-1], lb, rb);
        }
        for(path_type::iterator it=path.begin(); it!=path.end(); ++it) {
            size_type c_ = it->c;
            size_type lb_ = it->lb, rb_ = it->rb;
            lb = lb_; rb = lb_;
            size_type p2_ = it->p;
            while(c_ >= l) {
                for(size_type k = lb_; k < lb; ++k) {
                    report_maximal_exact_match(cst1, s2, c_, p2_, k, _4column, refdescr, startpos, maxdescrlen, p2_offset);
                }
                for(size_type k = rb; k < rb_; ++k) {
                    report_maximal_exact_match(cst1, s2, c_, p2_, k, _4column, refdescr, startpos, maxdescrlen, p2_offset);
                }
                lb = lb_; rb = rb_;
                node_type p = cst1.parent( cst1.node(lb_, rb_-1) );
                c_  = cst1.depth(p);
                lb_ = cst1.lb(p);
                rb_ = cst1.rb(p)+1;
            }
        }
        if(0==c) {
            --p2;
        } else {
            node_type p = cst1.parent( cst1.node(i, j-1) );
            c = cst1.depth(p);
            i = cst1.lb(p);
            j = cst1.rb(p)+1;
        }
    }
}

void answer_query(const string &query_fasta,
                  const tCST &cst1,
                  tCST::size_type min_len,
                  const vector<string> &refdescr,
                  const vector<long> &startpos,
                  long maxdescrlen,
                  bool rev_comp,
                  bool _4column,
                  bool nucleotides_only) {
    string meta, line;
    ifstream data(query_fasta.c_str());
    if(!data.is_open()) {
        cerr << "unable to open " << query_fasta << endl; exit(1);
    }

    // Collect meta data.
    while(!data.eof()) {
        getline(data, line); // Load one line at a time.
        if(line.length() == 0) {
            continue;
        }
        if(line[0] == '>') {
            long start = 1, end = line.length() - 1;
            trim(line, start, end);
            for(long i = start; i <= end; i++) {
                if( line[i] == ' ') break; // Behave like MUMmer 3 cut off meta after first space.
                meta += line[i];
            }
            cerr << "# " << meta << endl;
            break;
        }
    }

    string *P = new string;
    while(!data.eof()) {
        getline(data, line); // Load one line at a time.
        if(line.length() == 0) continue;
        long start = 0, end = line.length() - 1;
        // Meta tag line and start of a new sequence
        // Collect meta data.
        if(line[0] == '>') {
            if(meta != "") {
                cerr << "# P.length()=" << P->length() << endl;
                printf("> %s\n", meta.c_str());
                maximal_exact_matches(cst1, (const unsigned char*)P->c_str(), P->size(), min_len, _4column, refdescr, startpos, maxdescrlen, 0);
                if(rev_comp) {
                    reverse_complement(*P, nucleotides_only);
                    printf("> %s Reverse\n", meta.c_str());
                    maximal_exact_matches(cst1, (const unsigned char*)P->c_str(), P->size(), min_len, _4column, refdescr, startpos, maxdescrlen, 0);
                }
                delete P; P = new string;
                meta = "";
            }
            start = 1;
            trim(line, start, end);
            for(long i = start; i <= end; i++) {
                if(line[i] == ' ') break; // Behave like MUMmer 3 cut of meta after first space.
                meta += line[i];
            }
            cerr << "# " << meta << endl;
        }
        else { // Collect sequence data.
            trim(line, start,end);
            for(long i = start; i <= end; i++) {
                char c = std::tolower(line[i]);
                if(nucleotides_only) {
                    switch(c) {
                    case 'a': case 't': case 'g': case 'c': break;
                    default:
                        c = '~';
                    }
                }
                line[i] = c;
            }
            P->append(line, start, end-start+1);
        }
    }
    // Handle very last sequence.
    if(meta != "") {
        cerr << "# P.length()=" << P->length() << endl;
        printf("> %s\n", meta.c_str());
        maximal_exact_matches(cst1, (const unsigned char*)P->c_str(), P->size(), min_len, _4column, refdescr, startpos, maxdescrlen, 0);
        if(rev_comp) {
            reverse_complement(*P, nucleotides_only);
            printf("> %s Reverse\n", meta.c_str());
            maximal_exact_matches(cst1, (const unsigned char*)P->c_str(), P->size(), min_len, _4column, refdescr, startpos, maxdescrlen, 0);
        }
    }
    delete P;
}

void write_lock(int i) {
    ofstream lockfile("lock.txt", ios_base::trunc);
    lockfile<<i<<endl;
    lockfile.close();
}

int main(int argc, char* argv[]) {
    write_lock(0);
    int min_len = 20;
    //mum_t type = MEM;
    bool rev_comp = false, _4column = false, nucleotides_only = false;
    long maxdescrlen = 0;

    // Collect parameters from the command line.
    while (1) {
        static struct option long_options[] = {
            {"l", 1, 0, 0}, // 0
            //      {"mumreference", 0, 0, 0}, // 1
            {"b", 0, 0, 0}, // 2
            {"maxmatch", 0, 0, 0}, // 3
            //      {"mum", 0, 0, 0}, // 4
            //      {"mumcand", 0, 0, 0},  // 5
            {"F", 0, 0, 0}, // 6
            //      {"k", 1, 0, 0}, // 7
            //      {"threads", 1, 0, 0}, // 8
            {"n", 0, 0, 0}, // 9
            //      {"qthreads", 1, 0, 0}, // 10
            {0, 0, 0, 0}
        };
        int longindex = -1;
        int c = getopt_long_only(argc, argv, "", long_options, &longindex);
        if(c == -1) break; // Done parsing flags.
        else if(c == '?') { // If the user entered junk, let him know.
            cerr << "Invalid parameters." << endl;
            usage(argv[0]);
        }
        else {
            //mum_t type = MEM;
            // Branch on long options.
            cout << longindex << endl;
            switch(longindex) {
            case 0: min_len = atol(optarg); break;
                //      case 1: type = MAM; break;
            case 1: rev_comp = true;	break;
                //case 2: type = MEM; break;
                //      case 4: type = MUM; break;
                //      case 5: type = MAM; break;
            case 3: _4column = true; break;
                //      case 7: K = atoi(optarg); break;
                //      case 8: num_threads = atoi(optarg); break;
            case 4: nucleotides_only = true; break;
                //      case 10: query_threads = atoi(optarg) ; break;
            default: break;
            }
        }
    }
    if (argc - optind != 2) usage(argv[0]);

    string ref_fasta = argv[optind];
    string query_fasta = argv[optind+1];

    string ref;
    string query;

    vector<string> refdescr;
    vector<string> querydescr;
    vector<long> startpos;

    load_fasta(ref_fasta, ref, refdescr, startpos);
    for(size_t i=0;i<refdescr.size();i++) {
        cerr<<refdescr[i]<<endl;
        cerr<<startpos[i]<<endl;
    }
    
    // Get maximum query sequence description length.
    maxdescrlen = 0;
    for(long i = 0; i < (long)refdescr.size(); i++) {
        if(maxdescrlen < (long)refdescr[i].length())  maxdescrlen = refdescr[i].length();
    }

    if(startpos.size() > 1 ) _4column = true;

    { // free ref
        string dummy;
        dummy.swap(ref);
    }

    tCST cst;

    string file_name = ref_fasta+"_"+XSTR(BWTK)+".idx";
    std::cerr<<"# file_name "<<file_name<<std::endl;
    // create compressed suffix tree
    if(!load_from_file(cst,file_name)) {
        cerr<<"# create suffix tree of dna string "<<endl;
        startpos.clear(); refdescr.clear(); 
        cerr<<"ref_fasta = "<<ref_fasta<<endl;
        load_fasta(ref_fasta, ref, refdescr, startpos);
        construct_im(cst, ref, 1);
        cerr<<"# suffix tree created"<<endl;
        store_to_file(cst, file_name);
    } else {
        cerr<<"#load index from disk"<<endl;
    }
    cerr<<"# size of suffix tree in MB "<< (size_in_mega_bytes(cst))/(1<<20) <<endl;
    cerr<<"# alphabet size: "<< cst.csa.sigma - 1 << endl;
    for(int i = 1;i<cst.csa.sigma; ++i) {
        cerr<<"\""<<cst.csa.comp2char[i]<<"\" ";
    }
    cerr<<endl;
    {
        string dummy;
        dummy.swap(ref);
    }
    write_lock(1);
    stop_watch stopwatch;
    stopwatch.start();
    answer_query(query_fasta, cst, min_len, refdescr, startpos, maxdescrlen, rev_comp, _4column, nucleotides_only);
    stopwatch.stop();
    cerr<<"# Time for algorithm: ";
    cerr<<stopwatch.getUserTime();
    ofstream stopwatchfile("stopwatch.txt", ios_base::trunc);
    stopwatchfile<<stopwatch.getUserTime()<<endl;
    stopwatchfile.close();
    write_lock(0);
    cerr<<endl;
}

void usage(string prog) {
    cerr << "Usage: " << prog << " [options] <reference-file> <query-file>" << endl;
    cerr << "Implemented MUMmer v3 options:" << endl;
    //  cerr << "-mum           compute maximal matches that are unique in both sequences" << endl;
    //  cerr << "-mumreference  compute maximal matches that are unique in the reference-" << endl;
    //  cerr << "               sequence but not necessarily in the query-sequence (default)" << endl;
    //  cerr << "-mumcand       same as -mumreference" << endl;
    cerr << "-maxmatch      compute all maximal matches regardless of their uniqueness" << endl;
    cerr << "-l             set the minimum length of a match" << endl;
    cerr << "               if not set, the default value is 20" << endl;
    cerr << "-b             compute forward and reverse complement matches" << endl;
    cerr << "-F             force 4 column output format regardless of the number of" << endl;
    cerr << "               reference sequence inputs"  << endl;
    cerr << "-n             match only the characters a, c, g, or t" << endl;
    //  cerr << endl;
    //  cerr << "Additional options:" << endl;
    //  cerr << "-k             sampled suffix positions (one by default)" << endl;
    //  cerr << "-threads       number of threads to use for -maxmatch, only valid k > 1 " << endl;
    //  cerr << "-qthreads      number of threads to use for queries " << endl;
    //  cerr << endl;
    //  cerr << "Example usage:" << endl;
    //  cerr << endl;
    //  cerr << "./mummer -maxmatch -l 20 -b -n -k 3 -threads 3 query.fa ref.fa" << endl;
    //  cerr << "Find all maximal matches on forward and reverse strands" << endl;
    //  cerr << "of length 20 or greater, matching only a, c, t, or g." << endl;
    //  cerr << "Index every 3rd position in the ref.fa and use 3 threads to find MEMs." << endl;
    //  cerr << "Fastest method for one long query sequence." << endl;
    //  cerr << endl;
    //  cerr << "./mummer -maxmatch -l 20 -b -n -k 3 -qthreads 3 query.fa ref.fa" << endl;
    //  cerr << "Same as above, but now use a single thread for every query sequence in" << endl;
    //  cerr << "query.fa. Fastest for many small query sequences." << endl;
    exit(1);
}
