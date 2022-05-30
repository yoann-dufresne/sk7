#include <iostream>
#include <unordered_map>

#include "headers/FastaReader.hpp"
#include "headers/Kmer.hpp"
#include "headers/exampleHash.hpp"
#include "headers/Minimiser.hpp"

using namespace std;

int main(int argc, char*argv[]) {
    if (argc < 4) {
        cout << "Usage : ./sk7 file k m" << endl;
        return 0;
    }

    int k = atoi(argv[2]);
    int m = atoi(argv[3]);

    FILE* fasta = fopen(argv[1], "r");
    unordered_map<char, int> encoding;
    encoding['A'] = 0;
    encoding['C'] = 1;
    encoding['G'] = 2;
    encoding['T'] = 3;

    FastaReader fastareader = FastaReader(fasta, false);
    uint64_t filter = fastareader.set_filter(k);

    uint64_t current_value = 0;
    int cmpt  = 0;
    Minimiser minimiser(alpha, m);

    while (true) {
        char c = fastareader.read_fasta();
        if (c == EOF) break;
        cmpt++;
        current_value = ((current_value << 2) + encoding.at(c)) & filter;
        if (cmpt == k) {
            Kmer kmer = Kmer(current_value, k);
            minimiser.init(kmer);
            cout << "First kmer : " << kmer.toString() << " of value : " << kmer.getValue()
            << " its minimiser is : " << minimiser.toString()
            << " its pos is : " << minimiser.getPos()
            << " and value : " << minimiser.getValue() << endl;
        }
        if (cmpt > k) {
            Kmer kmer = Kmer(current_value, k);
            minimiser.fromNewEnd(kmer);
            cout << "kmer : " << kmer.toString() << " its value is : " << kmer.getValue()
            << " its minimiser is : " << minimiser.toString()
            << " at pos : " << minimiser.getPos()
            <<" of value : " << minimiser.getValue() << endl;
        }
    }

    fclose(fasta);

    return 0;
}
