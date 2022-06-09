#include <iostream>

#include "lest.hpp"
#include "Kmer.hpp"
#include "exampleHash.hpp"
#include "Minimiser.hpp"

using namespace std;
using namespace lest;

const test minimiser[] {
    CASE({"alphabetical order hash (A < C < T < G) : ") {

            Kmer testKmer1 = Kmer(0b0100101010, 5); // CATTT
            Kmer testKmer2 = Kmer(0b11111001011111, 7); // GGTCCGG
            Kmer testKmer3 = Kmer(0b101001000010000110, 9); // TTCAATACT
            Kmer testKmer4 = Kmer(0b0010000100011111011000, 11); // ATACACGGCTA
            Kmer testKmer5 = Kmer(0b00111000100000100110010011, 13); // AGTATAATCTCAG
            Kmer testKmer6 = Kmer(0b000011001001110000111001101111, 15); // AAGATCGAAGTCTGG
            Kmer testKmer7 = Kmer(0b1011101101100101, 8);// TGTGCTCC
            Kmer testKmer8 = Kmer(0b01001101101000101001, 10); // CAGCTTATTC
            Kmer testKmer9 = Kmer(0b000011000000110000010010, 12); // AAGAAAGAACAT
            Kmer testKmer10 = Kmer(0b1101010011011101110101111000, 14); // GCCAGCGCGCCGTA

            hashPos res1 = alpha(testKmer1, 3);
            hashPos expected1 = {0b001010 , 1}; //ATT
            hashPos res2 = alpha(testKmer2, 4);
            hashPos expected2 = {0b01011111 , 3}; // CCGG
            hashPos res3 = alpha(testKmer3, 5);
            hashPos expected3 = {0b0000100001 , 3}; // AATAC
            hashPos res4 = alpha(testKmer4, 6);
            hashPos expected4 = {0b000100011111 , 2}; // ACACGG
            hashPos res5 = alpha(testKmer5, 7);
            hashPos expected5 = {0b001001100100 , 5}; // AATCTCA
            hashPos res6 = alpha(testKmer6, 8);
            hashPos expected6 = {0b0000110010011100 , 0}; // AAGATCGA
            hashPos res7 = alpha(testKmer7, 4);
            hashPos expected7 = {0b01100101 , 4}; // CTCC
            hashPos res8 = alpha(testKmer8, 5);
            hashPos expected8 = {0b0011011010 , 1}; // AGCTT
            hashPos res9 = alpha(testKmer9, 6);
            hashPos expected9 = {0b000000110000 , 3}; // AAAGAA
            hashPos res10 = alpha(testKmer10, 7);
            hashPos expected10 = {0b00110111011101 , 3}; // AGCGCGC

            EXPECT(res1.hashValue == expected1.hashValue);
            EXPECT(res1.pos == expected1.pos);
            EXPECT(res2.hashValue == expected2.hashValue);
            EXPECT(res2.pos == expected2.pos);
            EXPECT(res3.hashValue == expected3.hashValue);
            EXPECT(res3.pos == expected3.pos);
            EXPECT(res4.hashValue == expected4.hashValue);
            EXPECT(res4.pos == expected4.pos);
            EXPECT(res5.hashValue == expected5.hashValue);
            EXPECT(res5.pos == expected5.pos);
            EXPECT(res6.hashValue == expected6.hashValue);
            EXPECT(res6.pos == expected6.pos);
            EXPECT(res7.hashValue == expected7.hashValue);
            EXPECT(res7.pos == expected7.pos);
            EXPECT(res8.hashValue == expected8.hashValue);
            EXPECT(res8.pos == expected8.pos);
            EXPECT(res9.hashValue == expected9.hashValue);
            EXPECT(res9.pos == expected9.pos);
            EXPECT(res10.hashValue == expected10.hashValue);
            EXPECT(res10.pos == expected10.pos);

        }},
         CASE("")
};

int main() {
    int failures(0);
    cout << "*** Test on alphabetical order ***" << endl;
    if ((failures = run(minimiser))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    return 0;
}