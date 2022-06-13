#include <iostream>

#include "lest.hpp"
#include "Kmer.hpp"
#include "exampleHash.hpp"
#include "Minimiser.hpp"
#include "utils.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-braces"

using namespace std;
using namespace lest;

const test minimiser[] {
    CASE("alphabetical order hash (A < C < T < G) : ") {

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

        },
         CASE("Minimiser from kmer: ") {
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

            Minimiser minimiser1 = Minimiser(alpha, 3, testKmer1); // ATT
            Minimiser minimiser2 = Minimiser(alpha, 4, testKmer2); // CCGG
            Minimiser minimiser3 = Minimiser(alpha, 5, testKmer3); // AATAC
            Minimiser minimiser4 = Minimiser(alpha, 6, testKmer4); // ACACGG
            Minimiser minimiser5 = Minimiser(alpha, 7, testKmer5); // AATCTCA
            Minimiser minimiser6 = Minimiser(alpha, 8, testKmer6); // AAGATCGA
            Minimiser minimiser7 = Minimiser(alpha, 4, testKmer7); // CTCC
            Minimiser minimiser8 = Minimiser(alpha, 5, testKmer8); // AGCTT
            Minimiser minimiser9 = Minimiser(alpha, 6, testKmer9); // AAAGAA
            Minimiser minimiser10 = Minimiser(alpha, 7, testKmer10); // AGCGCGC


            EXPECT(minimiser1.getValue() == (uint64_t) 0b001010);
            EXPECT(minimiser2.getValue() == (uint64_t) 0b01011111);
            EXPECT(minimiser3.getValue() == (uint64_t) 0b0000100001);
            EXPECT(minimiser4.getValue() == (uint64_t) 0b000100011111);
            EXPECT(minimiser5.getValue() == (uint64_t) 0b00001001100100);
            EXPECT(minimiser6.getValue() == (uint64_t) 0b0000110010011100);
            EXPECT(minimiser7.getValue() == (uint64_t) 0b01100101);
            EXPECT(minimiser8.getValue() == (uint64_t) 0b0011011010);
            EXPECT(minimiser9.getValue() == (uint64_t) 0b000000110000);
            EXPECT(minimiser10.getValue() == (uint64_t) 0b00110111011101);


         },
         CASE("Minimiser from previous one : ") {
            Kmer testKmer1 = Kmer(0b0100101010, 5); // CATTT
            Kmer testKmer2 = Kmer(0b0010101011, 5); // ATTTG
            Kmer testKmer3 = Kmer(0b1010101101, 5); // TTTGC
            Kmer testKmer4 = Kmer(0b1010110100, 5); // TTGCA
            Kmer testKmer5 = Kmer(0b1011010010, 5); // TGCAT

            Minimiser minimiser = Minimiser(alpha, 3);
            minimiser.init(testKmer1); //ATT
            EXPECT(minimiser.getValue() == (uint64_t) 0b001010);
            minimiser.fromNewEnd(testKmer2); // ATT
            EXPECT(minimiser.getValue() == (uint64_t) 0b001010);
            minimiser.fromNewEnd(testKmer3); // TTT
            EXPECT(minimiser.getValue() == (uint64_t) 0b101010);
            minimiser.fromNewEnd(testKmer4); // TTG
            EXPECT(minimiser.getValue() == (uint64_t) 0b101011);
            minimiser.fromNewEnd(testKmer5); // CAT
            EXPECT(minimiser.getValue() == (uint64_t) 0b010010);

    }
};

const test kmerTest[] {
    CASE("Get Subkmer") {
        Kmer testKmer1 = Kmer(0b0100101010, 5); // CATTT
        Kmer testKmer2 = Kmer(0b11111001011111, 7); // GGTCCGG
        Kmer testKmer3 = Kmer(0b101001000010000110, 9); // TTCAATACT
        Kmer testKmer4 = Kmer(0b0010000100011111011000, 11); // ATACACGGCTA

        EXPECT(testKmer1.getSubKmer(1, 3).getValue() == (uint64_t) 0b001010);
        EXPECT(testKmer1.getSubKmer(1, 3).getLength() == 3);
        EXPECT(testKmer2.getSubKmer(4, 6).getValue() == (uint64_t) 0b011111);
        EXPECT(testKmer2.getSubKmer(4, 6).getLength() == 3);
        EXPECT(testKmer3.getSubKmer(2, 7).getValue() == (uint64_t) 0b010000100001);
        EXPECT(testKmer3.getSubKmer(2, 7).getLength() == 6);
        EXPECT(testKmer4.getSubKmer(6, 8).getValue() == (uint64_t) 0b111101);
        EXPECT(testKmer4.getSubKmer(6, 8).getLength() == 3);

        },
        CASE("Remove Part") {
        Kmer testKmer1 = Kmer(0b0100101010, 5); // CATTT
        Kmer testKmer2 = Kmer(0b11111001011111, 7); // GGTCCGG
        Kmer testKmer3 = Kmer(0b101001000010000110, 9); // TTCAATACT
        Kmer testKmer4 = Kmer(0b0010000100011111011000, 11); // ATACACGGCTA
        EXPECT(testKmer1.removePart(1, 3).getValue() == (uint64_t) 0b0110);
        EXPECT(testKmer1.removePart(1, 3).getLength() == 2);
        EXPECT(testKmer2.removePart(4, 2).getValue() == (uint64_t) 0b1111100111);
        EXPECT(testKmer2.removePart(4, 2).getLength() == 5);
        EXPECT(testKmer3.removePart(2, 4).getValue() == (uint64_t) 0b1010000110);
        EXPECT(testKmer3.removePart(2, 4).getLength() == 5);
        EXPECT(testKmer4.removePart(0, 5).getValue() == (uint64_t) 0b011111011000);
        EXPECT(testKmer4.removePart(0, 5).getLength() == 6);
    }
};

const test utilsTest[] {
    CASE("InterleaveOrder") {
        Kmer testKmer1 = Kmer(0b0100101010, 5); // CATTT
        Kmer testKmer2 = Kmer(0b11111001011111, 7); // GGTCCGG
        Kmer testKmer3 = Kmer(0b101001000010000110, 9); // TTCAATACT
        Kmer testKmer4 = Kmer(0b0010000100011111011000, 11); // ATACACGGCTA
        Minimiser minimiser1 = Minimiser(alpha, 3, testKmer1);
        Minimiser minimiser2 = Minimiser(alpha, 4, testKmer2);
        Minimiser minimiser3 = Minimiser(alpha, 5, testKmer3);
        Minimiser minimiser4 = Minimiser(alpha, 6, testKmer4);
        Kmer withoutMinimiser1 = testKmer1.removePart(minimiser1.getPos(), 3);
        Kmer withoutMinimiser2 = testKmer2.removePart(minimiser2.getPos(), 4);
        Kmer withoutMinimiser3 = testKmer3.removePart(minimiser3.getPos(), 5);
        Kmer withoutMinimiser4 = testKmer4.removePart(minimiser4.getPos(), 6);
        uint64_t mask1 = interleavedOrder(withoutMinimiser1, minimiser1.getPos());
        uint64_t mask2 = interleavedOrder(withoutMinimiser2, minimiser2.getPos());
        uint64_t mask3 = interleavedOrder(withoutMinimiser3, minimiser2.getPos());
        uint64_t mask4 = interleavedOrder(withoutMinimiser4, minimiser4.getPos());

        EXPECT(withoutMinimiser1.getValue() == (uint64_t) 0b1001);
        EXPECT(withoutMinimiser2.getValue() == (uint64_t) 0b001000110011);
        EXPECT(withoutMinimiser3.getValue() == (uint64_t) 0b100100100010);
        EXPECT(withoutMinimiser4.getValue() == (uint64_t) 0b011010000000);

        }
};

int main() {
    int failures(0);
    cout << "*** Test on minimiser order ***" << endl;
    if ((failures = run(minimiser))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    cout << "*** Test on Kmer methods ***" << endl;
    if ((failures = run(kmerTest))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    return 0;
}
#pragma clang diagnostic pop