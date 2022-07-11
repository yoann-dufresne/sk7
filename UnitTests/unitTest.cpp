#include <iostream>

#include "lest.hpp"
#include "Bucket.hpp"

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

    CASE("Minimiser with reverse complement : ") {
        sk7::initLib(8, 3);
        Kmer testKmer1 = Kmer(0b0100101101110011); // CATGCGAG
        Kmer testKmer2 = Kmer(0b0111001111100111); // CGAGGTCG
        Kmer testKmer3 = Kmer(0b0111011101110111); // CGCGCGCG
        Kmer testKmer4 = Kmer(0b0001110001111110); // ACGACGGT
        Kmer testKmer5 = Kmer(0b0100011111011110); // CACGGCGT

        Minimiser minimiser1 = Minimiser(alpha, testKmer1);
        Minimiser minimiser2 = Minimiser(alpha, testKmer2);
        Minimiser minimiser3 = Minimiser(alpha, testKmer3);
        Minimiser minimiser4 = Minimiser(alpha, testKmer4);
        Minimiser minimiser5 = Minimiser(alpha, testKmer5);

        EXPECT(minimiser1.getValue() == (uint64_t) 0b001011);
        EXPECT(minimiser1.getPos() == 1);
        EXPECT(testKmer1 == Kmer(0b0100101101110011));

        EXPECT(minimiser2.getValue() == (uint64_t) 0b000101);
        EXPECT(minimiser2.getPos() == 2);
        EXPECT(testKmer2 == Kmer(0b0111000101100111));

        EXPECT(minimiser3.getValue() == (uint64_t) 0b011101);
        EXPECT(minimiser3.getPos() == 4);
        EXPECT(testKmer3 == Kmer(0b0111011101110111));

        EXPECT(minimiser4.getValue() == (uint64_t) 0b000101);
        EXPECT(minimiser4.getPos() == 0);
        EXPECT(testKmer4 == Kmer(0b0001011110011110));

        EXPECT(minimiser5.getValue() == (uint64_t) 0b000111);
        EXPECT(minimiser5.getPos() == 1);
        EXPECT(testKmer5 == Kmer(0b0100011111011110));

        sk7::initLib(5, 3);
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
        EXPECT(mask1 == (uint64_t) 0b1111);
        EXPECT(mask2 == (uint64_t) 0b001100110011);
        EXPECT(mask3 == (uint64_t) 0b111100110011);
        EXPECT(mask4 == (uint64_t) 0b111111111100);

        },
        CASE ("reorderValue") {
        EXPECT(reorderValue(0b01101100, 1,2) == (uint64_t) 0b10000111);
        EXPECT(reorderValue(0b01101101, 2,2) == (uint64_t) 0b01100111);
        EXPECT(reorderValue(0b01110010, 2,1) == (uint64_t) 0b10110100);
        EXPECT(reorderValue(0b0110, 1,1) == (uint64_t) 0b1001);
    }
};

const test superKmerTest[] {
     CASE("accessBits") {
        vector<uint8_t> testVector = {(uint8_t) 0b10001011, (uint8_t) 0b11001101, (uint8_t) 0b01001100};
        SuperKmer skTest(testVector);
        EXPECT(skTest.accessBits(0, 24) == (uint64_t) 0b100010111100110101001100);
        EXPECT(skTest.accessBits(8, 16) == (uint64_t) 0b11001101);
        EXPECT(skTest.accessBits(5, 10) == (uint64_t) 0b01111);
        EXPECT(skTest.accessBits(14, 20) == (uint64_t) 0b010100);
     },
     CASE("setBits") {
        SuperKmer SKTest = SuperKmer({0});
        SKTest.setBits(0, 2, 0b11);
        SKTest.setBits(2, 2, 0b10);
        SKTest.setBits(4, 4, 0b1001);
        SKTest.setBits(8, 8, 0b10010110);
        EXPECT(SKTest.accessBits(0, 16) == (uint64_t) 0b1110100110010110);
     },

     CASE("Intersection") {

        SuperKmer SK1 = SuperKmer({0b10010111, 0b00100000});
        SuperKmer SK2 = SuperKmer({0b01100111, 0b00000000});
        bool equal = (SK1 & SK2) == SuperKmer({0b01010111});
        EXPECT(equal);

        SuperKmer SK3 = SuperKmer({0b01100110, 0b11000000});
        SuperKmer SK4 = SuperKmer({0b10010100, 0b00100000});
        equal = (SK3 & SK4) == SuperKmer();
        EXPECT(equal);

        SuperKmer SK5 = SuperKmer({0b10011000, 0b00100000});
        SuperKmer SK6 = SuperKmer();
        equal = (SK5 & SK6) == SuperKmer();
        EXPECT(equal);

        SuperKmer SK7 = SuperKmer({0b10100110, 0b10100000});
        SuperKmer SK8 = SuperKmer({0b10101110, 0b10100000});
        equal = (SK7 & SK8) == SuperKmer({0b10000010, 0b00100000});
        EXPECT(equal);

        SuperKmer SK9 = SuperKmer({0b10000010, 0b00110000});
        SuperKmer SK10 = SuperKmer({0b00100100, 0b00000000});
        equal = (SK9 & SK10) == SuperKmer();
        EXPECT(equal);

        SuperKmer SK11 = SuperKmer({0b10100111, 0b11100000});
        SuperKmer SK12 = SuperKmer({0b10100111, 0b1101000});
        equal = (SK11 & SK12) == SuperKmer({0b10010111, 0b00100000});
        EXPECT(equal);

        SuperKmer SK13 = SuperKmer({0b01100111,  0b10000000});
        SuperKmer SK14 = SuperKmer({0b01010110});
        equal = (SK13 & SK14) == SuperKmer();
        EXPECT(equal);
     },
     CASE("Union") {
        SuperKmer SK1 = SuperKmer({0b01101011, 0b10000000});
        SuperKmer SK2 = SuperKmer({0b10011011, 0b00010000});
        SuperKmer SK3 = SuperKmer({0b10000011, 0b00010000});

        EXPECT((SK1 | SK2) == SuperKmer({0b10101011, 0b10010000}));
        EXPECT((SK1 | SK3) == SuperKmer({0b10101011, 0b10010000}));

        SuperKmer SK4 = SuperKmer({0b00101000, 0b10000000});
        SuperKmer SK5 = SuperKmer({0b01011011});

        EXPECT((SK4 | SK5) == SuperKmer({0b01101011, 0b10000000}));
     },
     CASE("split") {
        SuperKmer SK1 = SuperKmer({0b10101011, 0b10010000});
        SuperKmer SK2 = SuperKmer({0b01010110});
        SuperKmer SK3 = SuperKmer({0b01101011, 0b10000000});
        SuperKmer SK4 = SuperKmer({0b10000010, 0b00100000});
        SuperKmer SK5 = SuperKmer({0b00101100, 0b00000000});

        std::vector<SuperKmer> split1 = SK1.split();
        EXPECT(split1.size() == (uint64_t) 3);
        bool equal = split1.at(0) == SuperKmer({0b10000011, 0b00010000});
        EXPECT(equal);
        equal = split1.at(1) == SuperKmer({0b01011011}) ;
        EXPECT(equal);
        equal = split1.at(2) == SuperKmer({0b00101000, 0b10000000});
        EXPECT(equal);

        std::vector<SuperKmer> split2 = SK2.split();
        EXPECT(split2.size() == (uint64_t) 3);
        equal = split2.at(1) == SK2;
        EXPECT(equal);

        std::vector<SuperKmer> split3 = SK3.split();
        EXPECT(split3.size() == (uint64_t) 3);
        equal = split3.at(1) == SuperKmer({0b01011011});
        EXPECT(equal);
        equal = split3.at(2) == SuperKmer({0b00101000, 0b10000000});
        EXPECT(equal);

        std::vector<SuperKmer> split4 = SK4.split();
        EXPECT(split4.size() == (uint64_t) 3);
        equal = split4.at(0) == SK4;
        EXPECT(equal);

        std::vector<SuperKmer> split5 = SK5.split();
        EXPECT(split5.size() == (uint64_t) 3);
        equal = split5.at(2) == SK5;
        EXPECT(equal);
     },
     CASE("nonInterleaved value") {
        SuperKmer SK1 = SuperKmer({0b10101011, 0b10010000});
        std::vector<SuperKmer> split = SK1.split();
        EXPECT(split.at(0).nonInterleavedKmerValue() == (uint64_t) 0b0111000000);
        EXPECT(split.at(1).nonInterleavedKmerValue() == (uint64_t) 0b1100000010);
        EXPECT(split.at(2).nonInterleavedKmerValue() == (uint64_t) 0b0000001010);

     },
     CASE("CompareSK") {
        SuperKmer SK1 = SuperKmer({0b10101011, 0b10010000});
        SuperKmer SK2 = SuperKmer({0b01010110});
        SuperKmer SK4 = SuperKmer({0b01101011, 0b10000000});
        SuperKmer SK6 = SuperKmer({0b10000010, 0b00100000});
        SuperKmer SK7 = SuperKmer({0b00101100, 0b00000000});


        std::vector<SuperKmer::logic> expected = {SuperKmer::EQUAL, SuperKmer::EQUAL, SuperKmer::EQUAL};
        EXPECT(SuperKmer::compareSK(SK1, SK1) == expected);

        expected = {SuperKmer::INCOMPARABLE, SuperKmer::SUPERIOR, SuperKmer::INCOMPARABLE};
        EXPECT(SuperKmer::compareSK(SK1, SK2) == expected);

        expected = {SuperKmer::INCOMPARABLE, SuperKmer::EQUAL, SuperKmer::EQUAL};
        EXPECT(SuperKmer::compareSK(SK1, SK4) == expected);

        expected = {SuperKmer::SUPERIOR, SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE};
        EXPECT(SuperKmer::compareSK(SK1, SK6) == expected);

        expected = {SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE, SuperKmer::INFERIOR};
        EXPECT(SuperKmer::compareSK(SK1, SK7) == expected);

     },

     CASE("readKmer") {
         SuperKmer SK1 = SuperKmer({0b10101011, 0b10010000});
         SuperKmer SK2 = SuperKmer({0b01011000});
         SuperKmer SK3 = SuperKmer({0b10000011, 0b00110000});
         SuperKmer SK4 = SuperKmer({0b00101000, 0b00000000});

         EXPECT(SK1.readKmer(0).getValue() == (uint64_t) 0b10001000);
         EXPECT(SK1.readKmer(1).getValue() == (uint64_t) 0b1011);
         EXPECT(SK1.readKmer(2).getValue() == (uint64_t) 0b110001);

         EXPECT(SK2.readKmer(0).length == 0);
         EXPECT(SK2.readKmer(1).getValue() == (uint64_t) 0b1000);
         EXPECT(SK2.readKmer(2).length == 0);

         EXPECT(SK3.readKmer(0).length == 0);
         EXPECT(SK3.readKmer(1).length == 0);
         EXPECT(SK3.readKmer(2).getValue() == (uint64_t) 0b110011);

         EXPECT(SK4.readKmer(0).getValue() == (uint64_t) 0b10000000);
         EXPECT(SK4.readKmer(1).length == 0);
         EXPECT(SK4.readKmer(2).length == 0);
     }
};

const test bucketTest[] {
    CASE("find") {
        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));

        // The bucket contains : G|CT, CG|TT, GG|

        Kmer toSearch = Kmer(0b0000001000, 5); // AAATA
        Kmer toSearch2 = Kmer(0b1000000001, 5); // TAAAC
        Kmer toSearch3 = Kmer(0b1001000000, 5); // TCAAA
        Kmer toSearch4 = Kmer(0b0000000110, 5); // AAACT
        Kmer toSearch5 = Kmer(0b1100000010, 5) ; // GAAAT
        Kmer toSearch6 = Kmer(0b1111000000, 5); // GGAAA
        Kmer toSearch7 = Kmer(0b0111000000, 5); // CGAAA
        Kmer toSearch8 = Kmer(0b0000001010, 5); // AAATT
        Kmer toSearch9 = Kmer(0b1100000001, 5); // GAAAC

        int position;

        EXPECT(not bucket.find(toSearch, position)); EXPECT(position == 1);
        EXPECT(not bucket.find(toSearch2, position)); EXPECT(position == 0);
        EXPECT(not bucket.find(toSearch3, position)); EXPECT(position == 0);
        EXPECT(bucket.find(toSearch4, position)); EXPECT(position == 0);
        EXPECT(bucket.find(toSearch5, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch6, position)); EXPECT(position == 2);
        EXPECT(bucket.find(toSearch7, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch8, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch9, position)); EXPECT(position == 0);

    },
    CASE("isSorted") {

        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));

        Bucket spaced = Bucket(0);
        spaced.addToList(SuperKmer({0b10000010, 0b00100000}));
        spaced.addToList(SuperKmer({0b01011001}));
        spaced.addToList(SuperKmer({0b10000001, 0b00010000}));

        EXPECT(bucket.isSorted());
        EXPECT(not spaced.isSorted());
    }
        ,
    CASE("addKmer") {
        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));

        Kmer toAdd = Kmer(0b0000001000, 5); // AAATA
        Kmer toAdd2 = Kmer(0b1000000001, 5); // TAAAC
        Kmer toAdd3 = Kmer(0b1001000000, 5); // TCAAA
        Kmer toSearch4 = Kmer(0b0000000110, 5); // AAACT;
        Kmer toSearch5 = Kmer(0b1100000010, 5) ; // GAAAT
        Kmer toSearch6 = Kmer(0b1111000000, 5); // GGAAA
        Kmer toSearch7 = Kmer(0b0111000000, 5) ;// CGAAA
        Kmer toSearch8 = Kmer(0b0000001010, 5); // AAATT
        Kmer toSearch9 = Kmer(0b1100000001, 5); // GAAAC

        EXPECT(bucket.getListSize() == (uint64_t) 3);
        bucket.addKmer(toAdd);
        EXPECT(bucket.getListSize() == (uint64_t) 4);
        bucket.addKmer(toAdd2);
        EXPECT(bucket.getListSize() == (uint64_t) 5);
        bucket.addKmer(toAdd3);
        EXPECT(bucket.getListSize() == (uint64_t) 6);
        bucket.addKmer(toAdd);
        EXPECT(bucket.getListSize() == (uint64_t) 6);


        EXPECT(bucket.isSorted());

        int position;
        EXPECT(bucket.find(toAdd, position) ); EXPECT(position == 3);
        EXPECT(bucket.find(toAdd2, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toAdd3, position)); EXPECT(position == 0);
        EXPECT(bucket.find(toSearch4, position)); EXPECT(position == 2);
        EXPECT(bucket.find(toSearch5, position)); EXPECT(position == 4);
        EXPECT(bucket.find(toSearch6, position)); EXPECT(position == 5);
        EXPECT(bucket.find(toSearch7, position)); EXPECT(position == 4);
        EXPECT(bucket.find(toSearch8, position)); EXPECT(position == 4);
        EXPECT(bucket.find(toSearch9, position)); EXPECT(position == 2);
    },
    CASE("SKtoKmer") {
        Bucket bucket = Bucket(0);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b01101011, 0b01000000})).getValue() == (uint64_t) 0b110000001001);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b11110110, 0b10110010})).getValue() == (uint64_t) 0b101110000000011000);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b00101100, 0b00000000})).getValue() == (uint64_t) 0b0000001100);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b10000010, 0b00100000})).getValue() == (uint64_t) 0b1010000000);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b01010111, 0b00000000})).getValue() == (uint64_t) 0b1100000001);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b10101011, 0b10010000})).getValue() == (uint64_t) 0b01110000001010);
        EXPECT(bucket.SKtoKmer(SuperKmer({0b01011001})).getValue() == (uint64_t) 0b0100000010);
    },

    CASE("addSuperKmer") {
        /// SuperKmer 1 : G|CT
        uint8_t firstSection1 = 0b00100111; //1, 2, C, G
        uint8_t secondSection1 = 0b10000000; //T, _, _, _
        vector<TYPE> tab1 = vector<TYPE>();
        tab1.push_back(firstSection1);
        tab1.push_back(secondSection1);
        SuperKmer SK1(tab1);

        ///SuperKmer 2 : CG|TT
        uint8_t firstSection2 = 0b10101011; //2, 2, T, G
        uint8_t secondSection2 = 0b10010000; // T, C, _, _
        vector<TYPE> tab2 = vector<TYPE>();
        tab2.push_back(firstSection2);
        tab2.push_back(secondSection2);
        SuperKmer SK2(tab2);

        ///SuperKmer 3 : GG|
        uint8_t firstSection3 = 0b10000011; //2, 0, _, G
        uint8_t secondSection3 = 0b00110000;//_, G, _, _
        vector<TYPE> tab3 = vector<TYPE>();
        tab3.push_back(firstSection3);
        tab3.push_back(secondSection3);
        SuperKmer SK3(tab3);

        Bucket bucket(0);
        bucket.addSuperKmer(SK3);
        EXPECT(bucket.getListSize() == (uint64_t) 1);
        bucket.addSuperKmer(SK2);
        EXPECT(bucket.getListSize() == (uint64_t) 4);
        bucket.addSuperKmer(SK2);
        EXPECT(bucket.getListSize() == (uint64_t) 4);
        bucket.addSuperKmer(SK1);
        EXPECT(bucket.getListSize() == (uint64_t) 5);
        bucket.addSuperKmer(SK1);
        EXPECT(bucket.getListSize() == (uint64_t) 5);
        bucket.addSuperKmer(SK2);
        EXPECT(bucket.getListSize() == (uint64_t) 5);
        bucket.addSuperKmer(SK3);
        EXPECT(bucket.getListSize() == (uint64_t) 5);

        EXPECT(bucket.isSorted());
    },



    CASE("findNextOkPosition") {
        Bucket bucket = Bucket(0);
        bucket.addSuperKmer(SuperKmer({0b01010110}));
        bucket.addSuperKmer(SuperKmer({0b00101100, 0b01000000}));
        bucket.addSuperKmer(SuperKmer({0b01100111, 0b10000000}));
        bucket.addSuperKmer(SuperKmer({0b10000011, 0b00110000}));

        EXPECT(bucket.findNextOkPosition(SuperKmer({0b01011101}), bucket.getListCopy(), 0) == (uint64_t) 4);
        EXPECT(bucket.findNextOkPosition(SuperKmer({0b10000001, 0b00100000}), bucket.getListCopy(), 0) == (uint64_t) 1);
        EXPECT(bucket.findNextOkPosition(SuperKmer({0b00101100, 0b11000000}), bucket.getListCopy(), 1) == (uint64_t) 5);

    },
    CASE("Intersection") {
        Bucket bucket1 = Bucket(0);
        Bucket bucket2 = Bucket(0);
        Bucket bucket3 = Bucket(0);
        Bucket bucket4 = Bucket(0);
        Bucket bucket5 = Bucket(0);

        int position;

        bucket1.addToList(SuperKmer({0b00101000, 0b00000000}));
        bucket1.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket1.addToList(SuperKmer({0b10000011, 0b00110000}));
//        cout << "------ 1 ------" << endl;
//        bucket1.print();
        EXPECT(bucket1.isSorted());

        bucket2.addToList(SuperKmer({0b10100111, 0b10010000}));
        bucket2.addToList(SuperKmer({0b01011101}));
        bucket2.addToList(SuperKmer({0b00101100, 0b11000000}));
//        cout << "----- 2 ------" << endl;
//        bucket2.print();
        EXPECT(bucket2.isSorted());

//        cout << "---- 1 & 2 ----" << endl;
//        (bucket1 & bucket2).print();
        EXPECT((bucket1 & bucket2).getListSize() == (uint64_t) 1);
        EXPECT((bucket1 & bucket2).find(Kmer(0b0111000000), position));


        bucket3.addToList(SuperKmer({0b01010110}));
        bucket3.addToList(SuperKmer({0b10101001, 0b10010000}));
        bucket3.addToList(SuperKmer({0b10101010, 0b11100000}));
        bucket3.addSuperKmer(SuperKmer({0b01011011}));
        bucket3.addSuperKmer(SuperKmer({0b00101100, 0b01000000}));
//        cout << "------ 3 ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 & 3 --------" << endl;
//        (bucket1 & bucket3).print();
        EXPECT((bucket1 & bucket3).getListSize() == (uint64_t) 2);
        EXPECT((bucket1 & bucket3).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 & bucket3).find(Kmer(0b0000001010), position));


        bucket3.addSuperKmer(SuperKmer({0b10000011, 0b00110000}));
//        cout << "------ 3++ ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 & 3++ --------" << endl;
//        (bucket1 & bucket3).print();
        EXPECT((bucket1 & bucket3).getListSize() == (uint64_t) 3);
        EXPECT((bucket1 & bucket3).isSorted());
        EXPECT((bucket1 & bucket3).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 & bucket3).find(Kmer(0b0000001010), position));
        EXPECT((bucket1 & bucket3).find(Kmer(0b1111000000), position));

        bucket4.addToList(SuperKmer({0b10000001, 0b00000000}));
        bucket4.addToList(SuperKmer({0b00101000, 0b00000000}));
        bucket4.addToList(SuperKmer({0b01101001, 0b10000000}));
        bucket4.addToList(SuperKmer({0b01011011}));
//        cout << "------ 4 -------" << endl;
//        bucket4.print();
        EXPECT(bucket4.isSorted());

//        cout << "----- 1 & 4 -----" << endl;
//        (bucket1 & bucket4).print();
        EXPECT((bucket1 & bucket4).getListSize() == (uint64_t) 3);
        EXPECT((bucket1 & bucket4).isSorted());
        EXPECT((bucket1 & bucket4).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 & bucket4).find(Kmer(0b0000001000), position));
        EXPECT((bucket1 & bucket4).find(Kmer(0b0000001010), position));


        bucket5.addToList(SuperKmer({0b01011000}));
        bucket5.addToList(SuperKmer({0b00100100, 0b10000000}));
        bucket5.addToList(SuperKmer({0b00101000, 0b11000000}));
        bucket5.addToList(SuperKmer({0b10101100, 0b01100000}));
        EXPECT(bucket5.isSorted());
        EXPECT((bucket1 & bucket5).getListSize() == (uint64_t) 0);

        Bucket empty = Bucket(0);
        EXPECT((bucket1 & empty).getListSize() == (uint64_t) 0);


    },
};


int main() {
    int failures(0);
    sk7::initLib(5, 3);
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
    cout << "*** Test on utils function(s) ***" << endl;
    if((failures = run(utilsTest))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    cout << "*** Test on SuperKmer method(s) ***" << endl;
    if((failures = run(superKmerTest))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    cout << "*** Test on bucket methods ***" << endl;
    if((failures = run(bucketTest))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
    return 0;
}
