
#include "Bucket.hpp"
#include "lest.hpp"

using namespace std;
using namespace lest;

const test superKmerTest[] {
        CASE("split") {
            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});

            std::vector<SuperKmer> split1 = SK1.split();
            EXPECT(split1.size() == (uint64_t) 11);
            bool equal = split1.at(0) == SuperKmer({0b10100000, 0b00010010, 0b00010011, 0b00110000, 0b00110010, 0b00010000});
            EXPECT(equal);
            equal = split1.at(1) == SuperKmer({0b10010001, 0b11010010, 0b00010011, 0b00110000, 0b00110010, 0b00010000});
            EXPECT(equal);
            equal = split1.at(2) == SuperKmer({0b10000010, 0b11011010, 0b00010011, 0b00110000, 0b00110010});
            EXPECT(equal);
            equal = split1.at(3) == SuperKmer({0b01110011, 0b11011010, 0b01010011, 0b00110000, 0b00110000});
            EXPECT(equal);
            equal = split1.at(4) == SuperKmer({0b01100100, 0b11011010, 0b01010011, 0b00110000});
            EXPECT(equal);
            equal = split1.at(5) == SuperKmer({0b01010101, 0b11011010, 0b01010011, 0b00110000});
            EXPECT(equal);
            equal = split1.at(6) == SuperKmer({0b01000110, 0b11011010, 0b01010011, 0b00001000});
            EXPECT(equal);
            equal = split1.at(7) == SuperKmer({0b00110111, 0b11011010, 0b01010000, 0b00001000, 0b01000000});
            EXPECT(equal);
            equal = split1.at(8) == SuperKmer({0b00101000, 0b11011010, 0b01000000, 0b00001000, 0b01001100});
            EXPECT(equal);
            equal = split1.at(9) == SuperKmer({0b00011001, 0b11011000, 0b01000000, 0b00001000, 0b01001100, 0b11000000});
            EXPECT(equal);
            equal = split1.at(10) == SuperKmer({0b00001010, 0b11001000, 0b01000000, 0b00001000, 0b01001100, 0b11000000});
            EXPECT(equal);

            std::vector<SuperKmer> split2 = SK2.split();
            EXPECT(split2.size() == (uint64_t) 11);
            equal = split2.at(0) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(1) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(2) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(3) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(4) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(5) == SuperKmer({0b01010101, 0b01101001, 0b00000111, 0b01100000});
            EXPECT(equal);
            equal = split2.at(6) == SuperKmer({0b01000110, 0b01101001, 0b00000111, 0b01001000});
            EXPECT(equal);
            equal = split2.at(7) == SuperKmer({0b00110111, 0b01101001, 0b00000100, 0b01001000, 0b11000000});
            EXPECT(equal);
            equal = split2.at(8) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(9) == SuperKmer();
            EXPECT(equal);
            equal = split2.at(10) == SuperKmer();
            EXPECT(equal);

            },
        CASE("CompareSK") {
            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});
            SuperKmer SK3 = SuperKmer({0b10101010, 0b11011010, 0b01110011, 0b00111000, 0b01111110, 0b11010000});

            std::vector<SuperKmer::logic> expected = {SuperKmer::EQUAL, SuperKmer::EQUAL, SuperKmer::EQUAL,
                                                      SuperKmer::EQUAL, SuperKmer::EQUAL, SuperKmer::EQUAL,
                                                      SuperKmer::EQUAL, SuperKmer::EQUAL, SuperKmer::EQUAL,
                                                      SuperKmer::EQUAL, SuperKmer::EQUAL};
            EXPECT(SuperKmer::compareSK(SK1, SK1) == expected);

            expected = {SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE,
                        SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE, SuperKmer::SUPERIOR,
                        SuperKmer::SUPERIOR, SuperKmer::SUPERIOR, SuperKmer::INCOMPARABLE,
                        SuperKmer::INCOMPARABLE, SuperKmer::INCOMPARABLE};
            EXPECT(SuperKmer::compareSK(SK1, SK2) == expected);

            expected = {SuperKmer::INFERIOR, SuperKmer::INFERIOR, SuperKmer::INFERIOR,
                        SuperKmer::INFERIOR, SuperKmer::INFERIOR, SuperKmer::INFERIOR,
                        SuperKmer::INFERIOR, SuperKmer::INFERIOR, SuperKmer::EQUAL,
                        SuperKmer::EQUAL, SuperKmer::EQUAL};
            EXPECT(SuperKmer::compareSK(SK1, SK3) == expected);


        },

        CASE("readKmer") {
            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});

            EXPECT(SK1.readKmer(0).getValue() == (uint64_t) 0b1100100001000000000010000100110011000000);
            EXPECT(SK1.readKmer(1).getValue() == (uint64_t) 0b110110000100000000001000010011001100);
            EXPECT(SK1.readKmer(2).getValue() == (uint64_t) 0b11011010010000000000100001001100);
            EXPECT(SK1.readKmer(3).getValue() == (uint64_t) 0b1101101001010000000010000100);
            EXPECT(SK1.readKmer(4).getValue() == (uint64_t) 0b110110100101001100001000);
            EXPECT(SK1.readKmer(5).getValue() == (uint64_t) 0b11011010010100110011);
            EXPECT(SK1.readKmer(6).getValue() == (uint64_t) 0b110110100101001100110000);
            EXPECT(SK1.readKmer(7).getValue() == (uint64_t) 0b1101101001010011001100000011);
            EXPECT(SK1.readKmer(8).getValue() == (uint64_t) 0b11011010000100110011000000110010);
            EXPECT(SK1.readKmer(9).getValue() == (uint64_t) 0b110100100001001100110000001100100001);
            EXPECT(SK1.readKmer(10).getValue() == (uint64_t) 0b01001000010011001100000011001000010000);

            EXPECT(SK2.readKmer(0).length == 0);
            EXPECT(SK2.readKmer(1).length == 0);
            EXPECT(SK2.readKmer(2).length == 0);
            EXPECT(SK2.readKmer(3).getValue() == (uint64_t) 0b0110100100000100010010001100);
            EXPECT(SK2.readKmer(4).getValue() == (uint64_t) 0b011010010000011101001000);
            EXPECT(SK2.readKmer(5).getValue() == (uint64_t) 0b01101001000001110110);
            EXPECT(SK2.readKmer(6).length == 0);
            EXPECT(SK2.readKmer(7).length == 0);
            EXPECT(SK2.readKmer(8).length == 0);
            EXPECT(SK2.readKmer(9).length == 0);
            EXPECT(SK2.readKmer(10).length == 0);

        }
};

const test bucketTest[] {
        CASE("find") {
            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});
            SuperKmer SK3 = SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000});

            Bucket bucket(0);
            bucket.addToList(SK2);
            bucket.addToList(SK1);
            bucket.addToList(SK3);
            EXPECT(bucket.isSorted());

            Kmer toSearch1 = Kmer(0b10010010100000000000000110000101); // TCACT|CTACC
            Kmer toSearch2 = Kmer(0b110001110000000000011000010110); // GACG|CTACCT
            Kmer toSearch3 = Kmer(0b001111111011000000000011110100); // AGGGTG|GGCA
            Kmer toSearch4 = Kmer(0b011011001111011001000000000011); // CTGAGGCTC|G
            Kmer toSearch5 = Kmer(0b110001100000000000011000010110); // GACT|CTACCT
            Kmer toSearch6 = Kmer(0b001111111010000000000011110100); // AGGGTT|GGCA
            Kmer toSearch7 = Kmer(0b111101100100000000001110010000); // GGCTC|GTCAA
            Kmer toSearch8 = Kmer(0b000000000011100100001001111100); // |GTCAATCGGA
            Kmer toSearch9 = Kmer(0b000110000000000001100001011011); // ACT|CTACCTG

            int position;

            EXPECT(not bucket.find(toSearch1, position)); EXPECT(position == 0);
            EXPECT(not bucket.find(toSearch2, position)); EXPECT(position == 1);
            EXPECT(not bucket.find(toSearch3, position)); EXPECT(position == 3);
            EXPECT(bucket.find(toSearch4, position)); EXPECT(position == 1);
            EXPECT(bucket.find(toSearch5, position)); EXPECT(position == 0);
            EXPECT(bucket.find(toSearch6, position)); EXPECT(position == 2);
            EXPECT(bucket.find(toSearch7, position)); EXPECT(position == 1);
            EXPECT(bucket.find(toSearch8, position)); EXPECT(position == 1);
            EXPECT(bucket.find(toSearch9, position)); EXPECT(position == 0);

        },
        CASE("isSorted") {

            Bucket bucket(0);
            bucket.addToList(SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000}));
            bucket.addToList(SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000}));
            bucket.addToList(SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000}));

            Bucket spaced = Bucket(0);
            spaced.addToList(SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000}));
            spaced.addToList(SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000}));
            spaced.addToList(SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000}));
            spaced.addToList(SuperKmer({0b10100000, 0b00010010, 0b00010011, 0b00110000, 0b00110010, 0b00010000}));

            EXPECT(bucket.isSorted());
            EXPECT(not spaced.isSorted());
        }
        ,
        CASE("addKmer") {
            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});
            SuperKmer SK3 = SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000});

            Bucket bucket(0);
            bucket.addToList(SK2);
            bucket.addToList(SK1);
            bucket.addToList(SK3);
            EXPECT(bucket.isSorted());

            Kmer toAdd1 = Kmer(0b10010010100000000000000110000101); // TCACT|CTACC
            Kmer toAdd2 = Kmer(0b110001110000000000011000010110);   // GACG|CTACCT
            Kmer toAdd3 = Kmer(0b001111111011000000000011110100);   // AGGGTG|GGCA
            Kmer toAdd4 = Kmer(0b011011001111011001000000000011);   // CTGAGGCTC|G
            Kmer toAdd5 = Kmer(0b110001100000000000011000010110);   // GACT|CTACCT
            Kmer toAdd6 = Kmer(0b001111111010000000000011110100);   // AGGGTT|GGCA
            Kmer toAdd7 = Kmer(0b111101100100000000001110010000);   // GGCTC|GTCAA
            Kmer toAdd8 = Kmer(0b000000000011100100001001111100);   // |GTCAATCGGA
            Kmer toAdd9 = Kmer(0b000110000000000001100001011011);   // ACT|CTACCTG

            EXPECT(bucket.getListSize() == (uint64_t) 3);
            bucket.addKmer(toAdd1);
            EXPECT(bucket.getListSize() == (uint64_t) 4);
            bucket.addKmer(toAdd2);
            EXPECT(bucket.getListSize() == (uint64_t) 5);
            bucket.addKmer(toAdd3);
            EXPECT(bucket.getListSize() == (uint64_t) 6);
            bucket.addKmer(toAdd1);
            EXPECT(bucket.getListSize() == (uint64_t) 6);
            EXPECT(bucket.isSorted());

            int position;
            EXPECT(bucket.find(toAdd1, position) ); EXPECT(position == 0); // TCACT|CTACC
            EXPECT(bucket.find(toAdd2, position)); EXPECT(position == 2);  // GACG|CTACCT
            EXPECT(bucket.find(toAdd3, position)); EXPECT(position == 5);  // AGGGTG|GGCA
            EXPECT(bucket.find(toAdd4, position)); EXPECT(position == 3);  // CTGAGGCTC|G
            EXPECT(bucket.find(toAdd5, position)); EXPECT(position == 1);  // GACT|CTACCT
            EXPECT(bucket.find(toAdd6, position)); EXPECT(position == 4);  // AGGGTT|GGCA
            EXPECT(bucket.find(toAdd7, position)); EXPECT(position == 3);  // GGCTC|GTCAA
            EXPECT(bucket.find(toAdd8, position)); EXPECT(position == 3);  // |GTCAATCGGA
            EXPECT(bucket.find(toAdd9, position)); EXPECT(position == 1);  // ACT|CTACCTG
        },
        CASE("SKtoKmer") {

            Bucket bucket = Bucket(0);
            EXPECT(bucket.SKtoKmer(SuperKmer({0b01010101, 0b01001010, 0b00100100, 0b01010000})).getValue() == (uint64_t) 0b010010100000000000000110000101);
            EXPECT(bucket.SKtoKmer(SuperKmer({0b01000110, 0b01111001, 0b00000111, 0b01001000})).getValue() == (uint64_t) 0b110001110000000000011000010110);

        },

        CASE("addSuperKmer") {

            SuperKmer SK1 = SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000});
            SuperKmer SK2 = SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000});
            SuperKmer SK3 = SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000});

            Bucket bucket(0);

            bucket.addSuperKmer(SK3);
            EXPECT(bucket.getListSize() == (uint64_t) 1);
            bucket.addSuperKmer(SK2);
            EXPECT(bucket.getListSize() == (uint64_t) 4);
            bucket.addSuperKmer(SK2);
            EXPECT(bucket.getListSize() == (uint64_t) 4);
            bucket.addSuperKmer(SK1);
            EXPECT(bucket.getListSize() == (uint64_t) 15);
            bucket.addSuperKmer(SK1);
            EXPECT(bucket.getListSize() == (uint64_t) 15);
            bucket.addSuperKmer(SK2);
            EXPECT(bucket.getListSize() == (uint64_t) 15);
            bucket.addSuperKmer(SK3);
            EXPECT(bucket.getListSize() == (uint64_t) 15);

            EXPECT(bucket.isSorted());
        },

        CASE("Intersection") {
            Bucket bucket1 = Bucket(0);
            Bucket bucket2 = Bucket(0);
//            Bucket bucket3 = Bucket(0);
//            Bucket bucket4 = Bucket(0);
//            Bucket bucket5 = Bucket(0);

//            int position;

            bucket1.addToList(SuperKmer({0b01010111, 0b01101001, 0b00000111, 0b01101000, 0b11000000}));
            bucket1.addToList(SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000}));
            bucket1.addToList(SuperKmer({0b01100110, 0b11101110, 0b01110011, 0b00111000}));
            cout << "------ 1 ------" << endl;
            bucket1.print();
            EXPECT(bucket1.isSorted());

            bucket2.addToList(SuperKmer({0b01010111, 0b01001001, 0b00000111, 0b01101000, 0b11000000}));
            bucket2.addToList(SuperKmer({0b10101010, 0b11011010, 0b01010011, 0b00111000, 0b01111110, 0b11010000}));
            bucket2.addToList(SuperKmer({0b01100110, 0b11111110, 0b01110011, 0b00111000}));
            cout << "----- 2 ------" << endl;
            bucket2.print();
            EXPECT(bucket2.isSorted());

            cout << "---- 1 & 2 ----" << endl;
            (bucket1 & bucket2).print();
            EXPECT((bucket1 & bucket2).getListSize() == (uint64_t) 11);
            EXPECT((bucket1 & bucket2).isSorted());

            Bucket empty = Bucket(0);
            EXPECT((bucket1 & empty).getListSize() == (uint64_t) 0);


        },
};


int main() {
    int failures(0);
    sk7::initLib(15, 5);
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