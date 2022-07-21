#include "lest.hpp"
#include "Bucket.hpp"

#include <iostream>

using namespace lest;
using namespace std;

/**
 * reverse inf
 * @param kmer1 the first Kmer
 * @param kmer2 the second Kmer
 * @return kmer2 < kmer1
 */
bool newInf(const Kmer &kmer1, const Kmer &kmer2) {
    return kmer2.getValue() < kmer1.getValue();
}

const test bucketTest[] {
    CASE("find") {
        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));
        EXPECT(bucket.isSorted());


        // The bucket contains : GG|, CG|TT, G|CT

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

        EXPECT(not bucket.find(toSearch, position)); EXPECT(position == 2);
        EXPECT(not bucket.find(toSearch2, position)); EXPECT(position == 3);
        EXPECT(not bucket.find(toSearch3, position)); EXPECT(position == 2);
        EXPECT(bucket.find(toSearch4, position)); EXPECT(position == 2);
        EXPECT(bucket.find(toSearch5, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch6, position)); EXPECT(position == 0);
        EXPECT(bucket.find(toSearch7, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch8, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch9, position)); EXPECT(position == 2);

    },

    CASE("isSorted") {

        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));


        Bucket spaced = Bucket(0);
        spaced.addToList(SuperKmer({0b10000001, 0b00010000}));
        spaced.addToList(SuperKmer({0b01011001}));
        spaced.addToList(SuperKmer({0b10000010, 0b00100000}));

        EXPECT(bucket.isSorted());
        EXPECT(not spaced.isSorted());
    },

    CASE("addKmer") {
        Bucket bucket(0);
        bucket.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket.addToList(SuperKmer({0b01100111, 0b10000000}));

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
        EXPECT(bucket.find(toAdd, position)); EXPECT(position == 3);
        EXPECT(bucket.find(toAdd2, position)); EXPECT(position == 5);
        EXPECT(bucket.find(toAdd3, position)); EXPECT(position == 2);
        EXPECT(bucket.find(toSearch4, position)); EXPECT(position == 4);
        EXPECT(bucket.find(toSearch5, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch6, position)); EXPECT(position == 0);
        EXPECT(bucket.find(toSearch7, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch8, position)); EXPECT(position == 1);
        EXPECT(bucket.find(toSearch9, position)); EXPECT(position == 4);
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


    CASE("Intersection") {
        Bucket bucket1 = Bucket(0);
        Bucket bucket2 = Bucket(0);
        Bucket bucket3 = Bucket(0);
        Bucket bucket4 = Bucket(0);
        Bucket bucket5 = Bucket(0);

        int position;

        bucket1.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket1.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket1.addToList(SuperKmer({0b00101000, 0b00000000}));
//        cout << "------ 1 ------" << endl;
//        bucket1.print();
        EXPECT(bucket1.isSorted());

        bucket2.addToList(SuperKmer({0b00101100, 0b11000000}));
        bucket2.addToList(SuperKmer({0b01011101}));
        bucket2.addToList(SuperKmer({0b10100111, 0b10010000}));
//        cout << "----- 2 ------" << endl;
//        bucket2.print();
        EXPECT(bucket2.isSorted());

//        cout << "---- 1 & 2 ----" << endl;
//        (bucket1 & bucket2).print();
        EXPECT((bucket1 & bucket2).getListSize() == (uint64_t) 1);
        EXPECT((bucket1 & bucket2).find(Kmer(0b0111000000), position));


        bucket3.addToList(SuperKmer({0b10101010, 0b11100000}));
        bucket3.addToList(SuperKmer({0b10101001, 0b10010000}));
        bucket3.addToList(SuperKmer({0b01010110}));
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

        bucket4.addToList(SuperKmer({0b01011011}));
        bucket4.addToList(SuperKmer({0b01101001, 0b10000000}));
        bucket4.addToList(SuperKmer({0b00101000, 0b00000000}));
        bucket4.addToList(SuperKmer({0b10000001, 0b00000000}));
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


        bucket5.addToList(SuperKmer({0b10101100, 0b01100000}));
        bucket5.addToList(SuperKmer({0b00101000, 0b11000000}));
        bucket5.addToList(SuperKmer({0b00100100, 0b10000000}));
        bucket5.addToList(SuperKmer({0b01011000}));
        EXPECT(bucket5.isSorted());
        EXPECT((bucket1 & bucket5).getListSize() == (uint64_t) 0);

        Bucket empty = Bucket(0);
        EXPECT((bucket1 & empty).getListSize() == (uint64_t) 0);

    },

    CASE("Union") {
        Bucket bucket1 = Bucket(0);
        Bucket bucket2 = Bucket(0);
        Bucket bucket3 = Bucket(0);
        Bucket bucket4 = Bucket(0);
        Bucket bucket5 = Bucket(0);

        int position;

        bucket1.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket1.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket1.addToList(SuperKmer({0b00101000, 0b00000000}));
//        cout << "------ 1 ------" << endl;
//        bucket1.print();
        EXPECT(bucket1.isSorted());

        bucket2.addToList(SuperKmer({0b00101100, 0b11000000}));
        bucket2.addToList(SuperKmer({0b01011101}));
        bucket2.addToList(SuperKmer({0b10100111, 0b10010000}));
//        cout << "----- 2 ------" << endl;
//        bucket2.print();
        EXPECT(bucket2.isSorted());

//        cout << "---- 1 | 2 ----" << endl;
        Bucket test12 = bucket1 | bucket2;
//        test12.print();
        EXPECT(test12.getListSize() == (uint64_t) 9);
        EXPECT(test12.find(Kmer(0b0111000000), position));
        EXPECT(test12.isSorted());

        bucket3.addToList(SuperKmer({0b10101010, 0b11100000}));
        bucket3.addToList(SuperKmer({0b10101001, 0b10010000}));
        bucket3.addToList(SuperKmer({0b01010110}));
        bucket3.addSuperKmer(SuperKmer({0b01011011}));
        bucket3.addSuperKmer(SuperKmer({0b00101100, 0b01000000}));
//        cout << "------ 3 ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 | 3 --------" << endl;
//        (bucket1 | bucket3).print();
        EXPECT((bucket1 | bucket3).getListSize() == (uint64_t) 12);
        EXPECT((bucket1 | bucket3).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 | bucket3).find(Kmer(0b0000001010), position));


        bucket3.addSuperKmer(SuperKmer({0b10000011, 0b00110000}));
//        cout << "------ 3++ ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 | 3++ --------" << endl;
//        (bucket1 | bucket3).print();
        EXPECT((bucket1 | bucket3).getListSize() == (uint64_t) 12);
        EXPECT((bucket1 | bucket3).isSorted());
        EXPECT((bucket1 | bucket3).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 | bucket3).find(Kmer(0b0000001010), position));
        EXPECT((bucket1 | bucket3).find(Kmer(0b1111000000), position));

        bucket4.addToList(SuperKmer({0b01011011}));
        bucket4.addToList(SuperKmer({0b01101001, 0b10000000}));
        bucket4.addToList(SuperKmer({0b00101000, 0b00000000}));
        bucket4.addToList(SuperKmer({0b10000001, 0b00000000}));
//        cout << "------ 4 -------" << endl;
//        bucket4.print();
        EXPECT(bucket4.isSorted());

//        cout << "----- 1 | 4 -----" << endl;
//        (bucket1 | bucket4).print();
        EXPECT((bucket1 | bucket4).getListSize() == (uint64_t) 7);
        EXPECT((bucket1 | bucket4).isSorted());
        EXPECT((bucket1 | bucket4).find(Kmer(0b1100000010), position));
        EXPECT((bucket1 | bucket4).find(Kmer(0b0000001000), position));
        EXPECT((bucket1 | bucket4).find(Kmer(0b0000001010), position));


        bucket5.addToList(SuperKmer({0b10101100, 0b01100000}));
        bucket5.addToList(SuperKmer({0b00101000, 0b11000000}));
        bucket5.addToList(SuperKmer({0b00100100, 0b10000000}));
        bucket5.addToList(SuperKmer({0b01011000}));
//        cout << "------ 5 ------" << endl;
//        bucket5.print();
        EXPECT(bucket5.isSorted());

//        cout << "------ 1 | 5 ------" << endl;
        Bucket bucket15 = bucket1 | bucket5;
//        bucket15.print();
        EXPECT(bucket15.getListSize() == (uint64_t) 11);

        Bucket empty = Bucket(0);
//        (bucket1 | empty).print();
        EXPECT((bucket1 | empty).getListSize() == (uint64_t) 5);

    },

    CASE("Symmetrical difference") {
        Bucket bucket1 = Bucket(0);
        Bucket bucket2 = Bucket(0);
        Bucket bucket3 = Bucket(0);
        Bucket bucket4 = Bucket(0);
        Bucket bucket5 = Bucket(0);

        int position;

        bucket1.addToList(SuperKmer({0b10000011, 0b00110000}));
        bucket1.addToList(SuperKmer({0b10101011, 0b10010000}));
        bucket1.addToList(SuperKmer({0b00101000, 0b00000000}));
//        cout << "------ 1 ------" << endl;
//        bucket1.print();
        EXPECT(bucket1.isSorted());

        bucket2.addToList(SuperKmer({0b00101100, 0b11000000}));
        bucket2.addToList(SuperKmer({0b01011101}));
        bucket2.addToList(SuperKmer({0b10100111, 0b10010000}));
//        cout << "----- 2 ------" << endl;
//        bucket2.print();
        EXPECT(bucket2.isSorted());

//        cout << "---- 1 ^ 2 ----" << endl;
        Bucket test12 = bucket1 ^ bucket2;
//        test12.print();
        EXPECT(test12.getListSize() == (uint64_t) 8);
        EXPECT(not test12.find(Kmer(0b0111000000), position));
        EXPECT(test12.isSorted());

        bucket3.addToList(SuperKmer({0b10101010, 0b11100000}));
        bucket3.addToList(SuperKmer({0b10101001, 0b10010000}));
        bucket3.addToList(SuperKmer({0b01010110}));
        bucket3.addSuperKmer(SuperKmer({0b01011011}));
        bucket3.addSuperKmer(SuperKmer({0b00101100, 0b01000000}));
//        cout << "------ 3 ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 ^ 3 --------" << endl;
//        (bucket1  ^ bucket3).print();
        EXPECT((bucket1 ^ bucket3).getListSize() == (uint64_t) 10);
        EXPECT(not (bucket1 ^ bucket3).find(Kmer(0b1100000010), position));
        EXPECT(not (bucket1 ^ bucket3).find(Kmer(0b0000001010), position));


        bucket3.addSuperKmer(SuperKmer({0b10000011, 0b00110000}));
//        cout << "------ 3++ ------" << endl;
//        bucket3.print();
        EXPECT(bucket3.isSorted());

//        cout << "------ 1 ^ 3++ --------" << endl;
//        (bucket1  ^ bucket3).print();
        EXPECT((bucket1 ^ bucket3).getListSize() == (uint64_t) 9);
        EXPECT((bucket1 ^ bucket3).isSorted());
        EXPECT(not (bucket1 ^ bucket3).find(Kmer(0b1100000010), position));
        EXPECT(not (bucket1 ^ bucket3).find(Kmer(0b0000001010), position));
        EXPECT(not (bucket1 ^ bucket3).find(Kmer(0b1111000000), position));

        bucket4.addToList(SuperKmer({0b01011011}));
        bucket4.addToList(SuperKmer({0b01101001, 0b10000000}));
        bucket4.addToList(SuperKmer({0b00101000, 0b00000000}));
        bucket4.addToList(SuperKmer({0b10000001, 0b00000000}));
//        cout << "------ 4 -------" << endl;
//        bucket4.print();
        EXPECT(bucket4.isSorted());

//        cout << "----- 1 ^ 4 -----" << endl;
//        (bucket1 | bucket4).print();
        EXPECT((bucket1 ^ bucket4).getListSize() == (uint64_t) 4);
        EXPECT((bucket1 ^ bucket4).isSorted());
        EXPECT(not (bucket1 ^ bucket4).find(Kmer(0b1100000010), position));
        EXPECT(not (bucket1 ^ bucket4).find(Kmer(0b0000001000), position));
        EXPECT(not (bucket1 ^ bucket4).find(Kmer(0b0000001010), position));


        bucket5.addToList(SuperKmer({0b10101100, 0b01100000}));
        bucket5.addToList(SuperKmer({0b00101000, 0b11000000}));
        bucket5.addToList(SuperKmer({0b00100100, 0b10000000}));
        bucket5.addToList(SuperKmer({0b01011000}));
//        cout << "------ 5 ------" << endl;
//        bucket5.print();
        EXPECT(bucket5.isSorted());

//        cout << "------ 1 ^ 5 ------" << endl;
        Bucket bucket15 = bucket1 ^ bucket5;
//        bucket15.print();
        EXPECT(bucket15.getListSize() == (uint64_t) 11);

        Bucket empty = Bucket(0);
//        (bucket1 ^ empty).print();
        EXPECT((bucket1 ^ empty).getListSize() == (uint64_t) 5);

    }
};

int main() {
    int failures(0);
    sk7::initLib(5, 3, &newInf);
    cout << "*** Test on bucket methods ***" << endl;
    if((failures = run(bucketTest))) {
        return failures;
    }
    cout << "*** Passed ***" << endl;
}