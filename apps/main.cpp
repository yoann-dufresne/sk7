#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(5, 3);

    Bucket bucket1 = Bucket(0);
    Bucket bucket2 = Bucket(0);

    cout << "### 1 ###" << endl;
    bucket1.addToList(SuperKmer({0b00101000, 0b00000000}));
    bucket1.addToList(SuperKmer({0b10101011, 0b10010000}));
    bucket1.addToList(SuperKmer({0b10000011, 0b00110000}));
    bucket1.print();

//    cout << "no" << endl;
//    bucket1.getListCopy().at(0).extract(0).print();
//    bucket1.getListCopy().at(0).extract(1).print();
//    bucket1.getListCopy().at(0).extract(2).print();
//    bucket1.getListCopy().at(1).extract(0).print();
//    bucket1.getListCopy().at(1).extract(1).print();
//    bucket1.getListCopy().at(1).extract(2).print();
//    bucket1.getListCopy().at(2).extract(0).print();
//    bucket1.getListCopy().at(2).extract(1).print();
//    bucket1.getListCopy().at(2).extract(2).print();

    cout << "### 2 ###" << endl;
    bucket2.addToList(SuperKmer({0b10100111, 0b10010000}));
    bucket2.addToList(SuperKmer({0b01011101}));
    bucket2.addToList(SuperKmer({0b00101100, 0b11000000}));
    bucket2.print();

    cout << "#########" << endl;
    Bucket::chainedUnion(bucket1, bucket2).print();

    cout << "### 3 ###" << endl;
    Bucket bucket3 = Bucket(0);
    bucket3.addToList(SuperKmer({0b00101000, 0b00000000}));
    bucket3.addToList(SuperKmer({0b10101011, 0b10010000}));
    bucket3.addToList(SuperKmer({0b10000011, 0b00110000}));
    bucket3.print();
    if (not bucket3.isSorted()) cout << "bhncdazhdrfzi,rfot" << endl << endl << endl;

    cout << "### 4 ###" << endl;
    Bucket bucket4 = Bucket(0);
//    bucket4.addToList(SuperKmer({0b10000001, 0b00000000}));
    bucket4.addToList(SuperKmer({0b01010110}));
    bucket4.addToList(SuperKmer({0b10101001, 0b10010000}));
    bucket4.addToList(SuperKmer({0b10101010, 0b11100000}));
    bucket4.addToList(SuperKmer({0b00101100, 0b01000000}));
    bucket4.addToList(SuperKmer({0b01011011}));
    bucket4.print();
    if (not bucket4.isSorted()) cout << "bhncdazhdrfzi,rfot" << endl << endl << endl;

    cout << "#########" << endl;
    Bucket::chainedUnion(bucket3, bucket4).print();

    cout << "### 5 ###" << endl;
    Bucket bucket5 = Bucket(0);
    bucket5.addToList(SuperKmer({0b00100100, 0b00000000}));
    bucket5.addToList(SuperKmer({0b10100101, 0b10010000}));
    bucket5.addToList(SuperKmer({0b10000001, 0b00110000}));
    bucket5.print();
    if (not bucket5.isSorted()) cout << "bhncdazhdrfzi,rfot" << endl << endl << endl;

    cout << "#########" << endl;
    Bucket::chainedUnion(bucket5, bucket4).print();
//    cout << Bucket::chainedUnion(bucket5, bucket4).isSorted() << endl;

    return 0;
}
