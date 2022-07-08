#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(15, 5);

    SuperKmer SK1 = SuperKmer({0b01100100, 0b11111110, 0b01110011, 0b00110000});
    SK1.print();
    Kmer test = SK1.readKmer(6);
    cout << test.toString() << endl;

//    for (auto it : SK1.split()) {
//        it.print();
//    }
    return 0;
}
