#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(5, 3);

    SuperKmer SK = SuperKmer({0b01100111, 0b10000000});
    SK.print();
    cout << "res = " << bitset<16>(SK.readKmer(1).getValue()) << endl;
    return 0;
}
