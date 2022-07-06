#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(5, 3);

    SuperKmer SK = SuperKmer({0b00101000, 0b00000000});
    SK.print();
    cout << "res = " << bitset<16>(SK.readKmer(0).getValue()) << endl;
    return 0;
}
