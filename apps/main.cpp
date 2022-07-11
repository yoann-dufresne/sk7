#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(8, 3);

    Kmer test = Kmer(0b1100000010101001); // GAAATTTC
    cout << " init : " << test.toString() << endl;
    Minimiser a = Minimiser(alpha, test);
    cout << "after : " << test.toString() << endl;
    cout << "Minimiser : " << a.toString() << " pos = " << a.getPos() << endl;
    return 0;
}
