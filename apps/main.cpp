#include <iostream>

#include "Bucket.hpp"

using namespace std;

int main() {

    cout << "I'm main" << endl;
    sk7::initLib(5, 3);


    cout << reorderValue(12,  2,0) << endl;

    return 0;
}
