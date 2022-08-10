#include <iostream>

#include "Kff_scanner.hpp"

using namespace std;

// TTCATTTTTTATCACATA
// TGATCAATATATTAAACA

int main(int argc, char** argv) {

    cout << "I'm main" << endl;
    if (argc < 2) {
        cout << "usage : ./main <file_name>" << endl;
        return 1;
    }

    Kff_scanner scanner = Kff_scanner(argv[1], false, false, false, 10);

    cout << "here" << endl;
    BucketMap* read = scanner.readAll();

    for (auto &it : *read->map) {
        cout << "New Bucket minimizer = " << bitset<20>(it.first) << " of size : " << it.second.getListSize() <<  endl;
        it.second.print();
        if (not it.second.isSorted()) {
            cout << "fucked up" << endl;
            exit(1);
        }
        cout << "#######" << endl;
    }

    cout << "end" << endl;

    delete read;

    return 0;
}
