#include <iostream>

#include "Kff_scanner.hpp"

using namespace std;

int main(int argc, char** argv) {

    cout << "I'm main" << endl;
    if (argc < 2) {
        cout << "usage : ./main <file_name>" << endl;
        return 1;
    }

    Kff_scanner scanner = Kff_scanner(argv[1], false, false, false, 12);

    cout << "here" << endl;
    std::unordered_map<uint64_t, sk7::Bucket> read = scanner.readAll();

    uint64_t cmpt = 0;

    for (auto &it : read) {

        if (not it.second.isSorted()) {
            cout << bitset<24>(it.first) << endl;
            it.second.print();
            cout << "cmpt = " << cmpt << endl;
            exit(1);
        }
        cmpt++;
    }

    cout << "end" << endl;

    return 0;
}
