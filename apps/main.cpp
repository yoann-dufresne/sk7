#include <iostream>
#include <fcntl.h>

#include "Kff_scanner.hpp"

#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */

using namespace std;

void help() {
    std::cout << "Usage : ./apps/main [OPTIONS]\n"
                 "Options :\n"
                 "\t-h : show this message and exit\n"
                 "\t-i <input_path> : indicate the input path of the Kff file" << RED << " REQUIRED" << RESET << "\n"
                 "\t-k <kmer_file_path> : indicate the path to a file containing Kmers to test, 1 Kmer per line"
                 << RED << " REQUIRED" << RESET << "\n"
                 "\t-b : put if the input is already a bucketed kff file\n"
                 "\t-m <m> : specifies the wanted minimizer length," << RED << " REQUIRED" << RESET << " if -b not present\n"
                 "\t-s : the input is sorted\n" << std::endl;
    exit(0);
}

void search_random_kmer(const uint64_t &n, BucketMap* &map) {
    srand(time(nullptr));
    uint64_t upper_limit = (2 << (sk7::k * 2)) - 1;
    int position = 0;
    for (uint64_t i = 0; i < n; i++) {
        Kmer to_search = Kmer(rand() % upper_limit);
        map->find(to_search, position);
    }
}

Kmer read_kmer(const string &line) {
    string kmer = line.substr(0, line.find('\t'));
    uint64_t value = 0;
    for (char c : kmer) {
        value <<= 2;
        value += (c >> 1) & 0b11;
    }
    return Kmer(value);
}

void search_kmer_from_file(BucketMap* &map, char* kmer_file) {
    ifstream file;
    file.open(kmer_file);
    string line;
    int position;
    uint64_t i = 0;
    while(getline(file, line)) {
        cout << "###\nSearching : " << read_kmer(line).toString() << endl;
        if (not map->find(read_kmer(line), position)) {
            cout << "i = " << i << " not found = " << read_kmer(line).toString() << endl;
            cout << "position = " << position << endl;
            exit(1);
        }
        i++;
    }
    file.close();
}

int main(int argc, char** argv) {

    /// Set values
    bool bucketed = false;
    bool sorted = false;
    char* input_name;
    char* kmer_name;
    bool set_input = false;
    int m = -1;

    /// Parse options
    int c;
    while((c = getopt(argc, argv, "bshm:i:k:")) != -1) {
        switch (c) {
            case 'h':
                sk7::initLib(4, 1);
                cout << Kmer("ACGT").toString() << endl;
                help();
                break;
            case 'k':
                kmer_name = optarg;
                break;
            case 'i':
                input_name = optarg;
                set_input = true;
                break;
            case 'b':
                bucketed = true;
                break;
            case 's':
                bucketed = true;
                sorted = true;
                break;
            case 'm':
                m = atoi(optarg);
                break;
            case '?':
                if (optopt == 'm') {
                    cerr << "m option requires an argument" << endl;
                    exit(1);
                }
                else if (isprint(optopt)) {
                    cerr << "Unknown option" << endl;
                    exit(1);
                }
                else {
                    cerr << "Unknown character" << endl;
                }
                break;
            default:
                abort();
        }
    }

    /// Verify setup
    if (not set_input) {
        cerr << "Input file needs to be specified via -i <input_path>" << endl;
        exit(1);
    }

    if (not bucketed && m <= 0) {
        cerr << "An unbucketed input needs ab specified m via -m <m_value> where 1 <= m_value <= 31" << endl;
        exit(1);
    }

    if (kmer_name == nullptr) {
        cerr << "You need to specify a file with the k-mers to search" << endl;
        exit(1);
    }

    /// Launch the read
    Kff_scanner* scanner = new Kff_scanner(input_name, m, bucketed, sorted);
    BucketMap* read = scanner->readAll();
    delete scanner;

    for (auto &it : *read->map) {
        if (not it.second.isSorted()) {
            it.second.print();
            cout << "fucked up" << endl;
            exit(1);
        }
    }

    /// Random searches
    cout << "--- random searches ---" << endl;
    search_random_kmer(10, read);

    /// good searches
    cout << "--- deterministic searches ---" << endl;
    search_kmer_from_file(read, kmer_name);

    /// Free memory
    delete read;

    cout << "end" << endl;
    return 0;
}
