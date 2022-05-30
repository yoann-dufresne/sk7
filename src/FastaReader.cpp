#include <cstdio>
#include <cstdlib>
#include <bitset>
#include <string>

#include "headers/FastaReader.hpp"

/**
 * Construct a reader for the giver file.
 */
FastaReader::FastaReader(FILE* file, bool header)  {
    this->file = file;
    if (header) set_fasta();
}

/**
 * Set the file pointer to the start of the DNA sequence.
 */
void FastaReader::set_fasta(){
    int c;
    c = fgetc(file);
    while(c != EOF && c != '\n') {
        c = fgetc(file);
    }
    if (c == EOF) {
        perror("Bad file format");
        exit(1);
    }
}

/**
 * Read the next character in a fasta file excluding end of line and 'N' characters.
 * @return the next character of the sequence.
 */
char FastaReader::read_fasta() {
    int c = fgetc(file);
    while (c == '\n' || c == 'N') c = fgetc(file);
    return (char) c;
}

uint64_t FastaReader::set_filter(int k) {
    std::string str;
    str.resize(2*(k+1), '1');
    std::fill(str.begin(), str.begin() + 2, '0');
    std::bitset<64> filtre(str);
    return filtre.to_ullong();
}