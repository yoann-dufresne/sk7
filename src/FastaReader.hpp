#ifndef SUPERKMERS_FASTAREADER_HPP
#define SUPERKMERS_FASTAREADER_HPP


class FastaReader {
private:
    FILE* file;
    void set_fasta();

public:
    FastaReader(FILE* file, bool header);
    char read_fasta();
    uint64_t set_filter(int k);
};


#endif //SUPERKMERS_FASTAREADER_HPP
