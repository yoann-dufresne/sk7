#include "Kff_scanner.hpp"

using namespace std;

/**
 * Initialise a Kff_scanner
 * @param bucketed does the input file when through the bucket tool
 * @param compacted does the input file when through the compact tool
 * @param sorted does the input file when through the compact tool with the -s option
 * @param m for non bucketed files the m is needed for further work
 * @param path path to the kff file
 */
Kff_scanner::Kff_scanner(char *path, bool bucketed, bool compacted, bool sorted, uint8_t m) {

    this->bucketed = bucketed;
    this->compacted = bucketed & compacted;
    this->sorted = compacted & sorted;
    this->file_path = path;

    if (not this->bucketed) {
        this->m = m;
    }

    preparation();

    this->file = new Kff_file(this->file_path, "r");
    uint8_t * metadata = new uint8_t[file->metadata_size + 1];
    file->read_metadata(metadata);
    metadata[file->metadata_size] = '\0';
    delete[] metadata;

    if (file->footer == nullptr) {
        std::cerr << "No footer when one expected!" << std::endl;
        exit(1);
    }
    std::cout << "Footer discovered" << std::endl;

    // --- Global variable read ---
    Section_GV sgv = Section_GV(file);

    this->k = file->global_vars["k"];
    this->m = file->global_vars["m"];
    this->max = file->global_vars["max"];
    this->data_size = file->global_vars["data_size"];

    sk7::initLib(k , m);

}

/**
 * Destructor
 */
Kff_scanner::~Kff_scanner() {
    file->close();
    delete file;
//    std::filesystem::remove_all("tmp/");
}

/// Preparing file

void Kff_scanner::preparation() {

    std::filesystem::create_directory("tmp/");

    if (not bucketed) {
        auto bucket = new Bucket(this->file_path, "tmp/bucket_tmp", this->m);
        bucket->exec();
        delete bucket;
    }

    if (not compacted or not sorted) {
        auto compact = new Compact("tmp/bucket_tmp", "tmp/compact_tmp", true);
        compact->exec();
        delete compact;

        this->file_path = "tmp/compact_tmp";
    }

}


/// Reading compacted and sorted file

/**
 * Read the content of a kff file
 * @return a map of Bucket where the keys are the minimizer value
 */
BucketMap* Kff_scanner::readAll() {

    BucketMap* result = new BucketMap();
    sk7::Bucket bucket;
    std::pair<uint64_t , sk7::Bucket> pair;

    while(file->tellp() < file->end_position) {
        char section_name = file->read_section_type();
        if (section_name == 'm') {
            bucket = readMinimiserSection();
            pair = {bucket.minimiser, bucket};
            result->addBucket(bucket);
        }

        else if (not file->jump_next_section()) {
            break;
        }
    }
    return result;
}


uint64_t decodeMinimiser(uint8_t * encoded, size_t size);
SuperKmer decodeSuperKmer(uint8_t * encoded, size_t size, int64_t mini_pos);
/**
 * Read a minimiser section and build a bucket
 * @return a Bucket containing the SuperKmers in the section
 */
sk7::Bucket Kff_scanner::readMinimiserSection() {

    uint8_t * seq = new uint8_t[(max + k) / 8 + 1];
    memset(seq, 0, (max + k) / 8 + 1);
    uint8_t * data = new uint8_t[max * data_size];
    memset(data, 0, max * data_size);

    Section_Minimizer sm = Section_Minimizer(file);
    uint64_t minimiser = decodeMinimiser(sm.minimizer, m);
    sk7::Bucket result = sk7::Bucket(minimiser);


    for (uint64_t i=0 ; i<sm.nb_blocks ; i++) {
        uint64_t mini_pos;
        uint nb_kmers = sm.read_compacted_sequence_without_mini(seq, data, mini_pos);
        SuperKmer SK = decodeSuperKmer(seq, (nb_kmers + k - sk7::m - 1), (int64_t) mini_pos);
        result.addToList(SK);
    }


    delete[] seq;
    delete[] data;

    return result;
}


/// Decoding functions

uint64_t decodeMinimiser(uint8_t * encoded, size_t size) {
    uint64_t result = 0;

    // Decode the truncated first compacted 8 bits
    size_t remnant = size % 4;
    if (remnant > 0) {
        for (size_t j=0 ; j<remnant ; j++) {
            result = (result << 2) + (0b11 & (encoded[0] >> 2*((3 - j))));
        }
        encoded++;
    }

    // Decode all the 8 bits packed
    size_t nb_uint_used = size / 4;
    for (size_t i=0 ; i<nb_uint_used ; i++) {
        for (size_t j=0 ; j<4 ; j++) {
            result = (result << 2) + (0b11 & (encoded[i] >> (2*(3 - j)))) ;
        }
    }


    return result;
}

SuperKmer decodeSuperKmer(uint8_t * encoded, size_t size, int64_t mini_pos) {

    SuperKmer result = SuperKmer();

    int64_t prefixLen = mini_pos;
    int64_t suffixLen = size - prefixLen;
    int64_t maxLen = fmax(suffixLen, prefixLen);

    int64_t remnant = (4 - ((int64_t) size % 4)) % 4 * 2;

    // Set prefix and suffix
    result.setBits(0, sk7::fixBitSize, prefixLen);
    result.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);

    // Decode the SK
    int64_t i = 0;
    for (; i < maxLen; i++) {


        int64_t index = (i + mini_pos + remnant / 2) / 4;
        uint64_t shift = (4 - ((i + mini_pos + 1 + remnant / 2) % 4)) % 4;


        if (i < suffixLen) {
            result.setBits(2u * sk7::fixBitSize + 4 * i, 2u, ((encoded[index] >> ((shift) * 2)) & 0b11));
        }
        else {
            result.setBits(2u * sk7::fixBitSize + 4 * i, 2u, 0u);
        }

        index = (mini_pos - i - 1 + remnant / 2) / 4;
        shift = (4 - ((mini_pos - i + remnant / 2) % 4)) % 4;
        if (i < prefixLen) {
            result.setBits(2u * sk7::fixBitSize + 4 * (i + 1) - 2, 2u, ((encoded[index] >> (shift * 2)) & 0b11));
        }
        else
            result.setBits(2u * sk7::fixBitSize + 4 * (i + 1) - 2, 2u, 0);
    }

    return result;
}
