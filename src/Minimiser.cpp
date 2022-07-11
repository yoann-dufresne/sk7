#include "Minimiser.hpp"

/**
 * Constructor for the minimiser of a given kmer
 * @param hashFunction a function that gives an value to a mmer
 * @param length the wanted length of the minimiser
 * @param kmer the initial kmer
 */
Minimiser::Minimiser(hashPos (*hashFunction)(Kmer, ushort), ushort length, Kmer kmer): Minimiser(hashFunction, length) {
    init(kmer);
}

/**
 * Constructor for a Minimiser taking reverse complement in consideration
 * @param hashFunction a function that gives an value to a mmer
 * @param kmer the initial kmer
 */
Minimiser::Minimiser(hashPos (*hashFunction)(Kmer, ushort), Kmer &kmer) {

    this->hashFunction = hashFunction;
    this->length = sk7::m;

    Kmer reverse = kmer.reverseComplement();

    hashPos hash = hashFunction(kmer, sk7::m);
    hashPos hashRev = hashFunction(reverse, sk7::m);

    if (hash.hashValue < hashRev.hashValue) {
        this->value = hash.hashValue;
        this->pos = hash.pos;
        return;
    }
    if (hash.hashValue > hashRev.hashValue) {
        this->value = hashRev.hashValue;
        this->pos = hashRev.pos;
        kmer = reverse;
        return;
    }

    // same minimiser value
    int dist = abs(hash.pos - (sk7::k - 1) / 2);
    int distRev = abs(hashRev.pos - (sk7::k - 1) / 2);
    if (dist < distRev) {
        this->value = hash.hashValue;
        this->pos = hash.pos;
        return;
    }
    if (dist > distRev) {
        this->value = hashRev.hashValue;
        this->pos = hashRev.pos;
        kmer = reverse;
        return;
    }

    // position and value equivalent
    Kmer withoutMinimiser = kmer.removePart(hash.pos, sk7::m);
    interleavedOrder(withoutMinimiser, hash.pos);
    Kmer withoutMinimiserRev = reverse.removePart(hashRev.pos, sk7::m);
    interleavedOrder(withoutMinimiserRev, hash.pos);
    if (withoutMinimiser.getValue() < withoutMinimiserRev.getValue()) {
        this->value = hash.hashValue;
        this->pos = hash.pos;
        return;
    }

    if(withoutMinimiser.getValue() > withoutMinimiserRev.getValue()) {
        this->value = hashRev.hashValue;
        this->pos = hashRev.pos;
        kmer = reverse;
        return;
    }


    this->value = hash.hashValue;
    this->pos = hash.pos;

}


/**
 * Constructor for a minimiser for a sequence of kmer
 * @param hashFunction a function that gives an whole value to a mmer
 * @param length the wanted length of the minimiser
 */
Minimiser::Minimiser(hashPos (*hashFunction)(Kmer, ushort), ushort length) {
    this->hashFunction = hashFunction;
    this->length = length;
    this->value = 0;
    this->pos = 0;
}

/**
 * Initialise value with the hash of the given kmer
 * @param kmer the initializer
 */
void Minimiser::init(Kmer kmer) {
    hashPos retHash = hashFunction(kmer, this->length);
    this->value = retHash.hashValue;
    this->pos = retHash.pos;
}

/**
 * Found the new minimiser depending on the new end
 * @param kmer the current kmer to minimise
 */
void Minimiser::fromNewEnd(Kmer kmer) {
    if (pos == 0) { //Forced change of the minimiser
        this->init(kmer);
        return;
    }
    Kmer end = kmer.getSubKmer(kmer.getLength() - length, kmer.getLength() - 1);
    hashPos retHash = hashFunction(end, length);
    uint64_t end_val = retHash.hashValue;
    if (end_val < this->value) { //The new ending is the minimiser
        this->value = end_val;
        this->pos = kmer.getLength() - length;
    } else {
        this->pos--; //We keep the same minimiser but its position shift
    }
}

/**
 * Getter for the minimiser value
 * @return the minimiser's value
 */
uint64_t Minimiser::getValue() const {
    return this->value;
}

/**
 * Getter for the pos attribute
 * @return the minimiser's pos
 */
short Minimiser::getPos() const {
    return this->pos;
}

/**
 * Give a textual representation of a minimiser
 * @return a string on {'A', 'C', 'G', 'T'} that represents the minimiser
 */
std::string Minimiser::toString() const {
    uint64_t tmp = value;
    std::string res;
    res.reserve(this->length);
    ushort b = 4;
    while (tmp != 0) {
        uint64_t q = tmp / b;
        int r = tmp % b;
        switch (r) {
            case 0:
                res = 'A' + res;
                break;
            case 1:
                res = 'C' + res;
                break;
            case 2:
                res = 'T' + res;
                break;
            default:
                res = 'G' + res;
                break;
        }
        tmp = q;
    }
    while (res.length() < this->length) {
        res = "A" + res;
    }
    return res;
}



