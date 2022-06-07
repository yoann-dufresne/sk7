#include <iostream>
#include <unordered_map>
#include <cmath>

#include "headers/FastaReader.hpp"
#include "headers/Kmer.hpp"
#include "headers/exampleHash.hpp"
#include "headers/Minimiser.hpp"
#include "headers/utils.hpp"
#include "headers/SuperKmer.hpp"
#include "headers/Bucket.hpp"

using namespace std;

int main(int argc, char*argv[]) {
//    if (argc < 4) {
//        cout << "Usage : ./sk7 file k m" << endl;
//        return 0;
//    }
//
//    int k = atoi(argv[2]);
//    int m = atoi(argv[3]);

//    FILE* fasta = fopen(argv[1], "r");
//    unordered_map<char, int> encoding;
//    encoding['A'] = 0b00;
//    encoding['C'] = 0b01;
//    encoding['G'] = 0b10;
//    encoding['T'] = 0b11;

//    FastaReader fastareader = FastaReader(fasta, false);
//    uint64_t filter = fastareader.set_filter(k);
//
//    uint64_t current_value = 0;
//    int cmpt  = 0;
//    Minimiser minimiser(alpha, m);
//
//    while (true) {
//        char c = fastareader.read_fasta();
//        if (c == EOF) break;
//        cmpt++;
//        current_value = ((current_value << 2) + encoding.at(c)) & filter;
//        if (cmpt == k) {
//            Kmer kmer = Kmer(current_value, k);
//            minimiser.init(kmer);
//            cout << "First kmer : " << kmer.toString() << endl;
//            cout << "test interleaved : " << interleavedOrder(kmer, m).toString() << endl;
//
//        }
//        if (cmpt > k) {
//            Kmer kmer = Kmer(current_value, k);
//            minimiser.fromNewEnd(kmer);
//            cout << "kmer : " << kmer.toString() << endl;
//            cout << "test interleaved : " << interleavedOrder(kmer, m).toString() << endl;
//
//        }
//    }
//
//    fclose(fasta);

    /// SuperKmer 1 : G|CT
    uint8_t firstSection1 = 0b01100110; //1, 2, C, G
    uint8_t secondSection1 = 0b11000000; //T, _, _, _
    vector<TYPE> tab1 = vector<TYPE>();
    tab1.push_back(firstSection1);
    tab1.push_back(secondSection1);
    SuperKmer SK1(tab1);

    ///SuperKmer 2 : CG|TT
    uint8_t firstSection2 = 0b10101110; //2, 2, T, G
    uint8_t secondSection2 = 0b11010000; // T, C, _, _
    vector<TYPE> tab2 = vector<TYPE>();
    tab2.push_back(firstSection2);
    tab2.push_back(secondSection2);
    SuperKmer SK2(tab2);

    ///SuperKmer 3 : GG|
    uint8_t firstSection3 = 0b10000010; //2, 0, _, G
    uint8_t secondSection3 = 0b00100000;//_, G, _, _
    vector<TYPE> tab3 = vector<TYPE>();
    tab3.push_back(firstSection3);
    tab3.push_back(secondSection3);
    SuperKmer SK3(tab3);


    Bucket bucket(3, 0, 5);
    bucket.addToList(SK1);
    bucket.addToList(SK2);
    bucket.addToList(SK3);


    cout << "The bucket contains : " << " G|CT" << " CG|TT" << " GG|" << endl;
    Kmer toSearch = Kmer(0b0000001100, 5); // AAATA
    Kmer toSearch2 = Kmer(0b1100000001, 5); // TAAAC
    Kmer toSearch3 = Kmer(0b1101000000, 5); // TCAAA
    Kmer toSearch4 = Kmer(0b0000000111, 5); // AAACT;
    Kmer toSearch5 = Kmer(0b1000000011, 5) ; // GAAAT
    Kmer toSearch6 = Kmer(0b1010000000, 5); // GGAAA
    Kmer toSearch7 = Kmer(0b0110000000, 5) ;// CGAAA
    Kmer toSearch8 = Kmer(0b0000001111, 5); // AAATT
    Kmer toSearch9 = Kmer(0b1000000001, 5); // GAAAC
    Kmer toSearch10 = Kmer(0b00000011, 4); // AAAT
    Kmer toSearch11 = Kmer(0b00000010, 4); // AAAC
    Kmer toSearch12 = Kmer(0b10000000, 4); // GAAA
    Kmer toSearch13 = Kmer(0b01000000, 4); // CAAA
    Kmer toSearch14 = Kmer(0b11000000, 4); // TAAA
    Kmer toSearch15 = Kmer(0b00000000, 4); // AAAA
    Kmer toSearch16 = Kmer(0b00000001, 4); // AAAC
    Kmer toSearch17 = Kmer(0b011000000011, 6); // CGAAAT
    Kmer toSearch18 = Kmer(0b100000000111, 6); // GAAACT
    Kmer toSearch19 = Kmer(0b100000001111, 6); // GAAATT
    Kmer toSearch20 = Kmer(0b01100000001111, 7); // CGAAATT
    Kmer toSearch21 = Kmer(0b101000000011, 6); // GGAAAT
    Kmer toSearch22 = Kmer(0b100000001100, 6); // GAAATA
    Kmer toSearch23 = Kmer(0b011000000001111, 8); // CGAAAATT
    Kmer toSearch24 = Kmer(0b01100000000011, 7); // CGAAAAT
    Kmer toSearch25 = Kmer(0b01101100000011, 7); // CGTAAAT
    Kmer toSearch26 = Kmer(0b001010000000, 6); //AGGAAA



    cout << "On cherche 1 : " << toSearch.toString() << " : " << bucket.isIn(toSearch) << endl;
    cout << "On cherche 2 : " << toSearch2.toString() << " : " << bucket.isIn(toSearch2) << endl;
    cout << "On cherche 3 : " << toSearch3.toString() << " : " << bucket.isIn(toSearch3) << endl;
    cout << "On cherche 4 : " << toSearch4.toString() << " : " << bucket.isIn(toSearch4) << endl;
    cout << "On cherche 5 : " << toSearch5.toString() << " : " << bucket.isIn(toSearch5) << endl;
    cout << "On cherche 6 : " << toSearch6.toString() << " : " << bucket.isIn(toSearch6) << endl;
    cout << "On cherche 7 : " << toSearch7.toString() << " : " << bucket.isIn(toSearch7) << endl;
    cout << "On cherche 8 : " << toSearch8.toString() << " : " << bucket.isIn(toSearch8) << endl;
    cout << "On cherche 9 : " << toSearch9.toString() << " : " << bucket.isIn(toSearch9) << endl;
    cout << "On cherche 10 : " << toSearch10.toString() << " : " << bucket.isIn(toSearch10) << endl;
    cout << "On cherche 11 : " << toSearch11.toString() << " : " << bucket.isIn(toSearch11) << endl;
    cout << "On cherche 12 : " << toSearch12.toString() << " : " << bucket.isIn(toSearch12) << endl;
    cout << "On cherche 13 : " << toSearch13.toString() << " : " << bucket.isIn(toSearch13) << endl;
    cout << "On cherche 14 : " << toSearch14.toString() << " : " << bucket.isIn(toSearch14) << endl;
    cout << "On cherche 15 : " << toSearch15.toString() << " : " << bucket.isIn(toSearch15) << endl;
    cout << "On cherche 16 : " << toSearch16.toString() << " : " << bucket.isIn(toSearch16) << endl;
    cout << "On cherche 17 : " << toSearch17.toString() << " : " << bucket.isIn(toSearch17) << endl;
    cout << "On cherche 18 : " << toSearch18.toString() << " : " << bucket.isIn(toSearch18) << endl;
    cout << "On cherche 19 : " << toSearch19.toString() << " : " << bucket.isIn(toSearch19) << endl;
    cout << "On cherche 20 : " << toSearch20.toString() << " : " << bucket.isIn(toSearch20) << endl;
    cout << "On cherche 21 : " << toSearch21.toString() << " : " << bucket.isIn(toSearch21) << endl;
    cout << "On cherche 22 : " << toSearch22.toString() << " : " << bucket.isIn(toSearch22) << endl;
    cout << "On cherche 23 : " << toSearch23.toString() << " : " << bucket.isIn(toSearch23) << endl;
    cout << "On cherche 24 : " << toSearch24.toString() << " : " << bucket.isIn(toSearch24) << endl;
    cout << "On cherche 25 : " << toSearch25.toString() << " : " << bucket.isIn(toSearch25) << endl;
    cout << "On cherche 26 : " << toSearch26.toString() << " : " << bucket.isIn(toSearch26) << endl;


    return 0;
}
