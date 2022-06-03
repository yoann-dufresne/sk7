#include "../../src/headers/Bucket.hpp"

#include <vector>
#include <iostream>

using namespace std;

int main() {

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
    Kmer toSearch3 = Kmer(0b1101000000, 5); //TCAAA
    Kmer toSearch4 = Kmer(0b0000000111, 5); //AAACT;
    Kmer toSearch5 = Kmer(0b1000000011, 5) ; //GAAAT
    Kmer toSearch6 = Kmer(0b1010000000, 5); //GGAAA
    Kmer toSearch7 = Kmer(0b0110000000, 5) ;//CGAAA
    Kmer toSearch8 = Kmer(0b0000001111, 5); //AAATT
    Kmer toSearch9 = Kmer(0b1000000001, 5); //GAAAC
    Kmer toSearch10 = Kmer(0b00000011, 4); //AAAT
    Kmer toSearch11 = Kmer(0b00000010, 4); //AAAC

    cout << "On cherche 1 : " << toSearch.toString() << " : " << bucket.isIn(toSearch) << endl ;
    cout << "On cherche 2 : " << toSearch2.toString() << " : " <<bucket.isIn(toSearch2) << endl;
    cout << "On cherche 3 : " << toSearch3.toString() << " : " << bucket.isIn(toSearch3) << endl;
    cout << "On cherche 4 : " << toSearch4.toString() << " : " << bucket.isIn(toSearch4) << endl;
    cout << "On cherche 5 : " << toSearch5.toString() << " : " <<bucket.isIn(toSearch5) << endl;
    cout << "On cherche 6 : " << toSearch6.toString() << " : " << bucket.isIn(toSearch6) << endl;
    cout << "On cherche 7 : " << toSearch7.toString() << " : " << bucket.isIn(toSearch7) << endl;
    cout << "On cherche 8 : " << toSearch8.toString() << " : " << bucket.isIn(toSearch8) << endl;
    cout << "On cherche 9 : " << toSearch9.toString() << " : " << bucket.isIn(toSearch9) << endl;
    cout << "On cherche 10 : " << toSearch10.toString() << " : " << bucket.isIn(toSearch10) << endl;
    cout << "On cherche 11 : " << toSearch11.toString() << " : " << bucket.isIn(toSearch11) << endl;

    return 0;
}