#include <iostream>
#include "Bucket.hpp"

using namespace std;

int main(int argc, char*argv[]) {

    cout << "I'm main" << endl;
    /// SuperKmer 1 : G|CT
    uint8_t firstSection1 = 0b01100111; //1, 2, C, G
    uint8_t secondSection1 = 0b10000000; //T, _, _, _
    vector<TYPE> tab1 = vector<TYPE>();
    tab1.push_back(firstSection1);
    tab1.push_back(secondSection1);
    SuperKmer SK1(tab1);

    ///SuperKmer 2 : CG|TT
    uint8_t firstSection2 = 0b10101011; //2, 2, T, G
    uint8_t secondSection2 = 0b10010000; // T, C, _, _
    vector<TYPE> tab2 = vector<TYPE>();
    tab2.push_back(firstSection2);
    tab2.push_back(secondSection2);
    SuperKmer SK2(tab2);

    ///SuperKmer 3 : GG|
    uint8_t firstSection3 = 0b10000011; //2, 0, _, G
    uint8_t secondSection3 = 0b00110000;//_, G, _, _
    vector<TYPE> tab3 = vector<TYPE>();
    tab3.push_back(firstSection3);
    tab3.push_back(secondSection3);
    SuperKmer SK3(tab3);


    Bucket bucket(3, 0, 5);
    bucket.addToList(SK1);
    bucket.addToList(SK2);
    bucket.addToList(SK3);

    Kmer toSearch5 = Kmer(0b1100000010, 5) ; // GAAAT
    Kmer toSearch3 = Kmer(0b1001000000, 5); // TCAAA

    bucket.addKmer(toSearch5); // do nothing
    bucket.addKmer(toSearch3); // insert in 0

    int position;
    cout << "after add1 : " << bucket.find(toSearch3, position) << " at : " << position << endl;

    Kmer toSearch2 = Kmer(0b1000000001, 5); // TAAAC
    bucket.addKmer(toSearch2); // insert in 1
    cout << "after add 2: " << bucket.find(toSearch2, position) << " at : " << position << endl;
    cout << "after add 2: " << bucket.find(toSearch3, position) << " at : " << position << endl;

    Kmer toSearch4 = Kmer(0b1100000010, 5); // GAAAT -> TG
    cout << "after add 2: " << bucket.find(toSearch4, position) << " at : " << position << endl;

    Kmer toSearch6 = Kmer(0b1111000000, 5); // GGAAA -> _G_G
    cout << "after add 2: " << bucket.find(toSearch6, position) << " at : " << position << endl;

    return 0;
}
