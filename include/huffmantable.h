//
//  huffmantable.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef HUFFMANTABLE_H
#define HUFFMANTABLE_H

#include <map>
#include <vector>
#include "TypeDef.h"

using namespace std;

class HuffmanTable {
    
public:
    
    // ID for which Huffman table it is using (there are 32 possible Huffman tables in JPEG)
    unsigned char tableID;
    
    // table length from bitstream
    unsigned short tableSegmentLengthFromBitstream;
    
    // tableClass 0 is for DC, 1 is for AC
    unsigned char  tableClass;
    vector <unsigned int> codes; //ehufco
    vector <unsigned int> codeLengths; //ehufsi
    
    vector <unsigned char> number_of_codes_for_each_1to16;   //bits
    
    // The array of Huffman maps: (length, code) -> value
    std::map<huffKey, unsigned char> huffData;   //( ehufsi , ehufco ) huffval[256]
    
    
    HuffmanTable();
};

#endif // HUFFMANTABLE_H

//typedef struct
//{
//    unsigned char bits[17];  /* bits[k] stores the number of symbols whose codes
//                              are k bits long. bits[0] is unused */
//    unsigned char huffval[256]; /* The symbols in order of increasing code length */
//
//    /* encode tables */
//    unsigned int ehufco[256]; /* Huffman code for each symbol */
//    char ehufsi[256];         /* length of each code */
//
//    /* decoding tables (1st element of each array are unused */
//    unsigned int mincode[17];
//
//    /* NOTE: there is an error in IJG Huffman coding implementation.
//     When set htbl->maxcode[17] = 0xFFFFFL function jdhuff_tbl(), it actually changes the value of next member, valptr[].
//     To fix it, either change long maxcode[17] to long maxcode[18] or change the setting statement to htbl->maxcode[16] = 0xFFFFFL
//     */
//    long maxcode[18];
//    short valptr[17]; /*index of 1st symbol of different length in huffval[] */
//}HUFF_TBL;
