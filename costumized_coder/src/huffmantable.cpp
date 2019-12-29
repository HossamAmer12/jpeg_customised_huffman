//
//  huffmantable.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//
#if IS_RUN_MS_WINDOWS
#include "stdafx.h"
#endif

#include "huffmantable.h"


// Intialize codes with size equal 256 and all ones (8 bits, 2^8 possible values)
HuffmanTable::HuffmanTable(): codes(256, 0xFFFFFFFF), codeLengths(256, 1000), number_of_codes_for_each_1to16(16, 0)
{
    
    // Initial value of 1000 for code length is logic. Why?
    // Because code length can NEVER reach value of 1000. If we put number<0 sometimes this if statement can be true, and we can have problems
    //  if (currentDataLength>=huffmanTable.luminanceDCHuffmanCodeLength[i] && huffmanTable.luminanceDChuffmanCode[i] == data >> (currentDataLength-huffmanTable.luminanceDCHuffmanCodeLength[i])) {    
}


/*void HuffmanTable::deleteUnnecessaryData() {
 for (int i=0; i<tables.size(); i++)
 delete tables[i];
 tables.clear();
 }*/
