//
//  huffmanelementscount.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef HUFFMANELEMENTSCOUNT_H
#define HUFFMANELEMENTSCOUNT_H
#include <iostream>
#include <vector>
using namespace std;

class HuffmanElementsCount
{
public:
    int codeLength;
    vector<int>elementsCodedWithCodeLengthBits;
    HuffmanElementsCount();
};

#endif // HUFFMANELEMENTSCOUNT_H
