//
//  quantizationtable.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef QUANTIZATIONTABLE_H
#define QUANTIZATIONTABLE_H

#include <vector>

class QuantizationTable
{
    
public:
    int tableID;
	int tableLength;
    int quantizationTableData[8][8];

    // component if true, then luminance; if false, then chrominance
    QuantizationTable(bool writeFileProcess, bool component);
    
};

#endif // QUANTIZATIONTABLE_H
