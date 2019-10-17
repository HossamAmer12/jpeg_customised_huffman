//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//  component.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Copyright © 2018 Hossam Amer. All rights reserved.
//

#ifndef COMPONENT_H
#define COMPONENT_H
#include "quantizationtable.h"

class Component
{
public:
    int componentID;
    int HFactor; // component horizontal sampling factor
    int VFactor; // component vertical sampling factor
    int HScale, VScale;
    int componentQTableID;
    QuantizationTable *componentQuantizationTable;
//    short int m_DCT[65];	// DCT coeffients for every component
    short int* m_DCT;	// DCT coeffients for every component on the heap
    Component(int id, int HFactor, int VFactor, int QTableID, QuantizationTable &table);
//    ~Component();
};

#endif // COMPONENT_H
