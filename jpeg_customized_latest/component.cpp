//
//  component.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#if IS_RUN_MS_WINDOWS
#include "stdafx.h"
#endif


#include "component.h"

Component::Component(int id, int HFactor, int VFactor, int QTableID, QuantizationTable &table) : componentQuantizationTable(&table)
{
    componentQTableID = QTableID;
    this->componentID = id;
    this->HFactor = HFactor;
    this->VFactor = VFactor;
    m_DCT = new short int [65];
}


// Destructor
//Component::~Component() {
//    
//    if(m_DCT != nullptr)
//    {
//        delete [] m_DCT;
//        m_DCT = nullptr;
//    }
//}
