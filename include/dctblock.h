//
//  dctblock.h
//  TCM
//
//  Created by Hossam Amer on 2018-08-20.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef DCT_BLOCK_H
#define DCT_BLOCK_H

class DCTBlock
{
public:
    int m_DCT [8][8];	// DCT coeffients for every component on the heap
    DCTBlock(short int & values) {
        for(int i = 0 ; i < 8; ++i) {
            for(int j = 0; j < 8; ++j) {
                m_DCT[i][j] = values[i + j*8];
            }
        }
    }
    
    
    DCTBlock(int values [8][8]) {
        for(int i = 0 ; i < 8; ++i) {
            for(int j = 0; j < 8; ++j) {
                m_DCT[i][j] = values[i][j];
            }
        }
    }
};


#endif /* DCT_BLOCK_H */
