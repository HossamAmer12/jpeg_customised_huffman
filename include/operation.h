// Created by Yanbing Jiang

//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//
// operation.cpp
// operation includes:
// 1. Smoothing
// 2. Expansion
// 3. Content Adaptive DCT

#include "TypeDef.h"

#include <iostream>
#include <math.h>
#include <vector>      // std::vector
#include <algorithm>
using namespace std;

// DCT Operation
float C_dct(int u) {
	if (u == 0){
		return 1.0f / sqrt(8.0f);
	}else{
		return sqrt(2.0 / 8.0);
	}
}


double func_dct(int x, int y, const int block[8][8]) {
	const float PI = 3.14f;
	float sum = 0;
	for (int u = 0; u < 8; ++u){
		for (int v = 0; v < 8; ++v){
			sum += (block[u][v] * cosf(((2 * v + 1) * x * PI) / 16)  * cosf(((2 * u + 1) * y * PI) / 16));
		} // end inner loop
	} // end outer loop
	sum = C_dct(x) * C_dct(y) * sum;
	return sum;
}

void perform_dct(vector<vector<double> > &outBlock, int inBlock[8][8]) {
	for (int y = 0; y < 8; ++y){
		for (int x = 0; x < 8; ++x){
			outBlock[y][x] = func_dct(x, y, inBlock);
		} // end inner loop
	} // end outer loop

}
// Generate the test 8*8 for every iteration
void generate_block(int** image, int positionH, int positionV , int (&block)[8][8]){
	for(int m = positionH; m<=positionH+7; m++){
    	for(int n = positionV; n<=positionV+7; n++){
    		block[m-positionH][n-positionV] = image[m][n];
    	}
    }
}

// block value set
void block_value_set(vector<vector <int> > &tCoeff_Y, int positionH, int positionV, double value){
	int inblock[8][8] = { 0 };
	vector<vector <double> > outblock;
	outblock.resize(8, vector<double>(8, 0));
	for(int row = 0; row<8; row++){
		for(int col = 0; col<8; col++){
			inblock[row][col] = value;
		}
	}
	perform_dct(outblock, inblock);
	for (int row = positionH; row<positionH + 8; row++) {
		for (int col = positionV; col<positionV + 8; col++) {
			tCoeff_Y[row][col] = outblock[row- positionH][col- positionV];
		}
	}
}
// test whether a block is a black block
bool zero_block_test(const vector<int> count_outlier_list, int block_idx) {
	/*for (int test_row = 0; test_row<8; test_row++) {
		for (int test_col = 0; test_col<8; test_col++) {
			if (region[test_row][test_col] != 0) {
				return false;
			}
		}
	}
	return true;*/

	if (count_outlier_list.at(block_idx) <= TCM_OUTLIER_THRESHOLD) {
		return true;
	}
	else {
		return false;
	}
}

// test whether the region we are scanning is all zeros
void sum_addition(const vector<int> count_outlier_list, int block_idx, int (&sum)[8][8], int& counter, int (&region)[8][8]){
	if(zero_block_test(count_outlier_list, block_idx)){
		counter++; //counter that count the number of black blocks around the test block
	}else{
		for(int test_row = 0; test_row<8; test_row++){
			for(int test_col = 0; test_col<8; test_col++){
				sum[test_row][test_col] += region[test_row][test_col];
			}
		}
	}

	//Recover the region to zero
	for (int test_row = 0; test_row<8; test_row++) {
		for (int test_col = 0; test_col<8; test_col++) {
			region[test_row][test_col] = 0;
		}
	}
}



// Smoothing Operation
int smoothing(int block_idx, const vector<int> count_outlier_list, int height, int width, int currentComponent, unsigned char** m_YPicture_buffer, 
	unsigned char** m_CbPicture_buffer, unsigned char** m_CrPicture_buffer) {
	
	unsigned char** m_buffer;
	switch (currentComponent) {
	case COMPONENT_Y:
		m_buffer = m_YPicture_buffer;
		break;
	case COMPONENT_Cb:
		m_buffer = m_CbPicture_buffer;
		break;
	case COMPONENT_Cr:
		m_buffer = m_CrPicture_buffer;
		break;
	default:
		m_buffer = m_YPicture_buffer;
		break;
	} // end switch
	
	// Remove Noise for normal image
	int number_of_block_row = height / 8;
	int number_of_block_col = width / 8;
	int pos_of_block_row = block_idx / number_of_block_col;
	int pos_of_block_col = ceil(block_idx % number_of_block_col);

	int NW, N, NE, W, E, SW, S, SE;

	// NW
	if (block_idx < number_of_block_col | block_idx % number_of_block_col == 0) {
		NW = block_idx;
	}
	else {
		NW = block_idx - 1 - number_of_block_col;
	}

	// N
	if (block_idx < number_of_block_col) {
		N = block_idx;
	}
	else {
		N = block_idx - number_of_block_col;
	}

	// NE
	if (block_idx < number_of_block_col | (block_idx + 1) % number_of_block_col == 0) {
		NE = block_idx;
	}
	else {
		NE = block_idx + 1 - number_of_block_col;
	}

	// W
	if (block_idx % number_of_block_col == 0) {
		W = block_idx;
	}
	else {
		W = block_idx - 1;
	}

	// E
	if ((block_idx + 1) % number_of_block_col == 0) {
		E = block_idx;
	}
	else {
		E = block_idx + 1;
	}

	// SW
	if (block_idx >= (number_of_block_row - 1)*number_of_block_col || block_idx % number_of_block_col == 0) {
		SW = block_idx;
	}
	else {
		SW = block_idx - 1 + number_of_block_col;
	}

	//S
	if (block_idx >= (number_of_block_row - 1)*number_of_block_col) {
		S = block_idx;
	}
	else {
		S = block_idx + number_of_block_col;
	}

	//SE
	if (block_idx >= (number_of_block_row - 1)*number_of_block_col || (block_idx + 1) % number_of_block_col == 0) {
		SE = block_idx;
	}
	else {
		SE = block_idx + 1 + number_of_block_col;
	}

	int counter = 0;
	int sum[8][8] = { 0 };
	int current[8][8] = { 0 };

	// CURRENT
	for (int m = pos_of_block_row; m <= pos_of_block_row + 7; m++) {
		for (int n = pos_of_block_col; n <= pos_of_block_col + 7; n++) {
			current[m - pos_of_block_row][n - pos_of_block_col] = m_buffer[m][n];
		}
	}

	// NW
	int pos_of_block_row_NW = NW / number_of_block_col;
	int pos_of_block_col_NW = ceil(NW % number_of_block_col);
	int block[8][8] = { 0 };
	for (int m = pos_of_block_row_NW; m <= pos_of_block_row_NW + 7; m++) {
		for (int n = pos_of_block_col_NW; n <= pos_of_block_col_NW + 7; n++) {
			block[m - pos_of_block_row_NW][n - pos_of_block_col_NW] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, NW, sum, counter, block);

	// N
	int pos_of_block_row_N = N / number_of_block_col;
	int pos_of_block_col_N = ceil(N % number_of_block_col);
	for (int m = pos_of_block_row_N; m <= pos_of_block_row_N + 7; m++) {
		for (int n = pos_of_block_col_N; n <= pos_of_block_col_N + 7; n++) {
			block[m - pos_of_block_row_N][n - pos_of_block_col_N] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, N, sum, counter, block);

	// NE
	int pos_of_block_row_NE = NE / number_of_block_col;
	int pos_of_block_col_NE = ceil(NE % number_of_block_col);
	for (int m = pos_of_block_row_NE; m <= pos_of_block_row_NE + 7; m++) {
		for (int n = pos_of_block_col_NE; n <= pos_of_block_col_NE + 7; n++) {
			block[m - pos_of_block_row_NE][n - pos_of_block_col_NE] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, NE, sum, counter, block);

	// W
	int pos_of_block_row_W = W / number_of_block_col;
	int pos_of_block_col_W = ceil(W % number_of_block_col);
	for (int m = pos_of_block_row_W; m <= pos_of_block_row_W + 7; m++) {
		for (int n = pos_of_block_col_W; n <= pos_of_block_col_W + 7; n++) {
			block[m - pos_of_block_row_W][n - pos_of_block_col_W] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, W, sum, counter, block);

	// E
	int pos_of_block_row_E = E / number_of_block_col;
	int pos_of_block_col_E = ceil(E % number_of_block_col);
	for (int m = pos_of_block_row_E; m <= pos_of_block_row_E + 7; m++) {
		for (int n = pos_of_block_col_E; n <= pos_of_block_col_E + 7; n++) {
			block[m - pos_of_block_row_E][n - pos_of_block_col_E] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, E, sum, counter, block);

	// SW
	int pos_of_block_row_SW = SW / number_of_block_col;
	int pos_of_block_col_SW = ceil(SW % number_of_block_col);
	for (int m = pos_of_block_row_SW; m <= pos_of_block_row_SW + 7; m++) {
		for (int n = pos_of_block_col_SW; n <= pos_of_block_col_SW + 7; n++) {
			block[m - pos_of_block_row_SW][n - pos_of_block_col_SW] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, SW, sum, counter, block);

	// S
	int pos_of_block_row_S = S / number_of_block_col;
	int pos_of_block_col_S = ceil(S % number_of_block_col);
	for (int m = pos_of_block_row_S; m <= pos_of_block_row_S + 7; m++) {
		for (int n = pos_of_block_col_S; n <= pos_of_block_col_S + 7; n++) {
			block[m - pos_of_block_row_S][n - pos_of_block_col_S] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, S, sum, counter, block);

	// SE
	int pos_of_block_row_SE = SE / number_of_block_col;
	int pos_of_block_col_SE = ceil(SE % number_of_block_col);
	for (int m = pos_of_block_row_SE; m <= pos_of_block_row_SE + 7; m++) {
		for (int n = pos_of_block_col_SE; n <= pos_of_block_col_SE + 7; n++) {
			//cout << "m = " << m << ", n= " << n<< endl;
			block[m - pos_of_block_row_SE][n - pos_of_block_col_SE] = m_buffer[m][n];
		}
	}
	sum_addition(count_outlier_list, SE, sum, counter, block);
	
    // if current block is black block and less than surrounding block is black
    if (zero_block_test(count_outlier_list, block_idx) && counter <= 4){// & err >= threshold){
        // calculate mean
        int avg_sum = 0;
        for(int m = 0; m<8; m++){
            for(int n =0; n<8; n++){
            	avg_sum += sum[m][n];
            }
        }
        avg_sum = avg_sum / 512;

		return avg_sum * 8;
        // replace with the average value
		//block_value_set(tCoeff_Y, pos_of_block_row, pos_of_block_col, avg_sum);
    }

    // more than 4 black blocks around
    if(counter > 4){// & err >= threshold){
        // for(int test_row = i; test_row<i+8; test_row++){
        // 	for(int test_col = j; test_col<j+8; test_col++){
        // 		image[test_row][test_col] = 0;
        // 	}
        // }
		//block_value_set(tCoeff_Y, pos_of_block_row, pos_of_block_col, 0);
		return 0;
    }

    //if(zero_block_test(current) && counter <=4 && err >= threshold){
    //    // expansion(image, m, n, i, j, N, S, W, E);
    //}

    return 0;
}

void perform_content_adaptive_dct(int currentX_topleft, int currentY_topleft, int currentComponent, vector<vector <int> > &tCoeff_Y, vector<vector <int> > &tCoeff_Cb, vector<vector <int> > &tCoeff_Cr) {
	int cdc_start_location = 0;
	switch (currentComponent) {
	case COMPONENT_Y:
		cdc_start_location = 3;
		break;
	case COMPONENT_Cb:
		cdc_start_location = 1;
		break;
	case COMPONENT_Cr:
		cdc_start_location = 1;
		break;
	default:
		cdc_start_location = 3;
		break;

	} // end switch
	 //int cdc_start_location = 3;
	 for (int i = cdc_start_location; i < 8; ++i) {
		 for (int j = cdc_start_location; j < 8; ++j) {
			 switch (currentComponent) {
			 case COMPONENT_Y:
				 tCoeff_Y[currentY_topleft + j][currentX_topleft + i] = 0;
				 break;
			 case COMPONENT_Cb:
				 tCoeff_Cb[currentY_topleft + j][currentX_topleft + i] = 0;
				 break;
			 case COMPONENT_Cr:
				 tCoeff_Cr[currentY_topleft + j][currentX_topleft + i] = 0;
				 break;
			 default:
				 tCoeff_Y[currentY_topleft + j][currentX_topleft + i] = 0;
				 break;

			 } // end switch
		 } // end inner loop
	}// end outer loop
	
} // end method

// Expansion Technique that applied inside the smoothing technique
//void expansion(int** image, int rows, int cols, int current_row, int current_col, int N, int S, int W, int E){
//	int N2 = max(0, current_col-8);
//    int S2 = min(cols-8, current_col+8);
//    int W2 = max(0, current_row-8);
//    int E2 = min(rows-8, current_row+8);
//
//    // 3 of NW
//    int block[8][8]={0};
//    generate_block(image, W2, N2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, W, N2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, W2, N, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, W, N, 0);
//    		}
//    	}
//    }
//
//    // 3 of N
//    block[8][8]={0};
//    generate_block(image, current_row, N2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, W, N2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, E, N2, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, current_row, N, 0);
//    		}
//    	}
//    }
//
//    // 3 of NE
//    block[8][8]={0};
//    generate_block(image, E2, N2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, E, N2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, E2, N, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, E, N, 0);
//    		}
//    	}
//    }
//
//    // 3 of W
//    block[8][8]={0};
//    generate_block(image, W2, current_col, block);
//    if(zero_block_test(block)){
//    	generate_block(image, W2, N, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, W2, S, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, W, current_col, 0);
//    		}
//    	}
//    }
//
//    // 3 of E
//    block[8][8]={0};
//    generate_block(image, E2, current_col, block);
//    if(zero_block_test(block)){
//    	generate_block(image, E2, N, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, E2, S, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, E, current_col, 0);
//    		}
//    	}
//    }
//
//    // 3 of SW
//    block[8][8]={0};
//    generate_block(image, W2, S2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, W, S2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, W2, S, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, W, S, 0);
//    		}
//    	}
//    }
//
//    // 3 of S
//    block[8][8]={0};
//    generate_block(image, current_row, S2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, W, S2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, E, S2, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, current_row, S, 0);
//    		}
//    	}
//    }
//
//    // 3 of SE
//    block[8][8]={0};
//    generate_block(image, E2, S2, block);
//    if(zero_block_test(block)){
//    	generate_block(image, E, S2, block);
//    	if(zero_block_test(block)){
//    		generate_block(image, E2, S, block);
//    		if(zero_block_test(block)){
//    			block_value_set(image, E, S, 0);
//    		}
//    	}
//    }
//}
