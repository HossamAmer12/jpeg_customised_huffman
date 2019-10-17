//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.

// chi_square.h
// Used for checking the image similar based on histogram by calculating the chi-squared distance
#include <math.h>
#include <vector>      // std::vector
#include <algorithm>
#include <stdlib.h>
using namespace std;

// Calculate the Chi-Squared Distance for image similarity check
double distChiSq (int X[], int Y[]){
	int m = sizeof(X)/sizeof(X[0]);
	int n = sizeof(Y)/sizeof(Y[0]);	// size of the histm, should be 256

	double upper,bottom; // numeraror and denominator in single element
	double ele; //every single element in the summation
	double D; //summation

	for(int i = 0; i < n; i++){
		upper = X[i] - Y[i];
		bottom = X[i] + Y[i];
		ele = pow(s,2)/d;

		D += ele;
	}

	D /=2;
	return D;
}

void imageHist(int img2D[][], int& hist[]){
	// Change 2D image array to a 1D array for easy count purpose
	// Size of the 2D image (row*column)
	int row = sizeof (img2D) / sizeof(img2D[0]);
	int col = sizeof(img2D[0]) / sizeof(img2D[0][0]);

	// Initialization of img 1D array
	int img[row*col];
	for(int i = 0; i < row; ++i)
    	for(int j = 0; j < col; ++j)
    		img[i * col + j] = img2D[i][j];

    // sort 1D image array
    int size_img = sizeof(img)/sizeof(img[0]);
    int max_img = max_element(img, img + size_img);
    sort(img, img + size_img);
    int place = img[0];

    for(int i = 0; i < size_img; i++){
    	if (img[i] == place){
    		hist[place]++
    	}else{
    		place = img[i];
    		hist[place]++;
    	}
    }
    // Done count the frequency of the pixel value.
}