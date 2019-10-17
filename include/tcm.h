//
//  tcm.h
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef TCM_H
#define TCM_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;

//#define MaxAmp 65536
#define MaxAmp 256
#define LambdaDelta 0.1
#define MinLikelyhood  -1.e30
//(-(1<<30))
#define StartPointProb 0.1

class tcm {
    
    typedef struct {
        int count;
        double AcumAbsAmp;
        double AcumSampNum;
        double prob;
        double lambda;
        double likelyhood;
    } tBucket;
    
public:
    
//    static tBucket test[10];

    // TCM
    static double ComputeLambdaGivenYc(double Yc, double sumYi, double totalNumYi)
    {
        double c = sumYi/totalNumYi ;
        double lambda, lambda_old ;
        
        if( c/Yc>=0.95) return -1.0 ; // Too much away from Laplacian
        
        lambda_old = c ;
        
        lambda = c - Yc * (1.0 - 1.0/(1.0 - exp(-Yc/lambda_old)) ) ;
        
        { int k;
            for(k=0;k<10;k++)
            {
                lambda_old = lambda ;
                
                lambda = c - Yc * (1.0 - 1.0/(1.0 - exp(-Yc/lambda_old)) ) ;
            }
        }
        while( fabs(lambda - lambda_old) > LambdaDelta )
        {
            lambda_old = lambda ;
            
            lambda = c - Yc * (1.0 - 1.0/(1.0 - exp(-Yc/lambda_old)) ) ;
        }
        
        return lambda ;
        
    }
    
    
    static int FindStartPoint(tBucket *buck, int peak, int N)
    {
        int k ;
        for(k=peak;k>0;k--)
        {
            if( buck[k].count==0) continue ;
            
            if( buck[k].AcumSampNum < N* (1.0-StartPointProb)) break ;
        }
        
        if(buck[0].count>N/100 && buck[1].count>N/100 && buck[2].count>N/100  && buck[3].count>N/100)
        {
            if(k<3) k=3;
        }
        else if(buck[0].count>N/100 && buck[1].count>N/100 && buck[2].count>N/100 )
        {
            if(k<2) k=2;
        }
        else
        {
            if(k<1) k=1;
        }
        return k;
    }
    
    static void  ComputeLikelyhood(int thePoint, int N, tBucket *buck, int peak)
//    static void  ComputeLikelyhood(int thePoint, int N, tBucket buck [MaxAmp], int peak)
    {
        
        double N1 = buck[thePoint].AcumSampNum ;
        
        double N2 = N - N1 ;
        double Yc = thePoint ;
        
        double sumYi = buck[thePoint].AcumAbsAmp ;
        double totalNumYi = buck[thePoint].AcumSampNum ;
        
        double lambda =  ComputeLambdaGivenYc(Yc,  sumYi, totalNumYi) ;
        
        double prob = (double)N1 / (double)N ;
        
        if(lambda>0)
        {
            buck[thePoint].likelyhood = N2 * log(1- prob) + N1 * log(prob) - N2 * log( (peak - Yc)*2.0)
            -N1 * log(1-exp(-Yc/lambda))
            -N1 * log(2*lambda) - sumYi/lambda ;
            
            buck[thePoint].lambda = lambda ;
            buck[thePoint].prob = prob ;
        }
        else
        {
            buck[thePoint].likelyhood = -MinLikelyhood ;
            
            buck[thePoint].lambda = lambda ;
            buck[thePoint].prob = 1 ;
        }
    }
    
//    tBucket test[MaxAmp];
//    tBucket tcm::test[10];
    static void TCMprocessOneSequence(const vector<int> &C, int Len, int *peak, double *prob, double *lambda, double *Yc)
    {
        int k, MaxPos, StartPoint;
        double MaxLikelyhood, likelyhood ;
        tBucket buck[MaxAmp];
        
//        test[0].count = 0;
        
        *peak = 0 ;
        for(k=0;k<Len;k++)
        {
            int absv = abs(C[k]);
            if(absv> (*peak) ) *peak = absv ;
        }
        
        if( *peak ==0 || *peak>=MaxAmp  )
        {
            
            //	        printf("Input data either all zero or exceed %d (%d)\n", MaxAmp, *peak);
            
            *lambda = -1 ;
            *Yc = 0;
            *prob = 1 ;
            return ;
        }
        
        for(k=0;k<= (*peak);k++)
        {
            buck[ k ].count = 0 ;
            buck[ k ].AcumAbsAmp = 0 ;
            buck[ k ].AcumSampNum = 0 ;
        }
        
        for(k=0;k<Len;k++)
        {
            int cabs = abs(C[k]);
            buck[ cabs ].count ++ ;
        }
        buck[0].AcumSampNum = buck[0].count ;
        
        for(k=1;k<= (*peak);k++)
        {
            buck[ k ].AcumAbsAmp  = buck[ k-1 ].AcumAbsAmp + k * buck[ k ].count ;
            buck[ k ].AcumSampNum = buck[ k-1 ].AcumSampNum + buck[ k ].count ;
        }
        
        // saftey check
        if(buck[0].count < (buck[1].count>>1) || buck[1].count < (buck[2].count>>1))
        {
            *lambda = -1 ;
            *Yc = 0;
            *prob = 1 ;
            return ;
        }
        
        StartPoint = FindStartPoint(buck, *peak, Len);
        
        //    if( StartPoint<=0 )
        //    {
        //        *prob = 1 ;
        //        *lambda = -1 ;
        //        *Yc = 0 ;
        //        printf("Data distribution is not like Laplacian1\n");
        //        return ;
        //    }
        
//        ComputeLikelyhood(StartPoint, Len, buck, *peak)  ;
        ComputeLikelyhood(StartPoint, Len, buck, *peak)  ;
        
        MaxLikelyhood = buck[StartPoint].likelyhood ;
        
        MaxPos = StartPoint ;
        for(k=StartPoint+1; k<=(*peak); k++)
        {
            if(buck[k].count==0) continue ;
            
            ComputeLikelyhood(k, Len, buck, *peak)  ;
            
            likelyhood = buck[k].likelyhood ;
            
            if( likelyhood > MaxLikelyhood)
            {
                MaxPos = k ;
                MaxLikelyhood = likelyhood ;
            }
        }
        
        if( MaxLikelyhood > MinLikelyhood)
        {
            *prob = buck[MaxPos].prob ;
            *lambda = buck[MaxPos].lambda ;
            *Yc = MaxPos ;
        }
        else  // the search failed
        {
            *prob = 1 ;
            *lambda = -1 ;
            *Yc = 0 ;
            // QQQQ printf("Data distribution is not like Laplacian2\n");
        }
    } // end TCMprocessOneSequence
    
    
    static void count_outliers(vector<double> const & yc_array, const int total_block, const vector<vector<int>> &tCoeff_Y_AC, vector<int>& count_outlier_list ) {
        
        // Attempt of Variable d(shift) according to the magnitude of Yc
        /*vector<int> shift;
         vector<int> move;
         for (int ac_coeff = 1; ac_coeff < yc_array.size() - 1; ++ac_coeff) {
         for (int blk_index = 0; blk_index < total_block; ++blk_index) {
         move.push_back(abs(tCoeff_Y_AC[ac_coeff][blk_index] - yc_array.at(ac_coeff)));
         }
         sort(move.begin(), move.end());
         shift.push_back(move[62]);
         move.clear();
         }*/
        
        for (int blk_index = 0; blk_index < total_block; ++blk_index) {
            /*vector<int> move;
             for (int ac_coeff = 1; ac_coeff <= yc_array.size() - 1; ++ac_coeff) {
             move.push_back(abs(tCoeff_Y_AC[ac_coeff][blk_index] - yc_array.at(ac_coeff)));
             }
             sort(move.begin(), move.end());
             int shift= move[61];
             move.clear();
             */
            
            int outlier_flag = 0;
            int shift = 30;
            
            // Compare with yc score with accure block AC value and count in outlier_flag and store in the count_outlier_list
            for (int ac_coeff = 1; ac_coeff <= yc_array.size()-1; ++ac_coeff) {
                if (yc_array.at(ac_coeff)+ shift < tCoeff_Y_AC[ac_coeff][blk_index] || -yc_array.at(ac_coeff) - shift > tCoeff_Y_AC[ac_coeff][blk_index]) {
                    outlier_flag++;
                }
            }
            count_outlier_list.at(blk_index) = outlier_flag;
        }
        
    }// end count_outliers
    
    
};


#endif /* tcm_h */
