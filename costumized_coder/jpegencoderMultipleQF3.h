//
//  jpegencoderMultipleQF.h
//  jpeg_tcm
//
//  Created by Hossam Amer on 2018-11-04.
//  Copyright :copyright: 2018 Hossam Amer. All rights reserved.
//

#ifndef JPEG_ENCODER_MULTIPLE_QF_H
#define JPEG_ENCODER_MULTIPLE_QF_H
#define Enable_MultiEncoder 0
#define Enable_MultiDecoder 1
#include <thread>
#include <mutex>
#include "jpegdecoder.h"
// #include "jpegencoder.h"
#define S7S_debug 0
static double getPSNR(const cv::Mat & i1, const cv::Mat & i2, int comp_width, int comp_height)
{
    
    int d = CV_8U;
    Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);
    
    
    Mat s1;
    absdiff(I1, I2, s1);       // |I1 - I2|
#if S7S_debug
    cout << "here **" << endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            cout << int(s1.at<unsigned char>(i, j)) << "," << int(I1.at<unsigned char>(i, j)) << "," << int(I2.at<unsigned char>(i, j)) << " ";
        }
        cout << "\n";
    }
#endif
    s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
    s1 = s1.mul(s1);           // |I1 - I2|^2
    
    Scalar s = sum(s1);         // sum elements per channel
    
    double mse = (s[0]); // sum channels
#if S7S_debug
    cout << "MSE: " << mse << endl;
#endif
    if (mse <= 1e-10) // for small values return zero
    {
        return 999.99;
    }
    else
    {
        double psnr = 10.0 * log10(double(255 * 255 * double(comp_width * comp_height)) / mse);
        return psnr;
    }
}

Scalar getMSSIM(const Mat& i1, const Mat& i2)
{
    const double C1 = 6.5025, C2 = 58.5225;
    /***************************** INITS **********************************/
    int d = CV_32F;
    
    Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);
    
    Mat I2_2 = I2.mul(I2);        // I2^2
    Mat I1_2 = I1.mul(I1);        // I1^2
    Mat I1_I2 = I1.mul(I2);        // I1 * I2
    
    /*************************** END INITS **********************************/
    
    Mat mu1, mu2;   // PRELIMINARY COMPUTING
    GaussianBlur(I1, mu1, Size(11, 11), 1.5);
    GaussianBlur(I2, mu2, Size(11, 11), 1.5);
    
    Mat mu1_2 = mu1.mul(mu1);
    Mat mu2_2 = mu2.mul(mu2);
    Mat mu1_mu2 = mu1.mul(mu2);
    
    Mat sigma1_2, sigma2_2, sigma12;
    
    GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;
    
    GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;
    
    GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;
    
    ///////////////////////////////// FORMULA ////////////////////////////////
    Mat t1, t2, t3;
    
    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
    
    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
    
    Mat ssim_map;
    divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;
    
    Scalar mssim = mean(ssim_map); // mssim = average of ssim map
    return mssim;
}


// Convert double pointer to vector
vector<vector<unsigned char>> Arr2Vec(unsigned char** x, int jpegImageWidth, int jpegImageHeight)
{
    vector<vector<unsigned char>> temp;
    temp.resize(jpegImageHeight, vector<unsigned char>(jpegImageWidth, 0));
    
    for (int i = 0; i < jpegImageHeight; ++i) {
        for (int j = 0; j < jpegImageWidth; ++j) {
            temp[i][j] = x[i][j];
        }
    }
    
    return temp;
}

cv::Mat vec2mat(vector< vector<unsigned char> >& matrix) {
    
    // Copy data flag has to be true
    // [https://answers.opencv.org/question/96739/returning-a-mat-from-a-function/]
    cv::Mat matrix_CV(0, matrix.size(), cv::DataType<unsigned char>::type, true);
    
    for (unsigned int i = 0; i < matrix.size(); ++i)
    {
        // Make a temporary cv::Mat row and add to NewSamples _without_ data copy
        cv::Mat Sample(1, matrix[0].size(), cv::DataType<unsigned char>::type, matrix[i].data());
        matrix_CV.push_back(Sample);
    }
    return matrix_CV;
}

void reportMetrics(double x, double y, double z)
{
    cout << "X: " << x << "\nY: " << y << "\nZ: " << z << endl;
}


// parallel running
static bool runEncoderWithMultipleQF(std::string filename, std::string enc_path_to_files, std::string txt_path) {
    

    std::string main_dir = filename.substr(0, 52);
    std::string first_token = filename.substr(41);
    //    std::string main_dir = filename.substr(0,46);
    //    std::string first_token = filename.substr(45);
    //
    size_t found1 = first_token.find_last_of("/\\");
    std::string second_token = first_token.substr(0, found1 + 1);
    std::string third_token = first_token.substr(found1 + 1);
    size_t found2 = third_token.find_last_of(".");
    std::string fifth_token = third_token.substr(found2);
    std::string fourth_token = third_token.substr(0, found2);
    
    //    1 /shard-0/1/ILSVRC2012_val_00000001.JPEG
    //    2 /shard-0/1/
    //    3 ILSVRC2012_val_00000001.JPEG
    //    4 ILSVRC2012_val_00000001
    //    5 .JPEG
    //    6 ILSVRC2012_val_00000001
    //    std::cout <<"1 " + first_token << std::endl;
    //    std::cout <<"2 " + second_token << std::endl;
    //    std::cout <<"3 " + third_token << std::endl;
    //    std::cout <<"4 " + fourth_token << std::endl;
    //    std::cout <<"5 " + fifth_token << std::endl;
    //    std::cout <<"6 " + fourth_token << std::endl;
    
    
    std::string f1_yuv = filename;
    // jpeg_decoder Org_YUV(f1_yuv);
    // Org:
//    Mat image_org = imread(f1_yuv, CV_LOAD_IMAGE_COLOR);
    Mat image_org = imread(f1_yuv, IMREAD_UNCHANGED);
    Mat ycbcr_org;
    int nComponents = image_org.channels();
    int width_org =image_org.cols;
    int height_org =image_org.rows;
    cv::Mat ycbcr_channels[3];
    if (nComponents >1)
    {
        cv::cvtColor(image_org, ycbcr_org, cv::COLOR_BGR2YCrCb);
        cv::split(ycbcr_org, ycbcr_channels);
    }
    
    int nloop = 21;
//    cv::Mat ycbcr_channels_qf[21][3];
    // vector of image_qf
    std::vector<cv::Mat> image_qf(nloop);
    std::vector<std::vector<cv::Mat>> ycbcr_channels_qf;
    std::vector<cv::Mat> ycbcr_qf(nloop);
    std::vector<string> decoded_filename(nloop);
    
    ycbcr_channels_qf.resize(nloop, std::vector<cv::Mat> (nComponents));
    
    // Parallel version
    // number of threads from the given hardware
    const size_t nthreads = std::thread::hardware_concurrency();
    {
        // Pre loop
        //        std::cout <<"parallel (" << nthreads << " threads):" <<std::endl;
        std::vector<std::thread> threads(nthreads);
        std::mutex critical;
        
        // Create the threads
        for (int t = 0; t < nthreads; ++t)
        {
            threads[t] = std::thread(std::bind(
                                               [&](const int bi, const int ei, const int t)
                                               {
                                                   // loop over all items
                                                   for (int i = bi; i < ei; i++)
                                                   {
                                                       {


                                                           decoded_filename[i] = enc_path_to_files + second_token + fourth_token + "-QF-" + to_string(i*5) + fifth_token;
                                                           image_qf[i] = imread(decoded_filename[i], IMREAD_UNCHANGED);
                                                           if (nComponents >1)
                                                           {
                                                               cv::cvtColor(image_qf[i], ycbcr_qf[i], cv::COLOR_BGR2YCrCb);
                                                               cv::split(ycbcr_qf[i], ycbcr_channels_qf[i]);
                                                           }
                                                       }
                                                   }
                                                   
                                               },
                                               t * nloop / nthreads,
                                               (t + 1) == nthreads ? nloop : (t + 1) * nloop / nthreads,
                                               t)
                                     );
        }
        // Launch the threads:
        std::for_each(threads.begin(), threads.end(), [](std::thread& x) {x.join(); });
        
        
        
        // Calculate
        double psnr, psnr_Y, psnr_Cb, psnr_Cr;
        double mssim, mssim_Y, mssim_Cb, mssim_Cr;
        mssim_Y = mssim_Cb = mssim_Cr = 0;
        psnr_Y = psnr_Cb = psnr_Cr = 0;
        cv::Mat QF_mat_Cr_420, QF_mat_Cb_420, ORG_mat_Cr_420, ORG_mat_Cb_420;
        
//        for (int kDecoder = 0; kDecoder < nloop; ++kDecoder)
//        for (int kDecoder = 0; kDecoder < 1; ++kDecoder)
        for (int kDecoder = 0; kDecoder < nloop; ++kDecoder)
        {
            
            if (nComponents > 1)
            {
                cv::Size s = ycbcr_channels[0].size();
                int height = s.height;
                int width = s.width;
                
                
                psnr_Y = getPSNR(ycbcr_channels[0], ycbcr_channels_qf[kDecoder][0], width, height);
                Scalar final_mssim = getMSSIM(ycbcr_channels[0], ycbcr_channels_qf[kDecoder][0]);
                mssim_Y = final_mssim[0];
                // Org:
                // Cb component
                //         Calculate the quality metrics
                psnr_Cb = getPSNR(ycbcr_channels[1], ycbcr_channels_qf[kDecoder][1], width, height);
                psnr_Cr = getPSNR(ycbcr_channels[2], ycbcr_channels_qf[kDecoder][2], width, height);
                
                final_mssim = getMSSIM(ycbcr_channels[1], ycbcr_channels_qf[kDecoder][1]);
                mssim_Cb = final_mssim[0];
                final_mssim = getMSSIM(ycbcr_channels[2], ycbcr_channels_qf[kDecoder][2]);
                mssim_Cr = final_mssim[0];
                mssim = (6 * mssim_Y + mssim_Cb + mssim_Cr) / 8;
                psnr  = (6 * psnr_Y + psnr_Cb + psnr_Cr) / 8;
            }
            else
            {
                mssim = mssim_Y = getMSSIM(image_org, image_qf[kDecoder])[0];
                Scalar psnr_1   = getPSNR(image_org, image_qf[kDecoder], width_org, height_org);
                psnr = psnr_Y   = psnr_1[0];
            }
            
//             cout << kDecoder << "-" << psnr << ", " << mssim << endl;
//             cout << "qf: " << kDecoder << endl;
//             reportMetrics(psnr_Y, psnr_Cb, psnr_Cr);
//             reportMetrics(mssim_Y, mssim_Cb, mssim_Cr);
            
            
             // Write to out file:
             const int quality_factor = kDecoder * 5;
             ofstream myfile_txt;
             myfile_txt.open(txt_path + second_token + fourth_token + "-QF-" + to_string(quality_factor) + ".txt");
             myfile_txt << to_string(psnr) << "\t" << to_string(mssim) << "\n";
             myfile_txt.close();
            
            
            
        } // end kDecoder loop
        
        // Post loop:
        // ---------
        
        
        std::string encoded_filename = filename;
        
        ////// String Processing -- Get the file Name
        size_t found = encoded_filename.find_last_of("/\\");
        std::string filename_first_token = encoded_filename.substr(found + 1);
        found = filename_first_token.find_first_of(".");
        std::string filename_second_token = filename_first_token.substr(0, found);
        cout << "\n\nDone Execution; Output is: " << enc_path_to_files + second_token + fourth_token + "-QF-* .JPEG" << endl;
        cout << "\n\nDone Execution; Output is: " << txt_path + second_token + fourth_token + "-QF-* .txt" << endl;
    }
    
    return true;
    
}


#endif /* JPEG_ENCODER_MULTIPLE_QF_H */
