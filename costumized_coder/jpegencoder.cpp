//
//  jpegencoder.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-16
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//
#pragma warning (disable : 4996)
#include "jpegencoder.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "operation.h"
#include <math.h>
#include <cstring>
#include <filesystem>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

const double pi = 3.1415926535897932384626433832795;
using namespace cv;
namespace fs = std::experimental::filesystem;
using namespace std;

jpeg_encoder::jpeg_encoder(const std::string filename) : quantization_table_write_process_luminance(true, true), quantization_table_write_process_chrominance(true, false) {
    
    image_to_export_filename = filename;
	counter_bits_buffer = 0;
    // Set the number of bytes in the buffer to 0
    num_bytes_in_jpeg_enc_write_buffer = 0;

	//const float PI = 3.14f;
	//const double inv16 = 1.0 / 16.0;
	//for (int i = 0; i < 8; i++)
	//{
	//	for (int j = 0; j < 8; j++)
	//	{
	//		cosine_idct[j][i] = cosf((2.0 * i + 1) * j * PI * inv16);
	//	}
	//}
	cosine_idct.resize(8, vector<double>(8, 0));
	cosine_idct[0][0] = 1.000000;
	cosine_idct[1][0] = 0.980785;
	cosine_idct[2][0] = 0.923880;
	cosine_idct[3][0] = 0.831470;
	cosine_idct[4][0] = 0.707107;
	cosine_idct[5][0] = 0.555570;
	cosine_idct[6][0] = 0.382683;
	cosine_idct[7][0] = 0.195090;
	cosine_idct[0][1] = 1.000000;
	cosine_idct[1][1] = 0.831470;
	cosine_idct[2][1] = 0.382683;
	cosine_idct[3][1] = -0.195090;
	cosine_idct[4][1] = -0.707107;
	cosine_idct[5][1] = -0.980785;
	cosine_idct[6][1] = -0.923880;
	cosine_idct[7][1] = -0.555570;
	cosine_idct[0][2] = 1.000000;
	cosine_idct[1][2] = 0.555570;
	cosine_idct[2][2] = -0.382684;
	cosine_idct[3][2] = -0.980785;
	cosine_idct[4][2] = -0.707107;
	cosine_idct[5][2] = 0.195090;
	cosine_idct[6][2] = 0.923880;
	cosine_idct[7][2] = 0.831469;
	cosine_idct[0][3] = 1.000000;
	cosine_idct[1][3] = 0.195090;
	cosine_idct[2][3] = -0.923880;
	cosine_idct[3][3] = -0.555570;
	cosine_idct[4][3] = 0.707107;
	cosine_idct[5][3] = 0.831469;
	cosine_idct[6][3] = -0.382684;
	cosine_idct[7][3] = -0.980785;
	cosine_idct[0][4] = 1.000000;
	cosine_idct[1][4] = -0.195090;
	cosine_idct[2][4] = -0.923880;
	cosine_idct[3][4] = 0.555570;
	cosine_idct[4][4] = 0.707107;
	cosine_idct[5][4] = -0.831470;
	cosine_idct[6][4] = -0.382683;
	cosine_idct[7][4] = 0.980785;
	cosine_idct[0][5] = 1.000000;
	cosine_idct[1][5] = -0.555570;
	cosine_idct[2][5] = -0.382683;
	cosine_idct[3][5] = 0.980785;
	cosine_idct[4][5] = -0.707107;
	cosine_idct[5][5] = -0.195090;
	cosine_idct[6][5] = 0.923879;
	cosine_idct[7][5] = -0.831470;
	cosine_idct[0][6] = 1.000000;
	cosine_idct[1][6] = -0.831470;
	cosine_idct[2][6] = 0.382684;
	cosine_idct[3][6] = 0.195090;
	cosine_idct[4][6] = -0.707107;
	cosine_idct[5][6] = 0.980785;
	cosine_idct[6][6] = -0.923880;
	cosine_idct[7][6] = 0.555571;
	cosine_idct[0][7] = 1.000000;
	cosine_idct[1][7] = -0.980785;
	cosine_idct[2][7] = 0.923880;
	cosine_idct[3][7] = -0.831470;
	cosine_idct[4][7] = 0.707107;
	cosine_idct[5][7] = -0.555571;
	cosine_idct[6][7] = 0.382684;
	cosine_idct[7][7] = -0.195092;



}

// this is what we work with 
jpeg_encoder::jpeg_encoder(jpeg_decoder* processed_jpeg_decoder, std::string output_filename) : quantization_table_write_process_luminance(true, true), quantization_table_write_process_chrominance(true, false) {
    counter_bits_buffer = 0;
    // set the output file name
    image_to_export_filename = output_filename;
    
    // TODO: that's indeed not safe -- you need to create copy constructor or make it unique_ptr
    jpegDecoder = processed_jpeg_decoder;
    
    // TODO: image width and height must be multiples of 8. If they are not, then we have to set them up
    
    
    // set image with and height:
    image_to_export_width = processed_jpeg_decoder->get_image_width();
    image_to_export_height = processed_jpeg_decoder->get_image_height();
    
    // Set the upscaled image width and height parameters
    image_to_export_width_dct = processed_jpeg_decoder->upscale_width;
    image_to_export_height_dct = processed_jpeg_decoder->upscale_height;
    
    
    // boolean that test whether we padded the input picture
    if (image_to_export_width_dct != image_to_export_width)  padded_width = true;
    if (image_to_export_height != image_to_export_height_dct)  padded_height = true;
    
    // TODO: remove thiss
    total_block = image_to_export_width_dct * image_to_export_height_dct / 64;
    total_block_C = ceil(sqrt(total_block) / 2.0) * ceil(sqrt(total_block) / 2.0);
    
    // Set the counts to 0
    count_block_Y = count_block_Cb = count_block_Cr = 0;
    
    // Set the counts to 0
    count_block_Y_huffman =  0;
    count_block_Cb_huffman = 0;
    count_block_Cr_huffman = 0;
    
    // queue reading stuff
    g_nbits_in_reservoir = 0;
    // the queue of bits
    m_bit_buffer = 0;
    // queue reading stuff
    m_bits_in = 0;
    
    // Set the number of bytes in the buffer to 0
    num_bytes_in_jpeg_enc_write_buffer = 0;
    
    // Set the quality factor to the default value (-1)
    quality_factor = -1;
    
}


jpeg_encoder::jpeg_encoder(jpeg_decoder* processed_jpeg_decoder, std::string output_filename, int input_qf) : quantization_table_write_process_luminance(true, true), quantization_table_write_process_chrominance(true, false) {
    
    counter_bits_buffer = 0;
    // set the output file name
    image_to_export_filename = output_filename;
    
    // TODO: that's indeed not safe -- you need to create copy constructor or make it unique_ptr
    jpegDecoder = processed_jpeg_decoder;
    
    // TODO: image width and height must be multiples of 8. If they are not, then we have to set them up
    
    
    // set image with and height:
    image_to_export_width = processed_jpeg_decoder->get_image_width();
    image_to_export_height = processed_jpeg_decoder->get_image_height();
    
    // Set the upscaled image width and height parameters
    image_to_export_width_dct = processed_jpeg_decoder->upscale_width;
    image_to_export_height_dct = processed_jpeg_decoder->upscale_height;
    
    
    // boolean that test whether we padded the input picture
    if (image_to_export_width_dct != image_to_export_width)  padded_width = true;
    if (image_to_export_height != image_to_export_height_dct)  padded_height = true;
    
    // TODO: remove thiss
    total_block = image_to_export_width_dct * image_to_export_height_dct / 64;
    total_block_C = ceil(sqrt(total_block) / 2.0) * ceil(sqrt(total_block) / 2.0);
    
    // Set the counts to 0
    count_block_Y = count_block_Cb = count_block_Cr = 0;
    // queue reading stuff
    g_nbits_in_reservoir = 0;
    // the queue of bits
    m_bit_buffer = 0;
    // queue reading stuff
    m_bits_in = 0;
    
    // Set the number of bytes in the buffer to 0
    num_bytes_in_jpeg_enc_write_buffer = 0;
    
    // Set the quality factor
    quality_factor = input_qf;
    
}



void jpeg_encoder::copy_qTables() {
    
    int comp_size = static_cast<int>(jpegDecoder->components.size());
    
    //    // Copy the quantization tables
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            
#if IS_ENABLE_USE_QTABLE_FROM_PICTURE
            QuantizationTable * from = jpegDecoder->components[COMPONENT_Y].componentQuantizationTable;
            quantization_table_write_process_luminance.quantizationTableData[i][j] = from->quantizationTableData[i][j];
#endif
            if (quality_factor >= 0) {
                
                // The range should start from 1; quality factor cannot be bigger than 0xFFFF
                if(quality_factor == 0) quality_factor = 1;
                
                int S;
                if (quality_factor < 50)
                    S = floor(5000 / quality_factor);
                else
                    S = 200 - 2 * quality_factor;
                quantization_table_write_process_luminance.quantizationTableData[i][j] = floor((S*quantization_table_write_process_luminance.quantizationTableData[i][j] + 50) / 100);
                quantization_table_write_process_luminance.quantizationTableData[i][j] = std::min(quantization_table_write_process_luminance.quantizationTableData[i][j], (1 << jpegDecoder->jpegImageSamplePrecision) - 1);
                if (quantization_table_write_process_luminance.quantizationTableData[i][j] == 0) quantization_table_write_process_luminance.quantizationTableData[i][j] = 1;
                
            }
            
            
            // TODO: check if you only the QTable for Cb and not Cr
            if (comp_size > 1) {
                
#if IS_ENABLE_USE_QTABLE_FROM_PICTURE
                from = jpegDecoder->components[COMPONENT_Cb].componentQuantizationTable;
                quantization_table_write_process_chrominance.quantizationTableData[i][j] = from->quantizationTableData[i][j];
#endif
                if (quality_factor >= 1) {
                    int S;
                    if (quality_factor < 50)
                        S = floor(5000 / quality_factor);
                    else
                        S = 200 - 2 * quality_factor;
                    
                    quantization_table_write_process_chrominance.quantizationTableData[i][j] = floor((S*quantization_table_write_process_chrominance.quantizationTableData[i][j] + 50) / 100);
                    quantization_table_write_process_chrominance.quantizationTableData[i][j] = std::min(quantization_table_write_process_chrominance.quantizationTableData[i][j], (1 << jpegDecoder->jpegImageSamplePrecision) - 1);
                    
                    if (quantization_table_write_process_chrominance.quantizationTableData[i][j] == 0) quantization_table_write_process_chrominance.quantizationTableData[i][j] = 1;
                }
            }
            
        } // end inner
    } // end outer
    
} // end copy_qTables

// NEW to TCM: Apply TCM
void::jpeg_encoder::perform_TCM() {
    int peak; double prob, lambda;
    double yc;
    // Create the yc array and put them all zeros (Note: ignore DC)
    yc_array.resize(64, 0);
    count_outlier_list.resize(total_block, 0);
    
    for (int i = 1; i < 64; ++i) {
        // function of TCM
        tcm::TCMprocessOneSequence(jpegDecoder->tCoeff_Y_AC.at(i), total_block, &peak, &prob, &lambda, &yc);
        
        // store output yc:
        yc_array.at(i) = yc;
    }
    
    
    
    // TODO: remove this: visual output of yc score
#if DEBUGLEVEL > 20
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
            cout << "   " << yc_array.at(i + j * 8);
        
        cout << "\n";
    }
#endif
    
    tcm::count_outliers(yc_array, total_block, jpegDecoder->tCoeff_Y_AC, count_outlier_list);
    
    //for (int i = 0; i < 8; i++) {
    //    for (int j = 0; j < 8; j++) {
    //        cout << jpegDecoder->tCoeff_Y_AC[i * 8 + j][3200] << " ";
    //        if (j == 7) cout << endl;
    //    }
    //}
    
    // counter for wiped blocks
    int wiped_counter = 0;
    int counter_wiped_AC = 0;
    
    count_outlier_list_sort = count_outlier_list;
    sort(count_outlier_list_sort.begin(), count_outlier_list_sort.end());
    int min = count_outlier_list.at(0);
    int max = -1;
    for (int i = 0; i < count_outlier_list.size(); ++i) {
        if (count_outlier_list.at(i) < min) {
            min = count_outlier_list.at(i);
        }
        
        if (max < count_outlier_list.at(i)) {
            max = count_outlier_list.at(i);
        }
    }
    cout << "Max OBF: " << max << "; Min OBF: " << min<<endl;
    
    
    ////// String Processing -- Get the file Name
    //    std::string encoded_filename = image_to_export_filename;
    //    size_t found = encoded_filename.find_last_of("/\\");
    //    std::string filename_first_token = encoded_filename.substr(found+1);
    //    found = filename_first_token.find_first_of(".");
    //    std::string filename_second_token = filename_first_token.substr(0, found);
    //    string tcm_store = "";
    //    std::ostringstream oss;
    //    oss << "/Volumes/DATA/ml/tcm_out6/" << filename_second_token << "_" << max;
    //    tcm_store = oss.str();
    //    char* pYUVFileName = tcm_store.empty()? NULL: strdup(tcm_store.c_str());
    //    FILE* sastre_pFile = fopen (pYUVFileName, "w");
    //    fprintf(sastre_pFile, "%d\n", max);
    //    fclose(sastre_pFile);
    
    // cout << count_outlier_list_sort[floor(count_outlier_list_sort.size()*0.70)];
    
    // Beform Doing TCM Record the histrogram of tCoeff_Y for reference
    
    // _______________________________APPLICATION OF TCM _____________________________________
    int currentComponent = COMPONENT_Y;
    
    for (currentComponent = COMPONENT_Y; currentComponent < jpegDecoder->components.size(); ++currentComponent) {
        
        // How many Y's in the horizontal and vertical direction (2x2 is the usual case)
        int HFactor = jpegDecoder->components[currentComponent].HFactor, VFactor = jpegDecoder->components[currentComponent].VFactor;
        
        // These Y's should be scaled/repeated how many times horizonatally and vertically
        int HScale = jpegDecoder->components[currentComponent].HScale, VScale = jpegDecoder->components[currentComponent].VScale;
        
        // Note: use the upscaled width and height to store the DCT coefficientss
        int comp_height = ceil(1.0* image_to_export_height_dct / VScale);
        int comp_width = ceil(1.0* image_to_export_width_dct / HScale);
        int width = comp_width, height = comp_height;
        
        // Initialize the positions
        int currentX = 0;
        int currentY = 0;
        int currentBlockHFactor = 0;
        int currentBlockVFactor = 0;
        int smooth_value = 0;
        
        // loop on the entire picture
        for (uint i = 0; i < comp_height; i += 8) {
            for (uint j = 0; j < comp_width; j += 8) {
                // fetch the block
                int block_idx = (currentX / 8) + (currentY / 8) * (width / 8);
                
                //wiped counter
                if (count_outlier_list.at(block_idx) <= count_outlier_list_sort[floor(count_outlier_list_sort.size()*0.70)]) { //TCM_OUTLIER_THRESHOLD
                    wiped_counter++;
                }
                
                //if (block_idx >=0) {
                //    cout << block_idx << endl;
                //    //cout << &jpegDecoder->m_YPicture_buffer[383][511];
                //}
                //smooth_value = smoothing(block_idx, count_outlier_list, height, width, currentComponent, jpegDecoder->m_YPicture_buffer, jpegDecoder->m_CbPicture_buffer, jpegDecoder->m_CrPicture_buffer);
                uint y;
                for (y = 0; y < 8; ++y) {
                    bool done = false;
                    
                    if (currentX >= width) break;
                    if (currentY + y >= height) break;
                    
                    int picture_y = currentY + y;
                    
                    uint x;
                    for (x = 0; x < 8; ++x) {
                        
                        if (currentX + x >= width) break;
                        int realx = currentX + x;
                        
#if IS_ONLY_TCM
                        if (count_outlier_list.at(block_idx) == 0) { //<= TCM_OUTLIER_THRESHOLD
                            //cout << &jpegDecoder->m_YPicture_buffer[383][511];
                            if (!(x == 0 && y == 0)) {
                                switch (currentComponent) {
                                    case COMPONENT_Y:
                                        if (jpegDecoder->tCoeff_Y[picture_y][realx] != 0) counter_wiped_AC++;
                                        jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                                        break;
                                    case COMPONENT_Cb:
                                        jpegDecoder->tCoeff_Cb[picture_y][realx] = 0;
                                        break;
                                    case COMPONENT_Cr:
                                        jpegDecoder->tCoeff_Cr[picture_y][realx] = 0;
                                        break;
                                    default:
                                        jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                                        break;
                                        
                                } // end switch
                            }
                        }
#else
                        if (smooth_value == 0) {
                            //    counter++;
                            //    cout << counter << endl;
                            //    //jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                            perform_content_adaptive_dct(currentX, currentY, currentComponent, jpegDecoder->tCoeff_Y, jpegDecoder->tCoeff_Cb, jpegDecoder->tCoeff_Cr);
                            done = true;
                            break;
                        }
                        
#endif
                        
                        //if (count_outlier_list.at(block_idx) <= 1) {
                        //#if IS_ONLY_TCM
                        //                            jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                        //#endif
                        //perform_content_adaptive_dct(currentX, currentY, currentComponent, jpegDecoder->tCoeff_Y, jpegDecoder->tCoeff_Cb, jpegDecoder->tCoeff_Cr);
                        //                            if (smooth_value != 0) {
                        //                                if (x == 0 && y == 0) {
                        //                                    //jpegDecoder->tCoeff_Y[picture_y][realx] = smooth_value;
                        //                                }
                        //                                //jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                        //                            }
                        //                            else {
                        //                                //jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                        //                                perform_content_adaptive_dct(currentX, currentY, currentComponent, jpegDecoder->tCoeff_Y, jpegDecoder->tCoeff_Cb, jpegDecoder->tCoeff_Cr);
                        //                                done = true;
                        //                                break;
                        //                            }
                        //
                        //}
                        //                        else {
                        //                            if (smooth_value == 0) {
                        //                                //jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
                        //                                perform_content_adaptive_dct(currentX, currentY, currentComponent, jpegDecoder->tCoeff_Y, jpegDecoder->tCoeff_Cb, jpegDecoder->tCoeff_Cr);
                        //                                done = true;
                        //                                break;
                        //                        }
                        
                        //if (currentComponent == 0 && block_idx == 3200) { //|| (block_idx > 50 && block_idx <=94)
                        //    cout << "position X: " << currentX << "  position Y:" << currentY << endl;
                        //    jpegDecoder->tCoeff_Y[picture_y][realx] = -128;
                        //}
                    }
                    if (done) break;
                }
                
                // Adjust the indices:
                HScale = 1; // note: you will store the dct coefficients without any scaling
                VScale = 1; // note: algorithm below works, but you have to set the scales into 1x1
                currentX += 8 * HScale;
                currentBlockHFactor++;
                // you made a line of blocks
                if (currentBlockHFactor >= HFactor) {
                    
                    // restore the current X to its initial position and reset the counters
                    currentX -= 8 * HScale * HFactor;
                    currentBlockHFactor = 0;
                    
                    // go to next line
                    currentY += 8 * VScale;
                    currentBlockVFactor++;
                    
                    // you made a column of blocks
                    if (currentBlockVFactor >= VFactor) {
                        
                        // restore the current Y to its initial position and reset the counters
                        currentY -= 8 * VScale * VFactor;
                        currentBlockVFactor = 0;
                        currentX += 8 * HScale * HFactor;
                        if (currentX >= width) {
                            currentX = 0;
                            currentY += 8 * VScale * VFactor;
                        } // end if (currentBlockVFactor_dct[comp] >= VFactor)
                        
                    } // end  if (currentBlockVFactor_dct[comp] >= VFactor)
                } // end if(currentBlockHFactor_dct[comp] >= HFactor)
            }
        }
        
        /*for (int i = 0; i < jpegDecoder->tCoeff_Y.size(); ++i) {
         for (int j = 0; j < jpegDecoder->tCoeff_Y[i].size(); ++j) {
         cout << jpegDecoder->tCoeff_Y[i][j] << " ";
         }
         }*/
        
        //________________________________________________________________________________________
        
        // _______________________________APPLICATION OF SMOOTHING _______________________________
        
        // Initialize the positions
        //currentX = 0;
        //currentY = 0;
        //currentBlockHFactor = 0;
        //currentBlockVFactor = 0;
        
        //// loop on the entire picture
        //for (uint i = 0; i < comp_height; i += 8) {
        //    for (uint j = 0; j < comp_width; j += 8) {
        //        // fetch the block
        //        int block_idx = (currentX / 8) + (currentY / 8) * (width / 8);
        //        uint y;
        //        for (y = 0; y < 8; ++y) {
        
        //            if (currentX >= width) break;
        //            if (currentY + y >= height) break;
        
        //            int picture_y = currentY + y;
        
        //            uint x;
        //            for (x = 0; x < 8; ++x) {
        
        //                if (currentX + x >= width) break;
        //                int realx = currentX + x;
        
        //                if (count_outlier_list.at(block_idx) <= 1) {
        //                    jpegDecoder->tCoeff_Y[picture_y][realx] = 0;
        //                }
        //            }
        //        }
        
        //        // Adjust the indices:
        //        HScale = 1; // note: you will store the dct coefficients without any scaling
        //        VScale = 1; // note: algorithm below works, but you have to set the scales into 1x1
        //        currentX += 8 * HScale;
        //        currentBlockHFactor++;
        //        // you made a line of blocks
        //        if (currentBlockHFactor >= HFactor) {
        
        //            // restore the current X to its initial position and reset the counters
        //            currentX -= 8 * HScale * HFactor;
        //            currentBlockHFactor = 0;
        
        //            // go to next line
        //            currentY += 8 * VScale;
        //            currentBlockVFactor++;
        
        //            // you made a column of blocks
        //            if (currentBlockVFactor >= VFactor) {
        
        //                // restore the current Y to its initial position and reset the counters
        //                currentY -= 8 * VScale * VFactor;
        //                currentBlockVFactor = 0;
        //                currentX += 8 * HScale * HFactor;
        //                if (currentX >= width) {
        //                    currentX = 0;
        //                    currentY += 8 * VScale * VFactor;
        //                } // end if (currentBlockVFactor_dct[comp] >= VFactor)
        
        //            } // end  if (currentBlockVFactor_dct[comp] >= VFactor)
        //        } // end if(currentBlockHFactor_dct[comp] >= HFactor)
        //    }
        //}
        
        //__________________________________________________________________________
        
        
        
#if DEBUGLEVEL > 20
        cout << endl << "--------------------------" << endl;
        for (int i = 0; i < 64; ++i) {
            for (int j = 0; j < 64; ++j)
                cout << "   " << count_outlier_list.at(i + j * 64);
            
            cout << "\n";
        }
#endif
        
        
        //cout << "\nNumber of wiped blocks " << wiped_counter << " --- Percentage " << 100.0*(1.0*wiped_counter / total_block) << endl;
        // cout << "Percentage of Wiped Pixel " << 100.0*counter_wiped_AC/(image_to_export_height*image_to_export_width) << "% " << endl;
        // getchar();
    } // end loop on currentComponent
    
    
} // end perform_tcm

void::jpeg_encoder::cal_stat(vector<int> &zz, vector<int> & counts, int totalpair, vector<int>& dccounts, int count_block)
{
    int temp, nbits;
    int k, r;
    
    /* encode DC coefficient */
    int start_index = returnIndexInZigZagArray(count_block);
    const int dc_delta = zz.at(start_index);
//    temp=zz[0];
    temp = dc_delta;
    if (temp < 0)
        temp = -temp;   /* temp is the absulate value of input coefficient */
    
    /* find the number of bits needed for the magnitude of the coeffcient */
    nbits = 0;
    while(temp)
    {
        nbits++;
        temp >>= 1;
    }
    dccounts[nbits]++;
    
    /* encode AC coefficients (refer to Figure F.1.2.2 */
    r = 0;
//    for (k = 1 + start_index; k < 64 + start_index; k++)
    // Modified this one to start from the right index
    for (k = 1 + start_index; k < 64 + start_index; ++k)
    {
        if ((temp=zz[k])==0)
            r++;
        else
        {
            /* if run length is greater than 15, emit special code 0xF0 */
            while (r > 15)
            {
                counts[15<<4]++;
                totalpair = totalpair+1;
                r = r-16;
            }
            if (temp < 0)
                temp = -temp;
            
            /* find the number of bits needed for the magnitude of the coefficient */
            nbits = 1; /* Nonzero AC value  is at least 1 bits long */
            while (temp >>= 1)
                nbits++;
            
            /* Count Huffman symbol for (run, size) */
            counts[(r<<4) + nbits]++;
            totalpair = totalpair+1;
            r = 0;
        }
    }
    /* if the last coefficients are zero, emit an end-of-block code */
    if (r>0)
    {
        counts[0]++;
        totalpair = totalpair+1;
    }
}

//// zigZagArray variable will hold data after zigZag transform
////// image variable holds data for one of three components (YCbCr)
void::jpeg_encoder::perform_fdct(uint_8 ** image, vector<int> &zigZagArray, int quantizationTable[8][8], int currentComponent) {
    
    vector<vector<double> >block8x8(8, vector<double>(8));
    double previousDCCoefficient = 0; //In this function DCPM is performed so DC coeficient is generated with formula Diff=DCi-DCi-1
    double prev_previousDC = 0;// previous previous
    
    // How many Y's in the horizontal and vertical direction (2x2 is the usual case)
    int HFactor = jpegDecoder->components[currentComponent].HFactor, VFactor = jpegDecoder->components[currentComponent].VFactor;
    
    // These Y's should be scaled/repeated how many times horizonatally and vertically
    int HScale = jpegDecoder->components[currentComponent].HScale, VScale = jpegDecoder->components[currentComponent].VScale;
    
    // Note: use the upscaled width and height to store the DCT coefficientss
    int comp_height = ceil(1.0* image_to_export_height_dct / VScale);
    int comp_width = ceil(1.0* image_to_export_width_dct / HScale);
    int width = comp_width, height = comp_height;
    
    // Initialize the positions
    int currentX = 0;
    int currentY = 0;
    int currentBlockHFactor = 0;
    int currentBlockVFactor = 0;
    
    // Initialize the counts
    count_block_Y_huffman = 0;
    count_block_Cb_huffman = 0;
    count_block_Cr_huffman = 0;
    
#if DEBUGLEVEL > 20
    cout << "Now airing from the encoder ;) " << endl;
    cout << "Component: " << currentComponent << ", comp_with: " << comp_width <<
    ", comp_height: " << comp_height << endl;
#endif
    
    
    // loop on the entire picture
    for (uint i = 0; i < comp_height; i += 8) {
        for (uint j = 0; j < comp_width; j += 8) {
            // fetch the block
#if DEBUGLEVEL > 20
            cout << "Top Left: " << currentX << ", " << currentY << ", Bottom Right: " << currentX + 8 << ", " << currentY + 8 << endl;
#endif
            uint y;
            for (y = 0; y < 8; ++y) {
                
                if (currentX >= width) break;
                if (currentY + y >= height) break;
                
                int picture_y = currentY + y;
                
                uint x;
                for (x = 0; x < 8; ++x) {
                    
                    if (currentX + x >= width) break;
                    int realx = currentX + x;
                    
                    if (currentComponent == COMPONENT_Y)
                        block8x8[y][x] = jpegDecoder->tCoeff_Y[picture_y][realx]; // these are dct coeeficients
                    else if (currentComponent == COMPONENT_Cb)
                        block8x8[y][x] = jpegDecoder->tCoeff_Cb[picture_y][realx];
                    else
                        block8x8[y][x] = jpegDecoder->tCoeff_Cr[picture_y][realx];
                }
            }
            
            // DCPM for the DC element of the compononet (including quantization)
            for (int u = 0; u < 8; ++u) {
                for (int v = 0; v < 8; ++v) {
                    block8x8[u][v] = round(block8x8[u][v] / quantizationTable[u][v]);
                    
                    if (!u && !v) {
                        
                        // Store the DC coefficient
                        prev_previousDC = block8x8[u][v];
                        
#if DEBUGLEVEL > 20
                        cout << "Current DC coefficient: " << block8x8[u][v] << ", Prev: " << previousDCCoefficient << endl;
#endif
                        
                        block8x8[u][v] = block8x8[u][v] - previousDCCoefficient;
                        
                        // quantization
                        previousDCCoefficient = prev_previousDC; // set the previous DC
                    }
                }
            }
            // here I can just do inverse dct and dequantization and save the image. 






            // Perform ZigZag coding, and array is ready for Huffman coding
            ZigZagCoding(block8x8, zigZagArray);
            
            // Perform cal stats after zigzag, quantization, and dct
            if(currentComponent == COMPONENT_Y)
            {
                cal_stat(zigZagArray, ac_counts_Y, total, dc_counts_Y, count_block_Y_huffman);
                count_block_Y_huffman++;
            }
            else if(currentComponent == COMPONENT_Cb)
            {
                cal_stat(zigZagArray, ac_counts_CbCr, total, dc_counts_CbCr, count_block_Cb_huffman);
                count_block_Cb_huffman++;
            }
            else if(currentComponent == COMPONENT_Cr)
            {
                cal_stat(zigZagArray, ac_counts_CbCr, total, dc_counts_CbCr, count_block_Cr_huffman);
                count_block_Cr_huffman++;
            }
            
            
            // TODO: remove FANKOOSH
#if DEBUGLEVEL > 50
            static int fankoosh = 0;
            fankoosh++;
            if (currentComponent == COMPONENT_Y) {
                ofstream myfile;
                std::string path_to_files = "C:/Users/y77jiang/OneDrive - University of Waterloo/5e. TCM-Inception C++/jpeg_tcm/dataset/";
                std::string output_csv_name = path_to_files + "Lena_Y_enc.csv";
                myfile.open(output_csv_name, std::ofstream::out | std::ofstream::app);
                
                std::stringstream oss;
                std::size_t found = output_csv_name.find_last_of(".");
                std::string path_with_name = output_csv_name.substr(0, found);
                found = output_csv_name.find_last_of("/\\");
                std::string name_file_only = path_with_name.substr(found + 1);
                
                myfile << currentX << "-" << currentY << "\n";
                int k = 0;
                for (int i = 0; i < 8; ++i) {
                    for (int j = 0; j < 8; ++j) {
                        //                        myfile << block8x8[i][j] << ",";
                        int val = zigZagArray.at(zigZagArray.size() - 64 + k++);
                        myfile << val << ",";
                        //                        myfile << block8x8[i][j] << ",";
                    }
                    
                    //                    if( (i + 1) < 8){
                    myfile << "\n";
                    //                    }
                }
                myfile.close();
            }
#endif
            
            // Adjust the indices:
            HScale = 1; // note: you will store the dct coefficients without any scaling
            VScale = 1; // note: algorithm below works, but you have to set the scales into 1x1
            currentX += 8 * HScale;
            currentBlockHFactor++;
            // you made a line of blocks
            if (currentBlockHFactor >= HFactor) {
                
                // restore the current X to its initial position and reset the counters
                currentX -= 8 * HScale * HFactor;
                currentBlockHFactor = 0;
                
                // go to next line
                currentY += 8 * VScale;
                currentBlockVFactor++;
                
                // you made a column of blocks
                if (currentBlockVFactor >= VFactor) {
                    
                    // restore the current Y to its initial position and reset the counters
                    currentY -= 8 * VScale * VFactor;
                    currentBlockVFactor = 0;
                    currentX += 8 * HScale * HFactor;
                    if (currentX >= width) {
                        currentX = 0;
                        currentY += 8 * VScale * VFactor;
                    } // end if (currentBlockVFactor_dct[comp] >= VFactor)
                    
                } // end  if (currentBlockVFactor_dct[comp] >= VFactor)
            } // end if(currentBlockHFactor_dct[comp] >= HFactor)
        }
    }
    
} // end perform_fdct



// DCT2 - Discrete Cosine Transform - vector output
void jpeg_encoder::DCT2(int ** input, vector<vector<double>> &output)
{
    cout << "Test in DCT..." << endl << endl;
    double ALPHA, BETA;
    int row = 8;
    int col = 8;
    int u = 0;
    int v = 0;
    int i = 0;
    int j = 0;
    
    for (u = 0; u < row; u++)
    {
        for (v = 0; v < col; v++)
        {
            if (u == 0)
            {
                ALPHA = sqrt(1.0 / row);
            }
            else {
                ALPHA = sqrt(2.0 / row);
            }
            
            if (v == 0)
            {
                BETA = sqrt(1.0 / col);
            }
            else {
                BETA = sqrt(2.0 / col);
            }
            
            double tmp = 0.0;
            for (i = 0; i < row; i++) {
                for (j = 0; j < col; j++) {
                    tmp += *((int*)input + col * i + j) * cos((2 * i + 1)*u*pi / (2.0 * row)) * cos((2 * j + 1)*v*pi / (2.0 * col));
                }
            }
            output[u][v] = ALPHA * BETA * tmp;
        }
    }
    
    /*cout << "the result of dct:" << endl;
     for (int m = 0; m < row; m++) {
     for (int n = 0; n < col; n++) {
     cout << setw(8) << output[m][n] << " \t";
     }
     cout << endl;
     }*/
}

float jpeg_encoder::C_dct(int u) {
    
    if (u == 0)
    {
        return 1.0f / sqrt(8.0f);
    }
    else
    {
        return sqrt(2.0 / 8.0);
    }
}


double jpeg_encoder::func_dct(int x, int y, const int block[8][8]) {
    const float PI = 3.14f;
    float sum = 0;
    for (int u = 0; u < 8; ++u)
    {
        for (int v = 0; v < 8; ++v)
        {
            sum += (block[u][v] * cosf(((2 * v + 1) * x * PI) / 16)  * cosf(((2 * u + 1) * y * PI) / 16));
        } // end inner loop
    } // end outer loop
    
    sum = C_dct(x) * C_dct(y) * sum;
    return sum;
    
}

void jpeg_encoder::perform_dct(vector<vector<double> > &outBlock, int inBlock[8][8]) {
    
    for (int y = 0; y < 8; ++y)
    {
        for (int x = 0; x < 8; ++x)
        {
            outBlock[y][x] = func_dct(x, y, inBlock);
            
        } // end inner loop
    } // end outer loop
    
}


char jpeg_encoder::getCategoryOfDCTCoefficient(int x) {
    
    if (x == 0)
        return 0;
    else if (x == -1 || x == 1)
        return 1;
    else if ((x >= -3 && x <= -2) || (x >= 2 && x <= 3))
        return 2;
    else if ((x >= -7 && x <= -4) || (x >= 4 && x <= 7))
        return 3;
    else if ((x >= -15 && x <= -8) || (x >= 8 && x <= 15))
        return 4;
    else if ((x >= -31 && x <= -16) || (x >= 16 && x <= 31))
        return 5;
    else if ((x >= -63 && x <= -32) || (x >= 32 && x <= 63))
        return 6;
    else if ((x >= -127 && x <= -64) || (x >= 64 && x <= 127))
        return 7;
    else if ((x >= -255 && x <= -128) || (x >= 128 && x <= 255))
        return 8;
    else if ((x >= -511 && x <= -256) || (x >= 256 && x <= 511))
        return 9;
    else if ((x >= -1023 && x <= -512) || (x >= 512 && x <= 1023))
        return 10;
    else if ((x >= -2047 && x <= -1024) || (x >= 1024 && x <= 2047))
        return 11;
    else if ((x >= -4095 && x <= -2048) || (x >= 2048 && x <= 4095))
        return 12;
    else if ((x >= -8191 && x <= -4096) || (x >= 4096 && x <= 8191))
        return 13;
    else if ((x >= -16383 && x <= -8192) || (x >= 8192 && x <= 16383))
        return 14;
    else if ((x >= -32767 && x <= -16384) || (x >= 16384 && x <= 32767))
        return 15;
    else
        return 0;
    
} // end getCategoryOfDCTCoefficient


void jpeg_encoder::ZigZagCoding(int block8x8[8][8], vector<char>&zigZagArray) {
    //k- is zigZagArray index, i,j are index of matrix
    zigZagArray.push_back(block8x8[0][0]);//Take the first element
    int i = 0, j = 1;//Define index for matrix
    while (1) {
        while (j != 0 && i != 7) {//Going upside down until j!=0
            zigZagArray.push_back(block8x8[i][j]);
            i = i + 1;
            j = j - 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take the edge element
        
        if (i<7)//If not last row, increment i
            i = i + 1;
        
        else if (i == 7)//If we hit the last row, we go right one place
            j = j + 1;
        
        
        while (i != 0 && j != 7) {//Going bottom up
            zigZagArray.push_back(block8x8[i][j]);
            i = i - 1;
            j = j + 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take edge element
        if (j<7)//If we didn't hit the edge, increment j
            j = j + 1;
        
        else if (j == 7)//If we hit the last element, go down one place
            i = i + 1;
        
        if (i >= 7 && j >= 7)//If we hit last element matrix[8][8] exit
            break;
    }
} // end ZigZagCoding


void jpeg_encoder::ZigZagCoding(vector<vector<int> > &block8x8, vector<char>&zigZagArray) {
    //k- is zigZagArray index, i,j are index of matrix
    zigZagArray.push_back(block8x8[0][0]);//Take the first element
    int i = 0, j = 1;//Define index for matrix
    while (1) {
        while (j != 0 && i != 7) {//Going upside down until j!=0
            zigZagArray.push_back(block8x8[i][j]);
            i = i + 1;
            j = j - 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take the edge element
        
        if (i<7)//If not last row, increment i
            i = i + 1;
        
        else if (i == 7)//If we hit the last row, we go right one place
            j = j + 1;
        
        
        while (i != 0 && j != 7) {//Going bottom up
            zigZagArray.push_back(block8x8[i][j]);
            i = i - 1;
            j = j + 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take edge element
        if (j<7)//If we didn't hit the edge, increment j
            j = j + 1;
        
        else if (j == 7)//If we hit the last element, go down one place
            i = i + 1;
        
        if (i >= 7 && j >= 7)//If we hit last element matrix[8][8] exit
            break;
    }
} // end ZigZagCoding


void jpeg_encoder::ZigZagCoding(double block8x8[8][8], vector<char>&zigZagArray) {
    //k- is zigZagArray index, i,j are index of matrix
    zigZagArray.push_back(block8x8[0][0]);//Take the first element
    int i = 0, j = 1;//Define index for matrix
    while (1) {
        while (j != 0 && i != 7) {//Going upside down until j!=0
            zigZagArray.push_back(block8x8[i][j]);
            i = i + 1;
            j = j - 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take the edge element
        
        if (i<7)//If not last row, increment i
            i = i + 1;
        
        else if (i == 7)//If we hit the last row, we go right one place
            j = j + 1;
        
        
        while (i != 0 && j != 7) {//Going bottom up
            zigZagArray.push_back(block8x8[i][j]);
            i = i - 1;
            j = j + 1;
        }
        zigZagArray.push_back(block8x8[i][j]);//Take edge element
        if (j<7)//If we didn't hit the edge, increment j
            j = j + 1;
        
        else if (j == 7)//If we hit the last element, go down one place
            i = i + 1;
        
        if (i >= 7 && j >= 7)//If we hit last element matrix[8][8] exit
            break;
    }
} // end ZigZagCoding


void jpeg_encoder::ZigZagCoding(vector<vector<double> > &block8x8, vector<int>&zigZagArray) {
    
    //k- is zigZagArray index, i,j are index of matrix
    zigZagArray.push_back(static_cast<int>(block8x8[0][0]));//Take the first element
    int i = 0, j = 1;//Define index for matrix
    while (1) {
        while (j != 0 && i != 7) {//Going upside down until j!=0
            zigZagArray.push_back(static_cast<int>(block8x8[i][j]));
            i = i + 1;
            j = j - 1;
        }
        zigZagArray.push_back(static_cast<int>(block8x8[i][j]));//Take the edge element
        
        if (i<7)//If not last row, increment i
            i = i + 1;
        
        else if (i == 7)//If we hit the last row, we go right one place
            j = j + 1;
        
        
        while (i != 0 && j != 7) {//Going bottom up
            zigZagArray.push_back(static_cast<int>(block8x8[i][j]));
            i = i - 1;
            j = j + 1;
        }
        zigZagArray.push_back(static_cast<int>(block8x8[i][j]));//Take edge element
        if (j<7)//If we didn't hit the edge, increment j
            j = j + 1;
        
        else if (j == 7)//If we hit the last element, go down one place
            i = i + 1;
        
        if (i >= 7 && j >= 7)//If we hit last element matrix[8][8] exit
            break;
    }
} // end Zigzag coding

void jpeg_encoder::write_baseline_dct_info(ofstream &file) {
    
    /*SOFn: Start of frame marker - Marks the beginning of the frame parameters. The subscript n identifies whether
     *the encoding process is baseline sequential, extended sequential, progressive, or lossless, as well as which
     *entropy encoding procedure is used.*/
    /*For Baseline DCT process, SOFn frame marker is FFC0*/
    /*Here is general scheme of SOFn block of data*/
    /*FFC0  Lf  P  Y   X   Nf  [C1 H1 V1 Tq1] [C2 H2 V2 Tq2] ......[Cn Hn Vn Tqn]
     *Number of bits    16   16  8  16  16  8    8  4  4   8
     *FFC0 is start of Baseline DCT marker
     *Lf - frame header length
     *P - precision
     *Y - Height
     *X - Width
     *Nf - Number of components (For our case is 3 (YCbCr))
     *C1 - Component ID
     *H1 - Horisontal sampling factor (usind in chroma subsamping)
     *V1 - Vertical sampling factor (usind in chroma subsamping)
     *Tq1- Quantization table ID used for C1 component*/
    
    emit_marker(M_SOF0, file);
    
    // Lf = 17 bytes
    int comp_size = static_cast<int>(jpegDecoder->components.size());
    
    if (comp_size > 1) {
        emit_word(0x11, file);
    }
    else {
        emit_word(0x0B, file);
    }
    
    // Precision
    emit_byte(0x08, file);
    
    // Y is 16 bits long. I keep height in int variable
    // so I need to do some calculations
    // Height and then width:
    emit_word((image_to_export_height & 0x0000FFFF), file);
    
    
    
    // X is 16 bits long. I keep height in int variable
    // so I need to do some calculations
    // Height and then width:
    emit_word((image_to_export_width & 0x0000FFFF), file);
    
    // number of components (Nf)
    emit_byte(jpegDecoder->components.size(), file);
    
    
    for (uint iComponent = 0; iComponent < jpegDecoder->components.size(); ++iComponent) {
        
        // Emit the component ID
        emit_byte(1 + iComponent, file);
        
        // Emit the subsampling
        int hFactor = jpegDecoder->components[iComponent].HFactor;
        int vFactor = jpegDecoder->components[iComponent].VFactor;
        uint_8 subsampling_byte = ((hFactor & 0xF) << 4)
        | (vFactor & 0xF);
        emit_byte(subsampling_byte, file);
        
        // Emit the quantization table ID
        QuantizationTable* qTable = jpegDecoder->components[iComponent].componentQuantizationTable;
        emit_byte(qTable->tableID, file);
        
    }
    
}


// copy & paste the header
// Writes the header part from the original picture until SOS

void jpeg_encoder::writeHeaderFromOriginalPicture(ofstream &file) {
    // Init words in header buffer
    num_bytes_in_jpeg_enc_write_buffer = 0;
    
    // previous marker
    prev_marker = 0;
    
    // current marker
    marker = 0;
    string orginal_fileName = jpegDecoder->jpeg_filename;
    FILE * fp_enc = fopen(orginal_fileName.c_str(), "rb");
    
    // counter to check the remaining length of the reserved application header length
    counter_FFEX = jpegDecoder->application_size + 2;
    
    if (fp_enc) {
        while (parseSegEnc(fp_enc, file));
        fclose(fp_enc);
    }
    else {
        perror("JPEG write encoder error");
    }
    
    //    cout << "Before flushing the number of bytes are " << num_bytes_in_jpeg_enc_write_buffer << endl;
    //    num_bytes_in_jpeg_enc_write_buffer--;
    //    flush_jpeg_enc_buffer(file);
    
#if IS_JPEG_ENCODER_WRITE_FAST
    // Seek to avoid 0xFFDA written twice for SOS marker
    num_bytes_in_jpeg_enc_write_buffer--;
    // Your number of bytes in the buffer cannot be less than zero under any case
    if(num_bytes_in_jpeg_enc_write_buffer < 0) {
        num_bytes_in_jpeg_enc_write_buffer = 0;
    }
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
    
    write_baseline_dct_info(file);
    
#if CUSTOMIZED_HUFFMAN_TABLE
    build_customized_huffman_tables();
#else
    build_default_huffman_tables();
#endif

    write_default_huffman_tables(file);
    
#elif IS_ENABLE_CUSTOMIZED_HUFF_TABLES
    write_baseline_dct_info(file);
    build_customized_huffman_tables();
    write_default_huffman_tables(file);
    
#else
    if (progressive_Huff_Format) {
        write_baseline_dct_info(file);
        build_default_huffman_tables();
        build_customized_huffman_tables();
        write_default_huffman_tables(file);
    }
#endif
    
    
#else
    // Seek to avoid 0xFFDA written twice for SOS marker
    num_bytes_in_jpeg_enc_write_buffer--;
    
#if DEBUGLEVEL > 30
    printf("2- Byte last is %X \n", jpeg_enc_write_buffer[num_bytes_in_jpeg_enc_write_buffer]);
    cout << "Current number of bytes " << num_bytes_in_jpeg_enc_write_buffer << endl;
#endif
    flush_jpeg_enc_buffer(file);
    
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
    write_baseline_dct_info(file);
//    build_default_huffman_tables();
    build_customized_huffman_tables();
    write_default_huffman_tables(file);
#else
    if (progressive_Huff_Format) {
        write_baseline_dct_info(file);
//        build_default_huffman_tables();
        build_customized_huffman_tables();
        write_default_huffman_tables(file);
    }
#endif
    
    
#endif
    
}

int jpeg_encoder::parseSegEnc(FILE * fp, ofstream &file) {
    int comp_size = static_cast<int>(jpegDecoder->components.size());
    if (!fp) {
        printf("File failed to open.\n");
        return JPEG_SEG_ERR;
    }
    
    // Read a byte
    prev_marker = marker;
    marker = fgetc(fp);
    uint_16 real_marker = prev_marker << 8 | marker;
    
#if DEBUGLEVEL > 40
    long fpos = ftell(fp);
    printf("Reading marker: %X at %d \n", marker, fpos - 1);
#endif
    
    switch (real_marker) {
            // For progressive we stop here to deal with it with our default sequential writing
#if IS_DEFAULT_QTABLE
        case 0xFFDB:
            
            if (counter_FFEX > 0) {
                add_byte_to_jpeg_enc_buffer(marker, file);
                counter_FFEX--;
            }
            else {
                if (counter_FFDB == 0) {
                    // Write the DQT:
                    emit_DQT(file);
                    counter_FFDB++;
                    uint_8 marker_skip = marker;
                    uint_16 real_marker_skip = -1;
                    int number_to_skip = -1;
                    
                    if (comp_size > 1) {
                        number_to_skip = SKIP_BYTES_Q_FACTOR_EXP_RGB;
                    }
                    else {
                        number_to_skip = SKIP_BYTES_Q_FACTOR_EXP_BLKANDWHITE;
                    }
                    
                    for(int counter_skip_bytes = 0; counter_skip_bytes <= number_to_skip; counter_skip_bytes++){
                        uint_8 prev_marker_skip = marker_skip;
                        marker_skip = fgetc(fp);
                        real_marker_skip = prev_marker_skip << 8 | marker_skip;
                        //cout <<std::hex<< real_marker_skip << endl;
                        //cout << ftell(fp) << endl;
                    }
                }
            }
            break;
#endif
            
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
            // default huffman tables in the sequential and progressive modes
        case 0xFFC0:
#endif
        case 0xFFC2:
            if (counter_FFEX > 0) {
                add_byte_to_jpeg_enc_buffer(marker, file);
                counter_FFEX--;
                break;
            }
            else {
                return 0;
            }
            // SOS marker
        case 0xFFDA:
            if (counter_FFEX > 0) {
                add_byte_to_jpeg_enc_buffer(marker, file);
                counter_FFEX--;
                break;
            }
            else {
                return 0;
            }
        default:
            add_byte_to_jpeg_enc_buffer(marker, file);
            counter_FFEX--;
            break;
    }
    return JPEG_SEG_OK;
}


void jpeg_encoder::add_byte_to_jpeg_enc_buffer(uint_8 byte, ofstream &file) {
    
    
#if DEBUGLEVEL > 40
    static int countBytes = 0;
    countBytes++;
    
    int byte_debug = 600;
    if (countBytes >= byte_debug) {
        
        printf("New Input Byte: %X, ", byte);
        cout << "Number of Bytes in the buffer is: " << num_bytes_in_jpeg_enc_write_buffer << endl;
        
        //for (int i = byte_debug - 4; i <= num_bytes_in_jpeg_enc_write_buffer; ++i) {
        for (int i = byte_debug - 4; i <= byte_debug + 12; ++i) {
            printf(" %X, ", jpeg_enc_write_buffer[i]);
        }
        
        cout << "\n Done with showing the status of the buffer " << endl;
        getchar();
        
    }
#endif
    
    //    cout << "Number bytes: " << num_bytes_in_jpeg_enc_write_buffer << endl;
    if (num_bytes_in_jpeg_enc_write_buffer < JPEG_OUT_HEADER_SIZE) {
        jpeg_enc_write_buffer[num_bytes_in_jpeg_enc_write_buffer++] = byte;
        
        if (num_bytes_in_jpeg_enc_write_buffer >= JPEG_OUT_HEADER_SIZE) {
            write_jpeg_enc_buffer(file, JPEG_OUT_HEADER_SIZE);
        } // end inner if
    } // end outer if
    
}


void jpeg_encoder::write_jpeg_enc_buffer(ofstream &file, int numBytes) {
    if (numBytes > 0 && numBytes <= JPEG_OUT_HEADER_SIZE) {
        file.write((char*)jpeg_enc_write_buffer, numBytes);
        file.flush();
        num_bytes_in_jpeg_enc_write_buffer = 0;
    }
    else {
        cout << "[ERROR WRITE] numBytes is " << numBytes << ", which is not within the correct size " << endl;
    }
}

void jpeg_encoder::flush_jpeg_enc_buffer(ofstream &file) {
    
    // Stop if you are less than or equal 0
    if(num_bytes_in_jpeg_enc_write_buffer <= 0)
        return;
    
    if (num_bytes_in_jpeg_enc_write_buffer > 0 && num_bytes_in_jpeg_enc_write_buffer <= JPEG_OUT_HEADER_SIZE) {
        file.write((char*)jpeg_enc_write_buffer, num_bytes_in_jpeg_enc_write_buffer);
        file.flush();
        num_bytes_in_jpeg_enc_write_buffer = 0;
    }
    else {
        cout << "[ERROR FLUSH] numBytes is " << num_bytes_in_jpeg_enc_write_buffer << ", which is not within the correct size " << endl;
    }
}

void jpeg_encoder::writeQuantizationTablesInFile(ofstream &file, vector<char> &table, int tableID) {
    /*Quantization table-specification syntax:
     *                DQT    Lq   Pq  Tq  [Q0  Q1  Q2  ...  Q63]
     *Number of bits: 16     16    4  4    8   8   8          8
     *
     *
     *   DQT - Define Quantization Table marker - Marks the beginning of quantization table-specification parameters (FFDB)
     *   Lq - Quantization table definition length - Specifies the length of all quantization table parameters
     *
     *   Pq - Quantization table element precision - Specifies the precision of the Qk values. Value 0 indicates 8-bit Qk
     *   values; value 1 indicates 16-bit Qk values. Pq shall be zero for 8 bit sample precision
     *   Tq: Quantization table destination identifier - Specifies one of four possible destinations at the decoder into
     *   which the quantization table shall be installed.
     *   Qk: Quantization table element - Specifies the kth element out of 64 elements, where k is the index in the zigzag
     *    ordering of the DCT coefficients. The quantization elements shall be specified in zig-zag scan order.
     */
    char a = (char)0xFF;
    file.write((char*)&a, 1);
    a = (char)0xDB;
    file.write((char*)&a, 1);
    a = 0x00;
    file.write((char*)&a, 1);
    a = 0x43;
    file.write((char*)&a, 1);
    if (tableID == 0) {
        a = 0x00;
        file.write((char*)&a, 1);
    }
    else {
        a = 0x01;
        file.write((char*)&a, 1);
    }
    for (unsigned int i = 0; i<table.size(); i++) {
        file.write((char*)&table[i], 1);
    }
}


void jpeg_encoder::emit_DQT(ofstream &file) {
    /*Quantization table-specification syntax:
     *                DQT    Lq   Pq  Tq  [Q0  Q1  Q2  ...  Q63]
     *Number of bits: 16     16    4  4    8   8   8          8
     *
     *
     *   DQT - Define Quantization Table marker - Marks the beginning of quantization table-specification parameters (FFDB)
     *   Lq - Quantization table definition length - Specifies the length of all quantization table parameters
     *
     *   Pq - Quantization table element precision - Specifies the precision of the Qk values. Value 0 indicates 8-bit Qk
     *   values; value 1 indicates 16-bit Qk values. Pq shall be zero for 8 bit sample precision
     *   Tq: Quantization table destination identifier - Specifies one of four possible destinations at the decoder into
     *   which the quantization table shall be installed.
     *   Qk: Quantization table element - Specifies the kth element out of 64 elements, where k is the index in the zigzag
     *    ordering of the DCT coefficients. The quantization elements shall be specified in zig-zag scan order.
     */
    
    //emit_marker(M_DQT, file);
    add_byte_to_jpeg_enc_buffer(M_DQT, file);
    int comp_size = static_cast<int>(jpegDecoder->components.size());
    
    // Current counter
    int counter = COMPONENT_Y;
    
    // Two quantization tables at max for {Y}, {Cb, Cr}
    do
    {
        /*QuantizationTable * qTable = jpegDecoder->components[counter].componentQuantizationTable;
         int tableLength = qTable->tableLength;
         */
        int tableLength = 0x43;
        if (counter == 0) {
            add_byte_to_jpeg_enc_buffer(0, file);
            add_byte_to_jpeg_enc_buffer(tableLength, file);
        }
        else{
            add_byte_to_jpeg_enc_buffer(0xFF, file);
            add_byte_to_jpeg_enc_buffer(M_DQT, file);
            add_byte_to_jpeg_enc_buffer(0, file);
            add_byte_to_jpeg_enc_buffer(tableLength, file);
        }
        
        // Write the tableID
        uint_8 tableID = counter;
        add_byte_to_jpeg_enc_buffer(tableID, file);
        
        // Write the QTable values
        vector <char> qtable_vec;
        
        if (counter == COMPONENT_Y) {
            ZigZagCoding(quantization_table_write_process_luminance.quantizationTableData, qtable_vec);
        }
        else {
            ZigZagCoding(quantization_table_write_process_chrominance.quantizationTableData, qtable_vec);
        }
        
        for (int i = 0; i < qtable_vec.size(); ++i) {
            add_byte_to_jpeg_enc_buffer(qtable_vec.at(i), file);
        }
        qtable_vec.clear();
        
        // Increment the counter
        counter++;
    } while (counter < comp_size - 1);
    
} // end write_dqt


void jpeg_encoder::writeStartOfFileByteInFile(ofstream &file) {
    emit_marker(M_SOI, file);
}


// Write huffman tables from the decoder:
void jpeg_encoder::write_default_huffman_tables(ofstream &file) {
    int comp_size = static_cast<int>(jpegDecoder->components.size());
    // First, calculate the table length (Lh)
    uint_16 table_length = 0;
    
    // Right most 4 bits of the byte
    uint_8  tableID = 0; // Specifies one of component: 0 for luminance and 1 for chrominance
    
    
    // Left most 4 bits of the byte
    uint_8  tableClass = 0; // Specifies is it DC element or AC element of table. 0-DC element 1-AC element
    uint_8  huffmanTableOptions = 0; // I will decompose this in two nibbles
    
    // Read components and their huffman tables
    for (int i = 0; i < default_huffmanTables.size(); ++i) {
        if (comp_size <= 1 && (i == 1 || i == 3) ) {
            continue;
        }
        
        const HuffmanTable* hTable = default_huffmanTables.at(i);
        
        // Emit the DHT marker
        emit_marker(M_DHT, file);
        
        // table length:
        table_length = hTable->tableSegmentLengthFromBitstream;
        emit_word(table_length, file);
        
        
        tableID = hTable->tableID;
        tableClass = hTable->tableClass;
        huffmanTableOptions = (tableClass << 4) | (tableID);
        emit_byte(huffmanTableOptions, file);
        
        // Next 16 bytes are number of elements coded with 1-16 bits
        for (int j = 0; j < 16; ++j)
        {
            emit_byte(hTable->number_of_codes_for_each_1to16[j], file);
        }
        
        // Remaining bytes are the data values to be mapped
        // Build the Huffman map of (length, code) -> value
        // Once the map has been built, emit it out
        std::map<huffKey, uint_8>::const_iterator iter_new;
        for (iter_new = hTable->huffData.begin(); iter_new != hTable->huffData.end(); ++iter_new) {
            
            // Print Code - Its Length : Equivalent letter
            uint_8 element = iter_new->second;
            emit_byte(element, file);
        }
    }
}


void jpeg_encoder::writeEOFMarker(ofstream &file) {
    
    emit_marker(M_EOI, file);
    
#if    IS_JPEG_ENCODER_WRITE_FAST
#if DEBUGLEVEL > 30
    cout << "FLUSHING AT THE END OF MARKER" << endl;
#endif
    
    flush_jpeg_enc_buffer(file);
#endif
    
}

void jpeg_encoder::emit_start_markers(ofstream &file) {
    
    // Write Start of File marker
    writeStartOfFileByteInFile(file);
    
    // Write the appliation signature:
    write_jfif_app0(file);
    
    // Write the DQT:
    emit_DQT(file);
    
    // Write the SOF:
    emit_sof(file);
    
} // end emit_start_markers


void jpeg_encoder::emit_sof(ofstream &file) {
    
    /*SOFn: Start of frame marker - Marks the beginning of the frame parameters. The subscript n identifies whether
     *the encoding process is baseline sequential, extended sequential, progressive, or lossless, as well as which
     *entropy encoding procedure is used.*/
    /*For Baseline DCT process, SOFn frame marker is FFC0*/
    /*Here is general scheme of SOFn block of data*/
    /*FFC0  Lf  P  Y   X   Nf  [C1 H1 V1 Tq1] [C2 H2 V2 Tq2] ......[Cn Hn Vn Tqn]
     *Number of bits    16   16  8  16  16  8    8  4  4   8
     *FFC0 is start of Baseline DCT marker
     *Lf - frame header length
     *P - precision
     *Y - Height
     *X - Width
     *Nf - Number of components (For our case is 3 (YCbCr))
     *C1 - Component ID
     *H1 - Horisontal sampling factor (usind in chroma subsamping)
     *V1 - Vertical sampling factor (usind in chroma subsamping)
     *Tq1- Quantization table ID used for C1 component*/
    
    // TODO: getters and setters
    uint_8 m_num_components = static_cast<uint_8>(jpegDecoder->components.size());
    uint_8 precision = jpegDecoder->jpegImageSamplePrecision;
    uint_16 height = image_to_export_height;
    uint_16 width = image_to_export_width;
    
    emit_marker(M_SOF0, file);                           /* baseline */
    emit_word(3 * m_num_components + 2 + 5 + 1, file); // length of the header
    emit_byte(precision, file);                                  /* precision */
    emit_word(height, file);
    emit_word(width, file);
    emit_byte(m_num_components, file);
    for (int i = 0; i < m_num_components; i++) {
        emit_byte(static_cast<uint_8>(i + 1), file);                                   /* component ID     */
        
        uint_8 HFactor = static_cast<uint_8>(jpegDecoder->components.at(i).HFactor);
        uint_8 VFactor = static_cast<uint_8>(jpegDecoder->components.at(i).VFactor);
        uint_8 qTableID = static_cast<uint_8>(jpegDecoder->components.at(i).componentQuantizationTable->tableID);
        emit_byte((HFactor << 4) + VFactor, file);  /* h and v sampling */
        
        emit_byte(qTableID, file);     /* quant. table num */
        
        
    }
    
}


void jpeg_encoder::emit_sos(ofstream &file) {
    
    /*Scan header looks like:
     *
     *
     *                    SOS   Lh  Ns [Cs1 Td1 Ta1][Cs2 Td2 Ta2]...[Csn Tdn Tan]  Ss Se Ah Al
     *
     *Number of bits:     16    16  8    8   4    4   8   4   4       8   4   4    8   8  4  4
     *
     *SOS: Start of scan marker - Marks the beginning of the scan parameters.
     *Lh :Scan header length - Specifies the length of the scan header
     *Ns: Number of image components in scan - Specifies the number of source image components in the scan. The
     *value of Ns shall be equal to the number of sets of scan component specification parameters (Csj, Tdj, and Taj)
     *present in the scan header.
     *
     *Csj - Scan component selector - Selects which of the Nf image components specified in the frame parameters
     *shall be the jth component in the scan. Each Csj shall match one of the Ci values specified in the frame header,
     *and the ordering in the scan header shall follow the ordering in the frame header. If Ns > 1, the order of
     *interleaved components in the MCU is Cs1 first, Cs2 second, etc.
     *
     *Tdj: DC entropy coding table destination selector - Specifies one of four possible DC entropy coding table
     **destinations from which the entropy table needed for decoding of the DC coefficients of component Csj is
     *retrieved
     *
     **Taj: AC entropy coding table destination selector - Specifies one of four possible AC entropy coding table
     *destinations from which the entropy table needed for decoding of the AC coefficients of component Csj is
     *retrieved
     *
     *Ss: Start of spectral or predictor selection - In the DCT modes of operation, this parameter specifies the first
     *DCT coefficient in each block in zig-zag order which shall be coded in the scan. This parameter shall be set to
     *zero for the sequential DCT processes
     *
     *Se: End of spectral selection - Specifies the last DCT coefficient in each block in zig-zag order which shall be
     **coded in the scan. This parameter shall be set to 63 for the sequential DCT processes. In the lossless mode of
     *operations this parameter has no meaning. It shall be set to zero.
     *
     *Ah: Successive approximation bit position high - This parameter specifies the point transform used in the
     *preceding scan (i.e. successive approximation bit position low in the preceding scan) for the band of coefficients
     **specified by Ss and Se. This parameter shall be set to zero for the first scan of each band of coefficients. In the
     *lossless mode of operations this parameter has no meaning. It shall be set to zero.
     *
     *Al: Successive approximation bit position low or point transform - In the DCT modes of operation this
     ***parameter specifies the point transform, i.e. bit position low, used before coding the band of coefficients
     *specified by Ss and Se. This parameter shall be set to zero for the sequential DCT processes. In the lossless
     *mode of operations, this parameter specifies the point transform, Pt.
     */
    
    // TODO: getters and setters
    uint_8 m_num_components = static_cast<uint_8>(jpegDecoder->components.size());
    
    emit_marker(M_SOS, file); // marker
    emit_word(2 * m_num_components + 2 + 1 + 3, file); // header length
    emit_byte(m_num_components, file); // number of components
    for (int i = 0; i < m_num_components; i++) {
        emit_byte(static_cast<uint_8>(i + 1), file); // component ID
        
        // TODO: getters
        uint_8 tableDC, tableAC;
        
        
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
        if (i == COMPONENT_Y) {
            tableDC = default_huffmanTables.at(Y_DC_IDX)->tableID;
            tableAC = default_huffmanTables.at(Y_AC_IDX)->tableID;
        }
        else {
            tableDC = default_huffmanTables.at(CbCr_DC_IDX)->tableID;
            tableAC = default_huffmanTables.at(CbCr_AC_IDX)->tableID;
        }
        
#else
        if (!progressive_Huff_Format) {
            // Sequential Mode use original picture Huffman Table
            tableDC = jpegDecoder->componentTablesDC[i]->tableID;
            tableAC = jpegDecoder->componentTablesAC[i]->tableID;
        }
        else {
            // Sequential Mode use original picture Huffman Table
            if (i == COMPONENT_Y) {
                tableDC = default_huffmanTables.at(Y_DC_IDX)->tableID;
                tableAC = default_huffmanTables.at(Y_AC_IDX)->tableID;
            }
            else {
                tableDC = default_huffmanTables.at(CbCr_DC_IDX)->tableID;
                tableAC = default_huffmanTables.at(CbCr_AC_IDX)->tableID;
            }
        }
#endif
        
        if (i == 0) {
            emit_byte((tableDC << 4) + tableAC, file); // component huffman table (left part is DC, right part is AC)
            //            emit_byte((0 << 4) + 0, file); // component huffman table (left part is DC, right part is AC)
        }
        else {
            emit_byte((tableDC << 4) + tableAC, file); // component huffman table (left part is DC, right part is AC)
            //            emit_byte((1 << 4) + 1, file);
        }
    }
    
    // TODO: getters
    //uint_8 zigStart = jpegDecoder->zigZagStart;
    //uint_8 zigEnd = jpegDecoder->zigZagEnd;
    
    // Even progressive encoded as sequential method
    uint_8 zigStart = 0;
    uint_8 zigEnd = 63;
    uint_8 dummyByte = 0;
    
    emit_byte(zigStart, file);     /* spectral selection */
    emit_byte(zigEnd, file);
    emit_byte(dummyByte, file); // TODO: Bit approximation for progressive JPEG
}


void jpeg_encoder::write_jfif_app0(ofstream &file) {
    
    emit_marker(M_APP0, file);
    emit_word(2 + 4 + 1 + 2 + 1 + 2 + 2 + 1 + 1, file);
    emit_byte(0x4A, file); emit_byte(0x46, file); emit_byte(0x49, file); emit_byte(0x46, file); /* Identifier: ASCII "JFIF" */
    emit_byte(0, file);
    emit_byte(1, file);      /* Major version */
    emit_byte(1, file);      /* Minor version */
    emit_byte(0, file);      /* Density unit */
    emit_word(1, file);
    emit_word(1, file);
    emit_byte(0, file);      /* No thumbnail image */
    emit_byte(0, file);
}

void jpeg_encoder::encodeImageEntryPoint(vector<int> const &luminanceZigZagArray,
                                         vector<int> const &chrominanceCbZigZagArray, vector<int> const &chrominanceCrZigZagArray, ofstream &file) {
    
    
    // Y is the original hFactor and vFactor; others are subsampled
    int hFactor = jpegDecoder->components[COMPONENT_Y].HFactor;
    int vFactor = jpegDecoder->components[COMPONENT_Y].VFactor;
    
    
    // MCU: minimum coding unit
    int xstride_by_mcu = 8 * hFactor; // mcu width
    int ystride_by_mcu = 8 * vFactor; // mcu height
    
    // Just encode the image by 'macroblock' (size is 8x8, 8x16, or 16x16)
    for (int y = 0; y < image_to_export_height; y += ystride_by_mcu)
    {
        for (int x = 0; x < image_to_export_width; x += xstride_by_mcu)
        {
            encode_mcu(luminanceZigZagArray, chrominanceCbZigZagArray, chrominanceCrZigZagArray, hFactor, vFactor, x, y, file);
        }
        
    }
    
    // append the last few bits that are left
    if (m_bits_in > 0 && m_bits_in < 8) {
        uint_8 byte_last = m_bit_buffer >> (32 - m_bits_in);
        byte_last <<= (8 - m_bits_in);
        byte_last |= ((1 << (8 - m_bits_in)) - 1);
        
        emit_byte(byte_last, file);
        if (byte_last == 0xFF) {
            emit_byte(0, file); // Byte Stuffing
        }
    }
}

void jpeg_encoder::encode_mcu(vector<int> const &luminanceZigZagArray, vector<int> const &chrominanceCbZigZagArray, vector<int> const &chrominanceCrZigZagArray, int componentWidth, int componentHeight, int start_x, int start_y, ofstream &file) {
    //cout << "height:" << componentHeight << " Width: " << componentWidth << endl;
    for (int y = 0; y < componentHeight; ++y)
    {
        for (int x = 0; x < componentWidth; ++x)
        {
            
            int luma_x = start_x + x * 8;
            int luma_y = start_y + y * 8;
            
#if PROFILE_JPEG_SAVE_PIC > 5
            auto sEncodeStartTime = std::chrono::high_resolution_clock::now();
#endif
            // encode block
            encode_block(luminanceZigZagArray, luma_x, luma_y, COMPONENT_Y, count_block_Y, file);
            count_block_Y++;
            
#if PROFILE_JPEG_SAVE_PIC > 20
            if(count_block_Y % 100 == 0) {
                auto sEncodeEndTime = std::chrono::high_resolution_clock::now();
                cout << "Encoding of 100 Y-blocks elapsed time: " <<  std::chrono::duration_cast<std::chrono::milliseconds>(sEncodeEndTime - sEncodeStartTime).count() << " milliseconds" << endl;
            }
#endif
            
            
        }
        
    }
    
    int numberOfComponents = static_cast<int>(jpegDecoder->components.size());
    // The rest of the components if they exist:
    for (int iComponent = 1; iComponent < numberOfComponents; ++iComponent)
    {
        
        int chroma_x = start_x / componentWidth;
        int chroma_y = start_y / componentHeight;
        
        // encode block
        if (iComponent == COMPONENT_Cb) {
            
            encode_block(chrominanceCbZigZagArray, chroma_x, chroma_y, COMPONENT_Cb, count_block_Cb, file);
            count_block_Cb++;
        }
        else if (iComponent == COMPONENT_Cr) {
            
            encode_block(chrominanceCrZigZagArray, chroma_x, chroma_y, COMPONENT_Cr, count_block_Cr, file);
            count_block_Cr++;
            
        }
    }
}

void jpeg_encoder::encode_block(vector<int> const &zigZagArray, int CurrentX, int CurrentY, int currentComponent, int count_block, ofstream &file) {
    
    // DC coding
#if DEBUGLEVEL > 20
    if (currentComponent == COMPONENT_Y)
        cout << "Encode block at X: " << CurrentX << ", Y: " << CurrentY << endl;
#endif
    
    //cout << total_block << endl;
    int start_index = returnIndexInZigZagArray(count_block);
    const int dc_delta = zigZagArray.at(start_index);
    
    
#if DEBUGLEVEL > 20
    if (currentComponent == COMPONENT_Y && CurrentX >= 224 && CurrentY == 0)
    {
        int init = start_index;
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                
                int idx = start_index + j + i * 8;
                cout << zigZagArray.at(idx) << ", ";
            }
            cout << "\n";
        }
        
    }
#endif
    
#if DEBUGLEVEL > 20
    cout << "DC Delta " << dc_delta << " at start index: " << start_index << "\n" << endl;
#endif
    
    // Encode DC delta coefficient
    // 1. Write the code for the DC delta from Huffman DC table of current component with the its corresponding codeLength
    // (table->code[bit_count(dc_delta)], table->codeLength[bit_count(dc_delta)]
    // 2. Write the signed int bits for the dc_delta itself write(dc_delta, nbits)
    HuffmanTable* dc_hTable;
    
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
    if (currentComponent == COMPONENT_Y) {
        dc_hTable = default_huffmanTables.at(Y_DC_IDX);
    }
    else {
        dc_hTable = default_huffmanTables.at(CbCr_DC_IDX);
    }
#elif IS_ENABLE_COSUTMIZE_HUFF_TABLES
    
#else
    if (!progressive_Huff_Format ) {
        dc_hTable = jpegDecoder->componentTablesDC[currentComponent];
    }
    else {
        if (currentComponent == COMPONENT_Y) {
            dc_hTable = default_huffmanTables.at(Y_DC_IDX);
        }
        else {
            dc_hTable = default_huffmanTables.at(CbCr_DC_IDX);
        }
    }
#endif
    
    const uint nbits = bit_count(dc_delta);
    // put_bits:
    put_bits(dc_hTable->codes[nbits], dc_hTable->codeLengths[nbits], file);
	if (Y_only)
	{
		if ( currentComponent == COMPONENT_Y )
			counter_bits_buffer += dc_hTable->codeLengths[nbits];
	}
	else
		counter_bits_buffer += dc_hTable->codeLengths[nbits];
	
	put_signed_int_bits(dc_delta, nbits, file);
    
	if (Y_only)
	{
		if (currentComponent == COMPONENT_Y)
			counter_bits_buffer += nbits;
	}
	else 
		counter_bits_buffer += nbits;
	// remove
    // cout << dc_hTable->codes[nbits] << "----" << dc_hTable->codeLengths[nbits] << endl;
    // cout << dc_delta << "====" << nbits << endl;
#if DEBUGLEVEL >20
    long pos = file.tellp();
    if (pos >= 5555522) {
        cout << pos << endl;
        cout << "# of bits:" << nbits << endl;
        cout << dc_hTable->codes[nbits] << "--" << dc_hTable->codeLengths[nbits] << endl;
        cout << "postion X:" << CurrentX << ", position Y:" << CurrentY << endl;
    }
#endif
    // ----------
    
    // Encode AC coefficients:
    HuffmanTable* ac_hTable;
    
    
#if IS_ENABLE_USE_DEFAULT_HUFF_TABLES
    if (currentComponent == COMPONENT_Y) {
        ac_hTable = default_huffmanTables.at(Y_AC_IDX);
    }
    else {
        ac_hTable = default_huffmanTables.at(CbCr_AC_IDX);
    }
#else
    
    if (!progressive_Huff_Format) {
        ac_hTable = jpegDecoder->componentTablesAC[currentComponent];
    }
    else {
        if (currentComponent == COMPONENT_Y) {
            ac_hTable = default_huffmanTables.at(Y_AC_IDX);
        }
        else {
            ac_hTable = default_huffmanTables.at(CbCr_AC_IDX);
        }
    }
#endif
    
    int run_len = 0;
    
    for (int i = start_index + 1; i < start_index + 64; i++) {
        const short ac_val = zigZagArray.at(i);
        
        if (ac_val == 0) {
            run_len++;
        }
        else {
            // 16 zeros case
            while (run_len >= 16)
            {
                // Write bits (Huffman code of 0xF0 and it's code Length 0XF0) for 16 zeros
                // put_bits:
                put_bits(ac_hTable->codes[0xF0], ac_hTable->codeLengths[0xF0], file);
                
				if (Y_only)
				{
					if (currentComponent == COMPONENT_Y)
						counter_bits_buffer += ac_hTable->codeLengths[0xF0];
				}
				else 
					counter_bits_buffer += ac_hTable->codeLengths[0xF0];
			
				run_len -= 16;
            }
            const uint nbits = bit_count(ac_val);
            const int code = (run_len << 4) + nbits; // left part is the run, and right part is the nBits of the next AC
            
            // Write HAC_codes[code], HAC_codeLength[code]
            // put_bits:
            put_bits(ac_hTable->codes[code], ac_hTable->codeLengths[code], file);
			if (Y_only)
			{
				if (currentComponent == COMPONENT_Y)
					counter_bits_buffer += ac_hTable->codeLengths[code];
			}
			else 
				counter_bits_buffer += ac_hTable->codeLengths[code];
			
			put_signed_int_bits(ac_val, nbits, file);
            
			if (Y_only)
			{
				if (currentComponent == COMPONENT_Y)
				{
					counter_bits_buffer += nbits;
					counter_bits_buffer += ac_hTable->codeLengths[code];
				}
			}
			else
			{
				counter_bits_buffer += nbits;
				counter_bits_buffer += ac_hTable->codeLengths[code];
			}
            run_len = 0;
#if DEBUGLEVEL >20
            long pos = file.tellp();
            if (pos >= 522) {
                cout << pos << endl;
                cout << "# of bits:" << nbits << endl;
                cout << "AC code " << std::hex << code << endl;
                cout << ac_hTable->codes[code] << "--" << ac_hTable->codeLengths[code] << endl;
                cout << std::dec << "postion X:" << CurrentX << ", position Y:" << CurrentY << endl;
            }
#endif
        }
        
    } // end for
    
    // If there is still a run
    if (run_len) {
        // Write the EOB or 000 code from AC huffman table HAC[0] code, HAC[0] codeLength
        // put_bits:
        // EOB
        put_bits(ac_hTable->codes[0], ac_hTable->codeLengths[0], file);
        
		if (Y_only)
		{
			if (currentComponent)
				counter_bits_buffer += ac_hTable->codeLengths[0];
		}
		else
		{
				counter_bits_buffer += ac_hTable->codeLengths[0];
		}

#if DEBUGLEVEL >20
        long pos = file.tellp();
        if (pos >= 522) {
            cout << pos << endl;
            cout << "# of bits:" << nbits << endl;
            //cout << "AC code " << std::hex << code << endl;
            cout << ac_hTable->codes[0] << "--" << ac_hTable->codeLengths[0] << endl;
            cout << std::dec << "postion X:" << CurrentX << ", position Y:" << CurrentY << endl;
        }
#endif
    }
    
}


int jpeg_encoder::returnIndexInZigZagArray(int count_block) {
    return 64 * count_block;
}

void jpeg_encoder::put_signed_int_bits(int bits, uint_32 bits_length, ofstream &file) {
    
    if (bits < 0) {
        bits--;
    }
    
    put_bits(bits & ((1 << bits_length) - 1), bits_length, file);
}

void jpeg_encoder::put_bits(uint bits, uint_32 bits_length, ofstream &file) {
    // Add the bits to your buffer
    uint_32 bits_cast = static_cast<uint_32>(bits);
    m_bit_buffer = m_bit_buffer | (bits_cast << (32 - (m_bits_in += bits_length)));
    
    while (m_bits_in >= 8) {
        uint_8 byte = static_cast<uint_8> (((m_bit_buffer >> 24) & 0xFF));
        
#if IS_JPEG_ENCODER_WRITE_FAST
        add_byte_to_jpeg_enc_buffer(byte, file);
#else
        emit_byte(byte, file);
#endif
        
        if (byte == 0xFF) {
#if IS_JPEG_ENCODER_WRITE_FAST
            add_byte_to_jpeg_enc_buffer(0, file);
#else
            emit_byte(0, file); // Byte Stuffing
#endif
        }
        // Left Shift the buffer one byte
        m_bit_buffer <<= 8;
        // Decrease the amount of m_bits_in  by 8
        m_bits_in -= 8;
    } // end while
}

void jpeg_encoder::build_default_huffman_tables() {
    for (int huffman_table_counter = 0; huffman_table_counter < 4; ++huffman_table_counter) { // Totally 4 tables
        HuffmanTable* table = 0; // initialize a Huffman table to process
        table = new HuffmanTable();
        if (huffman_table_counter < 2) { // DC condition, either Y or C
            // initialization
            if (huffman_table_counter == 0) table->tableID = 0;
            else table->tableID = 1; // tableClass 0 is for DC
            table->tableClass = 0;
            table->tableSegmentLengthFromBitstream = 0x1F;
            default_huffmanTables.push_back(table);
        }
        else { // AC condition, either Y or C
            if (huffman_table_counter == 2) table->tableID = 0;
            else table->tableID = 1; // tableClass 1 is for AC
            table->tableClass = 1;
            table->tableSegmentLengthFromBitstream = 0xB5;
            default_huffmanTables.push_back(table);
        }
        
        // DC total 12 categories in default table
        int category = 0;
        for (int i = 0; i < 16; ++i) {
            table->number_of_codes_for_each_1to16[i] = code_length_freq[huffman_table_counter][i];
        }
        
        uint_32 code = 0; // Huffman code which will be connected with element
        uint_8 element = 0; // Read element from file
        int counter_coeff = 0;
        // Remaining bytes are the data values to be mapped
        // Build the Huffman map of (length, code) -> value
        for (int i = 0; i < 16; ++i)
        {
            for (int j = 0; j < code_length_freq[huffman_table_counter][i]; ++j)
            {
                // element is the actual unique element in the alphabet
                if(huffman_table_counter == 0 || huffman_table_counter == 1) element = kDCSyms[counter_coeff++]; // DC's
                else if (huffman_table_counter == 2) element = kACSyms[0][counter_coeff++]; // Y-AC
                else element = kACSyms[1][counter_coeff++]; // C-AC
                
                table->codes[element] = code;
                table->codeLengths[element] = i + 1; // the length is at least 1 and at most 16
                // so far used for printing purposes
                // <Length, code> = element
                table->huffData[huffKey(i + 1, code)] = element;
                code++; //Elements on the same tree depth have code incremented by one
            }// end j
            
            // multiply by 2 for next iteration // shift
            code <<= 1;
        } // end i
    } //end huffman_table_counter
    
    
}

//////////////////////
//////// HOSSAM HOSSAMM HOSSAM HOSSAM HOSSAM
// Low-level helper functions.
template <class T> inline void clear_obj(T &obj)
{
    memset(&obj, 0, sizeof(obj));
}

// Radix sorts sym_freq[] array by 32-bit key m_key. Returns ptr to sorted values.
static inline jpeg_encoder::sym_freq *radix_sort_syms(uint num_syms, jpeg_encoder::sym_freq *pSyms0, jpeg_encoder::sym_freq *pSyms1)
{
    const uint cMaxPasses = 4;
    uint32 hist[256 * cMaxPasses]; clear_obj(hist);
    for (uint i = 0; i < num_syms; i++) {
        uint freq = pSyms0[i].m_key;
        hist[freq & 0xFF]++;
        hist[256 + ((freq >> 8) & 0xFF)]++;
        hist[256*2 + ((freq >> 16) & 0xFF)]++;
        hist[256*3 + ((freq >> 24) & 0xFF)]++;
    }
    jpeg_encoder::sym_freq *pCur_syms = pSyms0, *pNew_syms = pSyms1;
    uint total_passes = cMaxPasses;
    while ((total_passes > 1) && (num_syms == hist[(total_passes - 1) * 256])) {
        total_passes--;
    }
    for (uint pass_shift = 0, pass = 0; pass < total_passes; pass++, pass_shift += 8) {
        const uint32 *pHist = &hist[pass << 8];
        uint offsets[256], cur_ofs = 0;
        for (uint i = 0; i < 256; i++) {
            offsets[i] = cur_ofs;
            cur_ofs += pHist[i];
        }
        for (uint i = 0; i < num_syms; i++) {
            pNew_syms[offsets[(pCur_syms[i].m_key >> pass_shift) & 0xFF]++] = pCur_syms[i];
        }
        jpeg_encoder::sym_freq *t = pCur_syms;
        pCur_syms = pNew_syms;
        pNew_syms = t;
    }
    return pCur_syms;
}

// calculate_minimum_redundancy() originally written by: Alistair Moffat, alistair@cs.mu.oz.au, Jyrki Katajainen, jyrki@diku.dk, November 1996.
static void calculate_minimum_redundancy(jpeg_encoder::sym_freq *A, int n)
{
    int root, leaf, next, avbl, used, dpth;
    if (n==0) {
        return;
    } else if (n==1) {
        A[0].m_key = 1;
        return;
    }
    A[0].m_key += A[1].m_key;
    root = 0;
    leaf = 2;
    for (next=1; next < n-1; next++) {
        if (leaf>=n || A[root].m_key<A[leaf].m_key) {
            A[next].m_key = A[root].m_key;
            A[root++].m_key = next;
        } else {
            A[next].m_key = A[leaf++].m_key;
        }
        if (leaf>=n || (root<next && A[root].m_key<A[leaf].m_key)) {
            A[next].m_key += A[root].m_key;
            A[root++].m_key = next;
        } else {
            A[next].m_key += A[leaf++].m_key;
        }
    }
    A[n-2].m_key = 0;
    for (next=n-3; next>=0; next--) {
        A[next].m_key = A[A[next].m_key].m_key+1;
    }
    avbl = 1;
    used = dpth = 0;
    root = n-2;
    next = n-1;
    while (avbl>0) {
        while (root>=0 && (int)A[root].m_key==dpth) {
            used++;
            root--;
        }
        while (avbl>used) {
            A[next--].m_key = dpth;
            avbl--;
        }
        avbl = 2*used; dpth++; used = 0;
    }
}

// Limits canonical Huffman code table's max code size to max_code_size.
static void huffman_enforce_max_code_size(int *pNum_codes, int code_list_len, int max_code_size)
{
    int MAX_HUFF_CODESIZE = 32;
    if (code_list_len <= 1) {
        return;
    }
    
    for (int i = max_code_size + 1; i <= MAX_HUFF_CODESIZE; i++) {
        pNum_codes[max_code_size] += pNum_codes[i];
    }
    
    uint32 total = 0;
    for (int i = max_code_size; i > 0; i--) {
        total += (((uint32)pNum_codes[i]) << (max_code_size - i));
    }
    
    while (total != (1UL << max_code_size)) {
        pNum_codes[max_code_size]--;
        for (int i = max_code_size - 1; i > 0; i--) {
            if (pNum_codes[i]) {
                pNum_codes[i]--;
                pNum_codes[i + 1] += 2;
                break;
            }
        }
        total--;
    }
}

// Generates an optimized offman table.
void jpeg_encoder::huffman_table::optimize(int table_len)
{
    sym_freq syms0[MAX_HUFF_SYMBOLS], syms1[MAX_HUFF_SYMBOLS];
    syms0[0].m_key = 1; syms0[0].m_sym_index = 0;  // dummy symbol, assures that no valid code contains all 1's
    int num_used_syms = 1;
    for (int i = 0; i < table_len; i++)
        if (m_count[i]) {
            syms0[num_used_syms].m_key = m_count[i];
            syms0[num_used_syms++].m_sym_index = i + 1;
        }
    sym_freq *pSyms = radix_sort_syms(num_used_syms, syms0, syms1);
    calculate_minimum_redundancy(pSyms, num_used_syms);
    
    // Count the # of symbols of each code size.
    int num_codes[1 + MAX_HUFF_CODESIZE];
    clear_obj(num_codes);
    for (int i = 0; i < num_used_syms; i++) {
        num_codes[pSyms[i].m_key]++;
    }
    
    const uint JPGE_CODE_SIZE_LIMIT = 16; // the maximum possible size of a JPEG Huffman code (valid range is [9,16] - 9 vs. 8 because of the dummy symbol)
    huffman_enforce_max_code_size(num_codes, num_used_syms, JPGE_CODE_SIZE_LIMIT);
    
    // Compute m_huff_bits array, which contains the # of symbols per code size.
    clear_obj(m_bits);
    for (int i = 1; i <= (int)JPGE_CODE_SIZE_LIMIT; i++) {
        m_bits[i] = static_cast<uint8>(num_codes[i]);
    }
    
    // Remove the dummy symbol added above, which must be in largest bucket.
    for (int i = JPGE_CODE_SIZE_LIMIT; i >= 1; i--) {
        if (m_bits[i]) {
            m_bits[i]--;
            break;
        }
    }
    
    // Compute the m_huff_val array, which contains the symbol indices sorted by code size (smallest to largest).
    for (int i = num_used_syms - 1; i >= 1; i--) {
        m_val[num_used_syms - 1 - i] = static_cast<uint8>(pSyms[i].m_sym_index - 1);
    }
}


///// HOSSAM HOSSAM HOSSAM HOSSAM HOSSAM
//////////////////

// Compute the actual canonical Huffman codes/code sizes given the JPEG huff bits and val arrays.
int jpeg_encoder::huffman_table::compute()
{
    int last_p, si;
    uint8 huff_size[257];
    uint huff_code[257];
    uint code;
    
    int p = 0;
    for (char l = 1; l <= 16; l++)
        for (int i = 1; i <= m_bits[l]; i++) {
            huff_size[p++] = l;
        }
    
    huff_size[p] = 0; last_p = p; // write sentinel
    
    code = 0;
    si = huff_size[0];
    p = 0;
    
    while (huff_size[p]) {
        while (huff_size[p] == si) {
            huff_code[p++] = code++;
        }
        code <<= 1;
        si++;
    }

    memset(m_codes, 0, sizeof(m_codes[0])*256);
    memset(m_code_sizes, 0, sizeof(m_code_sizes[0])*256);
    for (p = 0; p < last_p; p++) {
        m_codes[m_val[p]]      = huff_code[p];
        m_code_sizes[m_val[p]] = huff_size[p];
    }

#if DEBUGLEVEL > 20
    cout << "INSIDE COMPUTE INSIDE COMPUTE" << endl;
    for( p  = 0 ; p  < last_p ; p ++)
    {
        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X\n",
               m_codes[m_val[p]], m_code_sizes[m_val[p]],  m_val[p]);
    }
    std::cout << "-----------" << std::endl;
#endif
    
    return last_p;
}


void jpeg_encoder::build_customized_huffman_tables()
{
    
    // DC LUMA:
    // -------
    
    // Create the last_p: the number of indexes that have code length greater than 0
    std::vector<int> last_p;
    last_p.resize(4, 0);
    huffman_table dc;
    
    // When you copy, you copy everything
    for(int i = 0; i < AC_LUM_CODES; ++i)
    {
        dc.m_count[i] = dc_counts_Y[i];
    }
    
    dc.optimize(DC_LUM_CODES);
    last_p[0] = dc.compute();
    
    // AC LUMA:
    // -------
    huffman_table ac;
    for(int i = 0; i < AC_LUM_CODES; ++i)
    {
        ac.m_count[i] = ac_counts_Y[i];
    }
    
    ac.optimize(AC_LUM_CODES);
    last_p[2] =  ac.compute();
    
    // DC, AC Chroma:
    // -------------
    huffman_table dc_chroma, ac_chroma;
    if(jpegDecoder->numberOfComponents > 1)
    {

        // When you copy, you copy everything
        for(int i = 0; i < AC_CHROMA_CODES; ++i)
        {
            dc_chroma.m_count[i] = dc_counts_CbCr[i];
        }
        
        for(int i = 0; i < AC_CHROMA_CODES; ++i)
        {
            ac_chroma.m_count[i] = ac_counts_CbCr[i];
        }
        
        dc_chroma.optimize(DC_CHROMA_CODES);
        last_p[1] = dc_chroma.compute();
        
        ac_chroma.optimize(AC_CHROMA_CODES);
        last_p[3] = ac_chroma.compute();
        
    } // end if
    
    ///////////// S7es ///////////

#if DEBUGLEVEL > 20
    std::cout << "compute values before copy "<< std::endl;
    for( int p  = 0 ; p  < last_p[1] ; p ++)
    {
        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X\n",
               ac.m_codes[ac.m_val[p]], ac.m_code_sizes[ac.m_val[p]],  ac.m_val[p]);
    }
    std::cout << "-----------" << std::endl;
#endif
    
    
    for (int huffman_table_counter = 0; huffman_table_counter < 4; ++huffman_table_counter) { // Totally 4 tables
        HuffmanTable* table = 0; // initialize a Huffman table to process
        table = new HuffmanTable();
        if (huffman_table_counter < 2) { // DC condition, either Y or C
            // initialization
            if (huffman_table_counter == 0) table->tableID = 0;
            else table->tableID = 1; // tableClass 0 is for DC
            table->tableClass = 0;
            table->tableSegmentLengthFromBitstream = 0x1F;
            default_huffmanTables.push_back(table);
        }
        else { // AC condition, either Y or C
            if (huffman_table_counter == 2) table->tableID = 0;
            else table->tableID = 1; // tableClass 1 is for AC
            table->tableClass = 1;
            table->tableSegmentLengthFromBitstream = 0xB5;
            default_huffmanTables.push_back(table);
        }
        // sep
        uint_32 code = 0; // Huffman code which will be connected with element
        uint_8 element = 0; // Read element from file
        int length = 0;
        
        for(int p  = 0; p  < last_p[huffman_table_counter]; ++p)
        {
            // DC LUMA
            if (huffman_table_counter == 0)
            {
                element = dc.m_val[p];
                code = dc.m_codes[element];
                length = dc.m_code_sizes[element];
            }
            // DC Chroma
            else if(huffman_table_counter == 1 && jpegDecoder->numberOfComponents > 1)
            {
                element = dc_chroma.m_val[p];
                code = dc_chroma.m_codes[element];
                length = dc_chroma.m_code_sizes[element];
                
            }
            // AC LUMA
            else if(huffman_table_counter == 2)
            {
                element = ac.m_val[p];
                code = ac.m_codes[element];
                length = ac.m_code_sizes[element];
                
            }
            // AC CHROMA
            else if(huffman_table_counter == 3 && jpegDecoder->numberOfComponents > 1)
            {
                element = ac_chroma.m_val[p];
                code = ac_chroma.m_codes[element];
                length = ac_chroma.m_code_sizes[element];
                
            }
            
            table->codeLengths[element] = length; // the length is at least 1 and at most 16
            table->huffData[huffKey(length, code)] = element;
            table->codes[element] = code;
            
        } // end p loop
        
        int segmentLength = 0;
        for(int k = 0; k < 16; ++k)
        {
            // DC LUMA
            if (huffman_table_counter == 0 ){
               table->number_of_codes_for_each_1to16[k] = dc.m_bits[k + 1];
               segmentLength += dc.m_bits[k + 1];
            }
            // DC Chroma
            else if(huffman_table_counter == 1 && jpegDecoder->numberOfComponents > 1)
            {
                table->number_of_codes_for_each_1to16[k] = dc_chroma.m_bits[k + 1];
                segmentLength += dc_chroma.m_bits[k + 1];
            }
            // AC LUMA
            else if(huffman_table_counter == 2)
            {
                table->number_of_codes_for_each_1to16[k] = ac.m_bits[k + 1];
                segmentLength += ac.m_bits[k + 1];
            }
            // AC CHROMA
            else if(huffman_table_counter == 3 && jpegDecoder->numberOfComponents > 1)
            {
             table->number_of_codes_for_each_1to16[k] = ac_chroma.m_bits[k + 1];
             segmentLength += ac_chroma.m_bits[k + 1];
            }
            
            
        } // end copy the number of codes for 1 to 16
        
        // Total length of the segment in huffman
        segmentLength += 2 + 1 + 16;
        table->tableSegmentLengthFromBitstream = segmentLength;

    } //end huffman_table_counter
    
    
#if DEBUGLEVEL > 20
    cout << "START FANKOOSH AGAIN " << endl;
    
    // Once the map has been built, print it out
    std::map<huffKey, uint_8>::iterator iter;
    HuffmanTable* table  = default_huffmanTables[2];

    if(table->tableClass)
    {
        printf("Huffman table ID #%02X class 01 (AC Table):\n", table->tableID);
    }
    else {
        printf("Huffman table ID #%02X class 00 (DC Table):\n", table->tableID);
    }

    for (iter = table->huffData.begin();iter != table->huffData.end(); ++iter) {

        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X\n",
               iter->first.second, iter->first.first, iter->second);
    }
    cout << "DONE FANKOOSH AGAIN" << endl;
    
#endif
    
} // end build_customized_huffman_table


bool jpeg_encoder::savePicture(){
    
    //    cout << "\n\n Encode Start " << endl;
    // Identify whether this is a progressive encoder or not
    progressive_Huff_Format = jpegDecoder->progressive_Huff_Format;
    
    
    // Copy the quantization tables
    //#if !IS_DEFAULT_QTABLE
    copy_qTables();
    //#endif
    
    //In these arrays I will keep my data after FDCT
    vector<int> luminanceZigZagArray;
    vector<int> chrominanceCbZigZagArray;
    vector<int> chrominanceCrZigZagArray;
    
    
    //Prosess every component with FDCT
    uint_8 ** luminance = jpegDecoder->m_YPicture_buffer;
    uint_8 ** chrominanceCb = jpegDecoder->m_CbPicture_buffer;
    uint_8 ** chrominanceCr = jpegDecoder->m_CrPicture_buffer;
    
    // NEW to TCM: Apply TCM
    //    perform_TCM();
    total = 0;
    
    // Initialize your counts
    int put_size = 256;
    ac_counts_Y.resize(put_size, 0);
    dc_counts_Y.resize(put_size, 0);
    ac_counts_CbCr.resize(put_size, 0);
    dc_counts_CbCr.resize(put_size, 0);
    
    //
    
    perform_fdct(luminance, luminanceZigZagArray, quantization_table_write_process_luminance.quantizationTableData, COMPONENT_Y);

	int numberofComponent = jpegDecoder->components.size();
    if (numberofComponent > 1)
    {
        total = 0;
        perform_fdct(chrominanceCb, chrominanceCbZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cb);
        perform_fdct(chrominanceCr, chrominanceCrZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cr);
    }


//    perform_fdct(luminance, luminanceZigZagArray, quantization_table_write_process_luminance.quantizationTableData, COMPONENT_Y);
//    int numberofComponent = jpegDecoder->components.size();
//    if (numberofComponent > 1)
//    {
//        total = 0;
//        perform_fdct(chrominanceCb, chrominanceCbZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cb);
//        perform_fdct(chrominanceCr, chrominanceCrZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cr);
//    }
    
    
#if PROFILE_JPEG_SAVE_PIC
    startTime = std::chrono::high_resolution_clock::now();
#endif
    
    // Here starts file writing process:
    string fileName = image_to_export_filename;
    ofstream output(fileName.c_str(), ios::out | ios::binary);
    
    
    // copy & paste header directly
    writeHeaderFromOriginalPicture(output);
    
    // exit(0);
    
    // SOS:
    emit_sos(output);
    encodeImageEntryPoint(luminanceZigZagArray, chrominanceCbZigZagArray, chrominanceCrZigZagArray, output);
    writeEOFMarker(output);
    image_bpp = counter_bits_buffer/double(image_to_export_height*image_to_export_width);
#if PROFILE_JPEG_SAVE_PIC
    endTime = std::chrono::high_resolution_clock::now();
    cout << "Encoding elapsed time: " <<  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " milliseconds" << endl;
#endif
    
    
    //    cout << "Encoder Done!!" << endl;
    return true;
} // end savePicture



double jpeg_encoder::writeyuv() {


	
	progressive_Huff_Format = jpegDecoder->progressive_Huff_Format;
	// Copy the quantization tables
	//#if !IS_DEFAULT_QTABLE
	copy_qTables();
	//#endif

	//In these arrays I will keep my data after FDCT
	vector<int> luminanceZigZagArray;
	vector<int> chrominanceCbZigZagArray;
	vector<int> chrominanceCrZigZagArray;


	//Prosess every component with FDCT
	uint_8** luminance = jpegDecoder->m_YPicture_buffer;
	uint_8** chrominanceCb = jpegDecoder->m_CbPicture_buffer;
	uint_8** chrominanceCr = jpegDecoder->m_CrPicture_buffer;
	
	//encodeImageEntryPoint(luminanceZigZagArray, chrominanceCbZigZagArray, chrominanceCrZigZagArray, output);
	//image_bpp = counter_bits_buffer / double(image_to_export_height * image_to_export_width);

	// NEW to TCM: Apply TCM
	//    perform_TCM();
	total = 0;

	// Initialize your counts
	int put_size = 256;
	ac_counts_Y.resize(put_size, 0);
	dc_counts_Y.resize(put_size, 0);
	ac_counts_CbCr.resize(put_size, 0);
	dc_counts_CbCr.resize(put_size, 0);

	//

	perform_fdct(luminance, luminanceZigZagArray, quantization_table_write_process_luminance.quantizationTableData, COMPONENT_Y);

	int numberofComponent = jpegDecoder->components.size();
	if (numberofComponent > 1)
	{
		total = 0;
		perform_fdct(chrominanceCb, chrominanceCbZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cb);
		perform_fdct(chrominanceCr, chrominanceCrZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cr);
	}

	int currentComponent = 0; 
	// How many Y's in the horizontal and vertical direction (2x2 is the usual case)
	int HFactor = jpegDecoder->components[currentComponent].HFactor, VFactor = jpegDecoder->components[currentComponent].VFactor;

	// These Y's should be scaled/repeated how many times horizonatally and vertically
	int HScale = jpegDecoder->components[currentComponent].HScale, VScale = jpegDecoder->components[currentComponent].VScale;

	// Note: use the upscaled width and height to store the DCT coefficientss
	int comp_height = ceil(1.0 * image_to_export_height_dct / VScale);
	int comp_width = ceil(1.0 * image_to_export_width_dct / HScale);
	int width = jpegDecoder->jpegImageWidth , height = jpegDecoder->jpegImageHeight;

	int currentX = 0;
	int currentY = 0;
	int currentBlockHFactor = 0;
	int currentBlockVFactor = 0;
	
	// here using luminanceZigZagArray , chrominanceCbZigZagArray , chrominanceCrZigZagArray that are quantized zigzag ordr of all image, we shoud
	// do decoding, each 64 values in this arrays are one zigzag order block, it is enough to do dzizag and dequantize and idct on each blok and write 
	// them into filr

	int row = 0; 
	int col = 0; 

	std::vector<std::vector<uint_8> >Y(comp_height, std::vector<uint_8>(comp_width) );
	std::vector<std::vector<uint_8> >CB(comp_height, std::vector<uint_8>(comp_width));
	std::vector<std::vector<uint_8> >CR(comp_height, std::vector<uint_8>(comp_width));
	
	// derive Y 
	write_buff(Y, luminanceZigZagArray, quantization_table_write_process_luminance.quantizationTableData , COMPONENT_Y );
	if (numberofComponent > 1)
	{
		write_buff(CB, chrominanceCbZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cb);
		write_buff(CR, chrominanceCrZigZagArray, quantization_table_write_process_chrominance.quantizationTableData, COMPONENT_Cr);
	}
	string fileName = image_to_export_filename;
	ofstream output(fileName.c_str(), ios::out | ios::binary);
	// i am generating the cb and cr as size of Y which may not be ok so check it later 
	//for (int comp = 0; comp < numberofComponent; comp++) {
	//
	//	for (unsigned int i = 0; i < comp_height; i++)
	//	{
	//		for (unsigned int j = 0; j < comp_width; j++)
	//		{
	//			switch (comp)
	//			{
	//			case COMPONENT_Y:
	//				output << Y[i][j];
	//				break; 
	//			case COMPONENT_Cb:
	//				output << CB[i][j];
	//				break;
	//			case COMPONENT_Cr:
	//				output << CR[i][j];
	//				break; 
	//			}
	//		}
	//	}

	//}
	
// writing just Y only 	not padded
	for (unsigned int i = 0; i < height; i++)
		{
			for (unsigned int j = 0; j < width; j++)
			{
					output << Y[i][j];
			}
		}

	// calcualate psnr 
	
	cv::Mat Y_org(height, width, CV_8UC1);
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			Y_org.at<uchar>(i, j) = luminance[i][j];
	//imshow("org", Y_org);
	//waitKey(0);
	cv::Mat Y_rec(height, width, CV_8UC1);
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j)
			Y_rec.at<uchar>(i, j) = Y[i][j];
	//imshow("rec ", Y_rec); 
	//waitKey(0);
	Mat s1;
	double y_psnr = getPSNR(Y_rec , Y_org , height , width );


	return  y_psnr;

}

void::jpeg_encoder::write_buff( vector<vector<uint_8>> & buffer, vector<int>& zigZagArray, int quantizationTable[8][8], int comp) {

	jpeg_decoder obj;
	int block_num = 0;  
	int data[64] = { 0 };
	int HFactor = jpegDecoder->components[comp].HFactor, VFactor = jpegDecoder->components[comp].VFactor;
	int HScale = jpegDecoder->components[comp].HScale, VScale = jpegDecoder->components[comp].VScale;
	int comp_height = ceil(1.0 * image_to_export_height_dct / VScale);
	int comp_width = ceil(1.0 * image_to_export_width_dct / HScale);
	double previousDCCoefficient = 0; //In this function DCPM is performed so DC coeficient is generated with formula Diff=DCi-DCi-1
	double prev_previousDC = 0;// previous previous


	// Initialize the positions
	int currentX = 0;
	int currentY = 0;
	int currentBlockHFactor = 0;
	int currentBlockVFactor = 0;
	int width = jpegDecoder->jpegImageWidth, height = jpegDecoder->jpegImageHeight;

	for (uint y_count = 0; y_count < comp_height; y_count += 8) { //comp_height
		for (uint x_count = 0; x_count < comp_width; x_count += 8) { //comp_width

			int data[64] = { 0 };

			for (int k = 0; k < 64; k++)
			{
				data[k] = zigZagArray.at(block_num * 64 + k);
			}
			int block8x8[8][8] = { 0 };
			int dataReshapedInto8x8[8][8] = { 0 };
			obj.reshapeArrayInto8x8(dataReshapedInto8x8, data);

			/*std::cout << "data reshape to 8 x 8 " << std::endl;
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
					std::cout << dataReshapedInto8x8[i][j] << " " << std::flush;
				std::cout << std::endl;
			}
*/
			obj.inverseZigZagScanning(block8x8, dataReshapedInto8x8);

			/*std::cout << "data reshape to 8 x 8  dezigzag " << std::endl;
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
					std::cout << block8x8[i][j] << " " << std::flush;
				std::cout << std::endl;
			}*/

			for (int u = 0; u < 8; ++u) {
				for (int v = 0; v < 8; ++v) {
					block8x8[u][v] = block8x8[u][v] * quantizationTable[u][v];
					if (!u && !v) {

						// Store the DC coefficient
						block8x8[u][v] = (block8x8[u][v] + previousDCCoefficient);
						previousDCCoefficient = block8x8[u][v]; // set the previous DC
					}
				}
			}
			/*	std::cout << "dequantized data and DPCM" << std::endl;
				for (int i = 0; i < 8; i++)
				{
					for (int j = 0; j < 8; j++)
						std::cout << block8x8[i][j] << " " << std::flush;
					std::cout << std::endl;
				}
	*/
	////  third do idct  
			int output_value_dct[8][8] = { 0 };
			jpegDecoder->perform_idct(output_value_dct, block8x8);
			/*
						std::cout << "IDCT" << std::endl;
						for (int i = 0; i < 8; i++)
						{
							for (int j = 0; j < 8; j++)
								std::cout << output_value_dct[i][j] << " " << std::flush;
							std::cout << std::endl;
						}
			*/
			for (int y = 0; y < 8; ++y) {
				for (int x = 0; x < 8; ++x)
				{
					// TODO: use switch instead of if
					// Only 8 and 12-bit is supported by DCT per JPEG standard
					if (jpegDecoder->jpegImageSamplePrecision == 8) {
						output_value_dct[x][y] = Clamp(output_value_dct[x][y] + 128);
					}
					else if (jpegDecoder->jpegImageSamplePrecision == 12) {
						output_value_dct[x][y] = Clamp(output_value_dct[x][y] + 2048);
					}
				} // end inner loop

			} // end outer loop

			for (int y = 0; y < 8; ++y) {

				if (currentX >= width) break;

				if (currentY + y >= height) break; 

				// Repeat each line VScale times
				for (int vfy = 0; vfy < VScale; ++vfy) {

					int picture_y = currentY + y * VScale + vfy;
					if (picture_y >= height) break;

					for (int x = 0; x < 8; ++x) {

						if (currentX + x * HScale >= width) break;

						// x and y are reverted compared to image
						int value = output_value_dct[y][x];
						
						// Repeat each pixel HScale times
						int realx = currentX + x * HScale;
						for (int i = 0; i < HScale; ++i) {

							// Set the pixel <Xcor: realX, Ycor: picture_y>

							buffer[picture_y][realx] = value;
							//std::cout << value << std::endl;
							if (++realx >= width) break;
						} // end for i
					} // end for x
				} // end for vfy
			} // end for y
			
			currentX += 8 * HScale;
			currentBlockHFactor++;
			if (currentBlockHFactor >= HFactor) {

				// restore the current X to its initial position and reset the counters
				currentX -= 8 * HScale * HFactor;
				currentBlockHFactor = 0;

				// go to next line
				currentY += 8 * VScale;
				currentBlockVFactor++;

				// you made a column of blocks
				if (currentBlockVFactor>= VFactor) {

					// restore the current Y to its initial position and reset the counters
					currentY -= 8 * VScale * VFactor;
					currentBlockVFactor = 0;

					currentX += 8 * HScale * HFactor;
					if (currentX >= width) {
						currentX = 0;
						currentY += 8 * VScale * VFactor;

					} // end if (currentX >= width)

				} // end  if (currentBlockVFactor[comp] >= VFactor)
			} // end if(currentBlockHFactor[comp] >= HFactor)


			block_num++;

		}// inner loop 
	}//outer loop 

}


float jpeg_encoder::C(int u)
{
	if (u == 0)
		return (1.0f / sqrtf(2));
	else
		return 1.0f;
} // end C

int jpeg_encoder::func(int x, int y, const int block[8][8])
{
	double sum = 0;
	for (int u = 0; u < 8; ++u)
	{
		for (int v = 0; v < 8; ++v)
		{
			sum += (C(u) * C(v)) * block[u][v] * cosine_idct[u][x] * cosine_idct[v][y];
		} // end inner loop
	} // end outer loop
	return (int)((1.0 / 4.0) * sum);
} // end func

void jpeg_encoder::perform_idct(int outBlock[8][8], const int inBlock[8][8])
{
	for (int y = 0; y < 8; ++y)
	{
		for (int x = 0; x < 8; ++x)
		{
			outBlock[x][y] = func(x, y, inBlock);
		} // end inner loop
	} // end outer loop
} // end perform_idct



Scalar jpeg_encoder::getMSSIM(const Mat& i1, const Mat& i2)
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


 double jpeg_encoder::getPSNR(const cv::Mat& i1, const cv::Mat& i2, int comp_width, int comp_height)
{

	int d = CV_32F;
	Mat I1, I2;
	i1.convertTo(I1, d);           // cannot calculate on one byte large values
	i2.convertTo(I2, d);


	Mat s1;
	absdiff(I1, I2, s1);       // |I1 - I2|
	s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
	s1 = s1.mul(s1);           // |I1 - I2|^2

	Scalar s = sum(s1);         // sum elements per channel
	double mse = (s[0]); // sum channels

	// cout << "MSE: " << mse << endl;
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

cv::Mat jpeg_encoder::vec2mat(vector< vector<double> >& matrix) {

	cv::Mat matrix_CV(0, matrix.size(), cv::DataType<double>::type);

	for (unsigned int i = 0; i < matrix.size(); ++i)
	{
		// Make a temporary cv::Mat row and add to NewSamples _without_ data copy
		cv::Mat Sample(1, matrix[0].size(), cv::DataType<double>::type, matrix[i].data());
		matrix_CV.push_back(Sample);
	}
	return matrix_CV;
}