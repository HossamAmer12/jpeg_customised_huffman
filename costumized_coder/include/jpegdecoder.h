//
//  jpegdecoder.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef JPEG_DECODER_H
#define JPEG_DECODER_H

#include <stdio.h>

#include <math.h>

#include "inttypes.h"
#include "component.h"
#include "huffmantable.h"

// Include the tcm.h for the encoder TCM computations
#include "tcm.h"

// for debugging LookNBits
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#if OPEN_CV_ENABLED
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#endif


// maximum number of color components e.g. 3 for Y'CbCr, 4 for CMYK
#define ETF_FORMAT_MAX_COMPONENTS  4 // 4 components for Adobe format

// Macro to read a 16-bit word from file
#define READ_WORD() ((fgetc(fp) << 8) | fgetc(fp))

// Segment parsing error codes
#define JPEG_SEG_ERR  0
#define JPEG_SEG_OK   1
#define JPEG_SEG_EOF -1

class jpeg_decoder {
    
public:
    
    // Constructor
    jpeg_decoder(std::string);
    
    // Destructor
    ~jpeg_decoder();
    
    // Y buffer for the entire picture:
    unsigned char** m_YPicture_buffer;
    
    // Cb buffer for the entire picture:
    unsigned char** m_CbPicture_buffer;
    
    // Cr buffer for the entire picture:
    unsigned char** m_CrPicture_buffer;
    
    // Vector of DCT blocks: (Y: the block index, X: dct coefficients)
    vector<vector <int> > tCoeff_Y;
    vector<vector <int> > tCoeff_Cb;
    vector<vector <int> > tCoeff_Cr;
    
    // Vector of the DCT coefficient categorized by the AC indexes
    vector<vector <int> > tCoeff_Y_AC;
    
    // convert JPG to BMP file; must be called after your create the jpeg decoder
    int convert_jpg_to_bmp(const std::string output_bmp_filename);
    
    
#if OPEN_CV_ENABLED
    // display Y channel via OpenCV (comp should be unsigned int > 0)
    int display_jpg_yuv(std::string windowName, uint_32 comp = COMPONENT_Y) const;
    
    // display RGB via OpenCV
    int display_jpg_bmp(const std::string windowName) const;
    
#endif
    
    // save Y channel to a CSV file; must be called after your create the jpeg decoder
    int write_yuv_from_jpg_in_csv(const std::string output_csv_name);
    
    // Saves the Y channel in YUV file format (can do YUV later)
    int write_yuv_from_jpg_in_yuv(const std::string output_yuv_filename);
    
    
    // write Y channel dct to a CSV file; must be called after your create the jpeg decoder
    int write_tcoeff_y_from_jpg_in_csv(const std::string output_csv_name, int comp = COMPONENT_Y);
    
    
    // getters:
    uint_16 get_image_width() const { return jpegImageWidth; }
    uint_16 get_image_height() const { return jpegImageHeight; }
    
    
    //private:
    
    // TODO: make getters
public:
    
    // The file to be read from, opened by constructor
    FILE *fp;
    
    // the input file name
    std::string jpeg_filename;
    
    // the output file name
    std::string jpeg_output_bmp_filename;
    
    // Names of the possible segments
    std::string segNames[64];
    
    // Temp space used after the IDCT to store each component
    // should not be declared as automatic storge - you will need it
    // to the end of the program
    unsigned char*		m_Y; // at most 64 coefficients and 4 Y components (256 elements; max 2x2 hFactor*vFactor, each block is 8x8 = 8x8x2x2)
    unsigned char*		m_Cr;  // at most 64 coefficients and 1 Cb component
    unsigned char*		m_Cb;  // at most 64 coefficients for 1 Cr component
    
    // Final Red Green Blue pixel data
    unsigned char*		m_rgb;
    // Internal Pointer use for colorspace conversion, do not modify it !!! (it points to the m_rgb pointer)
    unsigned char*		m_colourspace; // August 14: TODO: to be removed
    
    
    unsigned char precision; // Precison of elements (8 or 12 bits) // not used
    
    // Data needed for decoding
    vector <Component> components;
    vector <QuantizationTable> quantizationTables;
    vector <HuffmanTable*> huffmanTables;
    HuffmanTable* componentTablesDC[ETF_FORMAT_MAX_COMPONENTS]; // from format.h
    HuffmanTable* componentTablesAC[ETF_FORMAT_MAX_COMPONENTS];
    bool endOfFile;
    int previousDC[ETF_FORMAT_MAX_COMPONENTS];
    
    bool progressive_Huff_Format = false; // is progressive_Huff_Format
    bool losslessFormat; // is lossless format or lossy
    int zigZagStart, zigZagEnd;
    uint_8 ssa;
    int approximationH, approximationL;
    int *scanLineCache[ETF_FORMAT_MAX_COMPONENTS];
    
    
    // Data used by addBlock method
    //    unsigned char* rawImagePointers[ETF_FORMAT_MAX_COMPONENTS];
    int lineBytes; // how many bytes in a line of an image (YCbCr YCbCr...etc)
    unsigned int currentX[ETF_FORMAT_MAX_COMPONENTS]; // Coordinates for next block
    unsigned int currentY[ETF_FORMAT_MAX_COMPONENTS];
    int maxSample; // maximum sample value (2^bit_precision - 1)
    bool hasSubSampling;
    int currentBlockHFactor[ETF_FORMAT_MAX_COMPONENTS], currentBlockVFactor[ETF_FORMAT_MAX_COMPONENTS];
    
    // indices of next block for DCT blocks
    unsigned int currentX_dct[ETF_FORMAT_MAX_COMPONENTS]; // Coordinates for next block
    unsigned int currentY_dct[ETF_FORMAT_MAX_COMPONENTS];
    int currentBlockHFactor_dct[ETF_FORMAT_MAX_COMPONENTS], currentBlockVFactor_dct[ETF_FORMAT_MAX_COMPONENTS];
    
    
    // Data used by IDCT method
    double cosine[106];
    int coefficients[4096];
    
    double cosine_idct[8][8];
    
    // information about the picture
    uint_16 jpegImageWidth, jpegImageHeight;
    int upscale_height;
    int upscale_width;
    int mcu_width, mcu_height;
    int mcu_cols, mcu_rows;
    
    // CHANGES FOR PROGRESSIVE MODE -------------------
    uint_8  numberOfComponents;
    vector<uint_8> componentID;
    int application_size = 0;
    vector<vector <int> > data_DCT;
    int EOB_run = 0;
    int counter_scan_blockidx = 0;
    void final_process_progressive();
    
    int total_block_Y, total_block_C;
    int count_block_Y, count_block_Cb, count_block_Cr;
    int jpegImageSamplePrecision; // 8 or 12 bits (frame header)
    
    
    // Main loop
    void readFile();
    
    // parse segment
    int parseSeg();
    
    uint_16 readFrameHeader(uint_16 headerLengthFromBitStream);
    uint_16 readHuffmanTables(uint_16 tableLengthFromBitStream);
    void readArithmeticCoding(); // TODO: incomplete
    uint_16 readQuantizationTables(uint_16 tableLengthFromBitStream);
    void readComments(); // TODO: incomplete
    uint_16 readScanHeader(uint_16 headerLengthFromBitStream);
    
    // createRawImagePointers
    void createRawImagePointers();
    void deleteRawImagePointers();
    void deleteDCTComponentPointers();
    void deleteRawPictureBufferPointers();
    
    // TODO: Support intervleaving order of YCbCr pixels
    void readImageEntryPoint();
    void readProgressiveImageEntryPoint();
    void decode_mcu(int componentWidth, int componentHeight, int currentX, int currentY);
    // FOR PROGRESSIVE
    void decode_mcu_progressive(int componentWidth, int componentHeight, int currentX, int currentY, vector<uint_8> componentID);
    
    void process_huffmann_data_unit_progressive(int currentComponent, int currentX, int currentY);
    void process_huffmann_data_unit(int currentComponent, int currentX, int currentY);
    // Decodes single 8x8 block
    void decode_single_block(int offset, int stride, int currentComponent, uint_8* outputBuf,
                             int currentX = 0, int currentY = 0);
    
    
    // initializes the necessary positions and buffers for storing the Y component
    void initPositionsBuffersForPictureBuffer();
    void initCurrentPosition();
    // adds single 8x8 block to the Y picture buffer
    void addBlockSubsampling(int dataBlock[8][8], int currentComponent = COMPONENT_Y);
    
    
    void addtCoeffBlock(int dataBlock[8][8], int currentComponent = COMPONENT_Y, int currentX = 0, int currentY = 0);
    
    void reshapeArrayInto8x8(int outArray[8][8], const int inArray[64]);
    
    bool is_exist_in_huffman_codes(int code, int codeLength, int currentComponent, int& decodedValue, bool is_dc = true);
    
    // YCbCr color conversion into RGB24 methods:
    // Formula: http://www.mir.com/DMG/ycbcr.html
    void convert_yccb_to_rgb24(int y, int cb, int cr, int& r, int& g, int& b);
    
    
    //    void ycrcb_to_rgb24_image(int w, int h, int imgx, int imgy, int imgw, int imgh);
    void ycrcb_to_rgb24_image();
    
    
    // Save a buffer in 24bits Bitmap (.bmp) format
    void write_bmp24(const std::string output_bmp_filename, int width, int height, unsigned char* RGB);
    
    
    
    // zig-zag scanning methods
    void inverseZigZagCoding(int *array, int matrix[8][8]); // not used
    
    void inverseZigZagScanning(int outMatrix[8][8], const int inMatrix[8][8]); // used
    
    bool isEndOfFile() { return endOfFile; }
    
    
    // image transpose
    void transpose(int outMatrix[8][8], const int inMatrix[8][8]);
    
    // TODO: Speed up IDCT with look up tables using the "cosines" instance variables
    // idct functions:
    void perform_idct(int outBlock[8][8], const int inBlock[8][8]);
    int func(int x, int y, const int block[8][8]);
    float C(int u);
    
    // decoding functions from the other source
    void IDCT (int block8x8[8][8], int block2[8][8]);
    void multiplyWithQuantizationTable(int dataBlock[8][8], int currentComponent);
    
    
    // void logErrorPictures() {
       
    //         cerr << "Error using jpegDecoder: Problem in parsing." << endl;
    //         size_t found = jpeg_filename.find_last_of("/\\");
    //         std::string filename_before_slash = jpeg_filename.substr(0, found);
        
    //         // up one level
    //         found = filename_before_slash.find_last_of("/\\");
    //         std::string filename_before_slash_up = filename_before_slash.substr(0, found);
        
    //         std::string filename_first_token = jpeg_filename.substr(found+1);
    //         found = filename_first_token.find_first_of(".");
    //         std::string filename_second_token = filename_first_token.substr(0, found);
    //         string log_filename = "";
    //         std::ostringstream oss;
    //         oss << filename_before_slash_up << "/Other_Error_Pictures.txt";
    //         log_filename = oss.str();
    //         char* log_filename_char = log_filename.empty()? NULL: strdup(log_filename.c_str());
    //         FILE* fileId = fopen (log_filename_char, "a");
    //         fprintf(fileId, "%s\n", filename_second_token.c_str());
    //         fclose(fileId);
        
    //         cout << "Exiting after making a log at: " << log_filename << endl;

    //         // Create a text file inside the folder with the image name
    //         oss.clear(); oss.str("");
    //         found = filename_second_token.find_last_of("/\\");
    //         filename_second_token = filename_second_token.substr(found+1);
    //         oss << filename_before_slash << "/" << filename_second_token << "_other" << ".txt";
    //         log_filename = oss.str();
    //         log_filename_char = log_filename.empty()? NULL: strdup(log_filename.c_str());
    //         fileId = fopen (log_filename_char, "w");
    //         fprintf(fileId, "%s\n", filename_second_token.c_str());
    //         fclose(fileId);
        
    //         cout << "Exiting after making a log at: " << log_filename << endl;

        
    //         exit(0);
    // }
    
    // void logCMYKErrorPictures() {
        
    //     cerr << "Error using jpegDecoder: JPEG images with CMYK colorspace are not currently supported." << endl;
    //     size_t found = jpeg_filename.find_last_of("/\\");
    //     std::string filename_before_slash = jpeg_filename.substr(0, found);
        
        
    //     // up one level
    //     found = filename_before_slash.find_last_of("/\\");
    //     std::string filename_before_slash_up = filename_before_slash.substr(0, found);
        
    //     std::string filename_first_token = jpeg_filename.substr(found+1);
    //     found = filename_first_token.find_first_of(".");
    //     std::string filename_second_token = filename_first_token.substr(0, found);
    //     string log_filename = "";
    //     std::ostringstream oss;
    //     oss << filename_before_slash_up << "/CMYK_Error_Pictures.txt";
    //     log_filename = oss.str();
    //     char* log_filename_char = log_filename.empty()? NULL: strdup(log_filename.c_str());
    //     FILE* fileId = fopen (log_filename_char, "a");
    //     fprintf(fileId, "%s\n", filename_second_token.c_str());
    //     fclose(fileId);
        
    //     cout << "Exiting after making a log at: " << log_filename << endl;
        
    //     // Create a text file inside the folder with the image name
    //     oss.clear(); oss.str("");
    //     found = filename_second_token.find_last_of("/\\");
    //     filename_second_token = filename_second_token.substr(found+1);
    //     oss << filename_before_slash << "/" << filename_second_token << "_cmyk" << ".txt";
    //     log_filename = oss.str();
    //     log_filename_char = log_filename.empty()? NULL: strdup(log_filename.c_str());
    //     fileId = fopen (log_filename_char, "w");
    //     fprintf(fileId, "%s\n", filename_second_token.c_str());
    //     fclose(fileId);
        
    //     cout << "Exiting after making a log at: " << log_filename << endl;

        
    //     exit(0);
    // }
    
    // Reading bits uitilities (declare them inline for speed purposes)
    
    // queue length
    unsigned int g_nbits_in_reservoir = 0;
    // the queue of bits
    unsigned int g_reservoir = 0;
    
    // changes: add int NumofByte in function parameter
    inline void FillNBits(int& nbits_wanted)
    {
        while ((int) g_nbits_in_reservoir < nbits_wanted)
        {
            const unsigned char c = fgetc(fp); // read a byte
            g_reservoir <<= 8; // left shift your reservoir to empty space for the new byte
            
            // if what I read is 0xFF and the next byte is 0x00
            // move the cusor one more byte --> SKIP!
            
            // Current position
            fpos_t current_file_position;
            // Get the current position
            fgetpos(fp, &current_file_position);
            
            const unsigned char c2 = fgetc(fp);
            
            if (c == 0xFF && c2 == 0x00)
            {
#if DEBUGLEVEL > 20
                cout << "0xFF00 Skip is done! " << endl;
#endif
            }
            else
            {
                // Restore the position
                fsetpos(fp, &current_file_position);
                
#if DEBUGLEVEL > 20
                cout << "Restore position for this byte! " << endl;
#endif
                
            }
            
            g_reservoir |= c; // read the new character (byte) in the left most part of your resevoir
            g_nbits_in_reservoir += 8; // increment your length by 1 byte
        } // end while loop
    } // end fillNBits
    
    inline short GetNBits(int nbits_wanted)
    {
        FillNBits(nbits_wanted);
        
        // right shift because your new byte is always stored at the left most of the reservoir
        short result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
        
        // decrement your reservoir length
        g_nbits_in_reservoir -= (nbits_wanted);
        
        // make sure that the left most bits only contain the bytes that are not yet fetched
        g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
        
        /*
         // Could do the sign conversion here!
         if (result < (short)(1UL<<((nbits_wanted)-1)))
         {
         result = result + (short)(0xFFFFFFFFUL<<(nbits_wanted))+1;
         }
         */
        return result;
    }
    
    
    inline void printNBits(int nbits_wanted)
    {
        FillNBits(nbits_wanted);
        
        int result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
        std::cout << "Read bits: " << std::hex << result  << endl;
    }
    
    inline int LookNBits(int nbits_wanted)
    {
        FillNBits(nbits_wanted);
        
        int result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
        return result;
    }
    
    inline void SkipNBits(int& nbits_wanted)
    {
        FillNBits(nbits_wanted);
        
        // decrement your reservoir length
        g_nbits_in_reservoir -= (nbits_wanted);
        
        // make sure that the left most bits only contain the bytes that are not yet fetched
        g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
    }
    
    
    // Clamp our integer between 0 and 255
    inline unsigned char Clamp(int i, int threshold = 255)
    {
        if (i < 0)
            return 0;
        
        // Note: maxSample should normally be 255
        // TODO: threshold should be maxSample (8 bits precision: 2^8-1 or 12 bits precision: 2^12-1)
        else if (i > threshold)
            return 255;
        else
            return i;
    }
    
    int DetermineSign(int val, int nBits)
    {
        bool negative = val < ( 1 << ( nBits - 1) );
        
        if (negative)
        {
            // (-1 << (s)), makes the last bit a 1, so we have 1000,0000 for example for 8 bits
            
            val = val + (-1 << (nBits)) + 1;
        }
        
        // Else its unsigned, just return
        return val;
    }
    
    char g_bigBuf[1024] = {0};
    char* IntToBinary(int val, int bits)
    {
        for (int i=0; i<32; i++) g_bigBuf[i]='\0';
        
        int c = 0;
        for (int i=bits-1; i>=0; i--)
        {
            bool on = (val & (1<<i)) ? 1 : 0;
            g_bigBuf[c] = on ? '1' : '0';
            c++;
        }
        
        return &g_bigBuf[0];
    }
    
    
    
    // debugging functions:
    void dumpBlock(int block[8][8]);
    void printDCHuffmanCodes(int currentComponent);
    void printACHuffmanCodes(int currentComponent);
    void dumpDCTValues(short dct[64]);
    void dumpDecodedBlock(int val[8][8]);
    void dumpDecodedBlockInDecimal(int val[8][8]);
    void dumpQuantizationTable(int currentComponent);
    
    inline static uint bit_count(int temp1)
    {
        if (temp1 < 0) {
            temp1 = -temp1;
        }
        
        uint nbits = 0;
        while (temp1) {
            nbits++;
            temp1 >>= 1;
        }
        return nbits;
    }
    
};


#endif /* jpegdecoder_h */
