//
//  main.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-06-27.
//  Copyright © 2018 Hossam Amer. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <math.h>
//#include <ctime>
//#include <cv.h>
//#include <cxcore.h>
//#include <highgui.h>
//#include <opencv2/opencv.hpp>
//#include <opencv2/highgui/highgui.hpp>


#include "tcm.h"
#include "jpegdecoder.h"

#include "jpegencoder.h"


using namespace std;
//using namespace cv;


#define IS_MAIN_NEW    2


#include "jpegencoderMultipleQF.h"


#if IS_MAIN_NEW

#if IS_MAIN_NEW>1
int main(int argc, const char * argv[]) {
    // insert code here...
    
    
    if(argc < 4)
    {
        // Tell the user how to run the program
        std::cerr << "Number of arguments should be 3: <input_file_with_full_path> <output_folder> <output_txt_folder>" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }
    
    /// JPEG Stuff
    //    std::string path_to_files = "/Users/hossam.amer/7aS7aS_Works/work/my_Tools/jpeg_tcm/dataset/";
    
    // Input file:
    std::string filename = argv[1];
    
    // Ouptut folder path:
    std::string enc_path_to_files = argv[2];
    
    // Ouptut TXT folder path:
    std::string enc_path_txt = argv[3];
    // Input quality factor:
    //    std::string arg = argv[4];
    //    int quality_factor;
    //    try {
    //        std::size_t pos;
    //        quality_factor = std::stoi(arg, &pos);
    //        if (pos < arg.size()) {
    //            std::cerr << "Trailing characters after number: " << arg << '\n';
    //        }
    //    } catch (std::invalid_argument const &ex) {
    //        std::cerr << "Invalid number: " << arg << '\n';
    //        return 1;
    //    } catch (std::out_of_range const &ex) {
    //        std::cerr << "Number out of range: " << arg << '\n';
    //        return 1;
    //    }
    //
    // Quality factor experiment:
    ////////////////////////////////////////////////////////
    
    try {
        // Hossam: Save the input fileName
        std::string encoded_filename = filename;
        
        ////// String Processing -- Get the file Name
        size_t found = encoded_filename.find_last_of("/\\");
        std::string filename_first_token = encoded_filename.substr(found+1);
        found = filename_first_token.find_first_of(".");
        std::string filename_second_token = filename_first_token.substr(0, found);

        runEncoderWithMultipleQF(filename, enc_path_to_files, enc_path_txt);
        
        
        //    } catch (exception e) {
    } catch (const std::exception &e) {
        cerr << "Input the folder properly" << endl;
        std::cerr << "My error is: " << e.what() << endl;
        return 1;
        
    }

    return 0;
} // end main


#else
int main(int argc, const char * argv[])
{
    
    std::string path_to_files = "/Users/ahamsala/Documents/validation_original/shard-0/1/";
    std::string enc_path_to_files = "/Users/ahamsala/Documents/validation_generated_QF_1/";


    //    std::string filename = path_to_files + "ILSVRC2012_val_00000034.JPEG";
    std::string filename = path_to_files + "ILSVRC2012_val_00000004.JPEG";
    //    std::string filename = path_to_files + "ILSVRC2012_val_00000126.JPEG";
    
    jpeg_decoder test(filename);
    
    // Quality factor experiment:
    ////////////////////////////////////////////////////////
    // Hossam: Save the input fileName
    std::string encoded_filename = filename;
    
    ////// String Processing -- Get the file Name
    size_t found = encoded_filename.find_last_of("/\\");
    std::string filename_first_token = encoded_filename.substr(found+1);
    found = filename_first_token.find_first_of(".");
    std::string filename_second_token = filename_first_token.substr(0, found);
    

    //    std::string enc_path_to_files = "/Users/ahamsala/Documents/4.DataBase_1/JPEG/";
    // for each quality factor
    int end_quality_factor = 100;
    std::string enc_path_txt="/Users/ahamsala/Documents/validation_generated_QF_TXT_1/";
    runEncoderWithMultipleQF(filename, enc_path_to_files, enc_path_txt);
    
//        for (int quality_factor = 0; quality_factor <= end_quality_factor ; quality_factor += 5)
//    //    for (int quality_factor = QFACTOR; quality_factor <= 10 ; quality_factor += 10)
//        {
//            encoded_filename = enc_path_to_files + filename_second_token + "-QF-" + to_string(quality_factor) + filename_first_token.substr(found);
//            cout << "encoded Output: " << encoded_filename << endl;
//            // Encoded output name
//            jpeg_encoder enc(&test, encoded_filename, quality_factor);
//            enc.savePicture();
//            ///
//            ofstream myfile_txt;
//            cout<<enc_path_txt + filename_second_token + "-QF-" + to_string(quality_factor) +".txt"<<endl;
//            myfile_txt.open(enc_path_txt + filename_second_token + "-QF-" + to_string(quality_factor) +".txt");
//            myfile_txt << to_string(enc.image_bpp)+"\n";
//            myfile_txt.close();
//            ///
//        }
 
    return 0;
}
#endif
#endif
