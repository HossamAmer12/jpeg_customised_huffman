//
//  jpegencoderMultipleQF.h
//  jpeg_tcm
//
//  Created by Hossam Amer on 2018-11-04.
//  Copyright © 2018 Hossam Amer. All rights reserved.
//

#ifndef JPEG_ENCODER_MULTIPLE_QF_H
#define JPEG_ENCODER_MULTIPLE_QF_H

#include <thread>
#include <mutex>
#include <functional>   // std::bind

#include "jpegdecoder.h"
#include "jpegencoder.h"

// filename: input file name
// enc_path_to_files: path to store the encoded pictures
static void runEncoderWithMultipleQF(std::string filename, std::string enc_path_to_files, std::string enc_path_txt)
{
	cout << "hii" << std::endl;
	// Decode:
	cout << "Start Decode " << filename << endl;
	jpeg_decoder test(filename);
	//    cout<<filename<<endl;
	cout << "Done Decode " << filename << endl;


	// Number of loops is 11 (11 quality factors)
	int nloop = 21;


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
					// Encode:
					// inner loop
					{
						const int quality_factor = i * 5;
						// Hossam: Save the input fileName
						std::string encoded_filename = filename;

						////// String Processing -- Get the file Name
						size_t found = encoded_filename.find_last_of("/\\");
						std::string filename_first_token = encoded_filename.substr(found + 1);
						found = filename_first_token.find_first_of(".");
						std::string filename_second_token = filename_first_token.substr(0, found);
					
						// Update the full path for the encoded_file name
						encoded_filename = enc_path_to_files + "\\" + filename_second_token + "-QF-" + to_string(quality_factor) + ".yuv"; //filename_first_token.substr(found);

						//                                                           cout << "\nStart Encode " << encoded_filename << " @ " << quality_factor << endl;
						jpeg_encoder enc(&test, encoded_filename, quality_factor);
						//enc.savePicture();
						double psnr = enc.writeyuv();
						
						//ofstream myfile_txt;

						/*myfile_txt.open(enc_path_txt + "\\" + filename_second_token + "-QF-" + to_string(quality_factor) + ".txt");
						cout << "bpp is" <<to_string(enc.image_bpp) + "\n";
						cout << to_string(psnr) + "\n";
						myfile_txt << to_string(enc.image_bpp) + "\n";
						myfile_txt << to_string(psnr) + "\n";
						myfile_txt.close();
						cout << "Done Encode; Output is " << encoded_filename << endl;
						*/
						//cout << "Txt Gerneration " << enc_path_txt + "\\" + filename_second_token + "-QF-" + to_string(quality_factor) + ".txt" << endl;
						// (optional) make output critical
//                                                           std::lock_guard<std::mutex> lock(critical);
//                                                           std::cout << bi << " " << ei << " " << quality_factor <<std::endl;
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

		// Post loop:
		// ---------

		std::string encoded_filename = filename;

		////// String Processing -- Get the file Name
		size_t found = encoded_filename.find_last_of("/\\");
		std::string filename_first_token = encoded_filename.substr(found + 1);
		found = filename_first_token.find_first_of(".");
		std::string filename_second_token = filename_first_token.substr(0, found);
		cout << "\n\nDone Execution; Output is: " << enc_path_to_files + filename_second_token + "-QF-*" << endl;
	}

	return;
}


#endif /* JPEG_ENCODER_MULTIPLE_QF_H */