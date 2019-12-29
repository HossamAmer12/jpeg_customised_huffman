//
//  exportpicture.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-16.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//
#if IS_RUN_MS_WINDOWS
#include "stdafx.h"
#endif

#include "exportpicture.h"


ExportPicture::ExportPicture () :luminance(0,vector<int>(0,0)),
                            chrominanceCb(0,vector<int>(0,0)),
                            chrominanceCr(0,vector<int>(0,0)),
                            width(0),
                            height(0)  {
    
} // end constructor

ExportPicture::ExportPicture(const jpeg_decoder& processed_jpeg_decoder): luminance(0,vector<int>(0,0)),
chrominanceCb(0,vector<int>(0,0)),
chrominanceCr(0,vector<int>(0,0)),
width(processed_jpeg_decoder.get_image_width()),
height(processed_jpeg_decoder.get_image_height()) {
    
    
    luminance.resize(height);
    chrominanceCb.resize(height);
    chrominanceCr.resize(height);
    for(int i = 0; i < height; ++i){
        luminance[i].resize(width);
        chrominanceCb[i].resize(width);
        chrominanceCr[i].resize(width);
    }
    
    copy_data_to_vectors(processed_jpeg_decoder);
} // end constructor


void ExportPicture::copy_data_to_vectors(const jpeg_decoder& processed_jpeg_decoder) {
    
    for(int i = 0; i < height; ++i)
    {
        for(int j = 0; j < width; ++j) {
            luminance[i][j]     = processed_jpeg_decoder.m_YPicture_buffer[i][j];
            chrominanceCb[i][j] = processed_jpeg_decoder.m_CbPicture_buffer[i][j];
            chrominanceCr[i][j] = processed_jpeg_decoder.m_CrPicture_buffer[i][j];
        }
    }
    
} // end copy_data_to_vectors

void ExportPicture::delete_data_from_vectors() {
    luminance.clear();
    chrominanceCb.clear();
    chrominanceCr.clear();
} // end delete_data_from_vectors

