//
//  exportpicture.hpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-16.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//

#ifndef EXPORT_PICTURE_H
#define EXPORT_PICTURE_H

#include <stdio.h>
#include <vector>
#include "jpegdecoder.h"
#include "quantizationtable.h"


class ExportPicture {
public:    
    ExportPicture ();
    ExportPicture (const jpeg_decoder& process_jpeg_decoder);
    
    vector<vector<int> >luminance;
    vector<vector<int> >chrominanceCb;
    vector<vector<int> >chrominanceCr;
    int width;
    int height;

    // TODO: need to copy the componentQuantizationTables
    
    void copy_data_to_vectors(const jpeg_decoder& processed_jpeg_decoder);
    void delete_data_from_vectors();
    
};


#endif /* EXPORT_PICTURE_H */
