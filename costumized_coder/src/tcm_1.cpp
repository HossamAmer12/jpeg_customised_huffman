// #include "stdafx.h"
//Mat img1 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/cheetahGray.png", 0);
//Mat img1 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/gooseGray.png", 0);
//
//
//namedWindow( "Display window", CV_WINDOW_AUTOSIZE );// Create a window for display.
//
//
//if(img1.empty())
//{
//    fprintf(stderr, "failed to load input image\n");
//    return -1;
//}
//
//int width  = img1.cols;
//int height = img1.rows;
//int Len = width*height;
//
////perform TCM and get outlier Yc
//int peak;
//double prob, lambda;
//double Yc;
//
//
//short * data = img1.ptr<short>(); // ptr to 1st element
//TCMprocessOneSequence(data, Len, &peak, &prob, &lambda, &Yc);
//
//cout << "Image 1 " << width << ", " << height << endl;
//cout << "Yc Gray is: " << Yc << endl;
////
//
//// --- Second image
//
////    Mat img2 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/cheetahTCM.png", 0);
//Mat img2 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/gooseTCM.png", 0);
//int width2  = img2.cols;
//int height2 = img2.rows;
//int Len2 = width2*height2;
//
//
//int peak2;
//double prob2, lambda2;
//double Yc2;
//
//if(img2.empty())
//{
//    fprintf(stderr, "failed to load input image\n");
//    return -1;
//}
//
//short * dataTCM = img2.ptr<short>(); // ptr to 1st element
//TCMprocessOneSequence(dataTCM, Len2, &peak2, &prob2, &lambda2, &Yc2);
//
//cout << "Image 2 " << width2 << ", " << height2 << endl;
//cout << "Yc TCM is: " << Yc2 << endl;
//
////    imshow( "Gray" , img1 );
//imshow( "TCM" , img2 );
//waitKey(0);


// August 9, 2018 -- TCM Main
//    Mat img = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/gooseGray.png", 0);
//Mat img = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/cheetahGray.png", 0);
//
//if(img.empty())
//{
//    fprintf(stderr, "failed to load input image\n");
//    return -1;
//}
//
//// Copy
//size_t len1 = img.total() * img.elemSize(); // in bytes
//short *data = new short[len1]; // *not* float!
//memcpy(data, img.ptr<short>(), len1);
//
//int width  = img.cols;
//int height = img.rows;
//int Len = width*height;
//
////perform TCM and get outlier Yc
//int peak;
//double prob, lambda;
//double Yc;
//TCMprocessOneSequence(data, Len, &peak, &prob, &lambda, &Yc);
//
//cout << "Yc Gray Image is: " << Yc << endl;
//
///// -------------------------
//
////    Mat img2 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/gooseTCM.png", 0);
//Mat img2 = imread("/Users/hossam.amer/7aS7aS_Works/work/my_Tools/tcmCpp/TCM/TCM/cheetahTCM.png", 0);
//width  = img2.cols;
//height = img2.rows;
//Len = width*height;
//
//// Copy
//size_t len = img2.total() * img2.elemSize(); // in bytes
//short *cs = new short[len]; // *not* float!
//memcpy(cs, img2.ptr<short>(), len);
//
//
//double Yc2;
//TCMprocessOneSequence(cs, Len, &peak, &prob, &lambda, &Yc2);
//
//
//double difference = abs((Yc2-Yc));
//cout << "Yc TCM Image is: " << Yc2 << endl;
//cout << "Difference is: " << difference << endl;
//
//ostringstream convert;   // stream used for the conversion
//ostringstream convert2;   // stream used for the conversion
//convert << Yc << ", Distance = " << difference;      // insert the textual representation of 'Number' in the characters in the stream
//convert2 << Yc2 << ", Distance = " << difference;      // insert the textual representation of 'Nu
//
//
//string s1 = "Gray Scale Image - Yc = "  + convert.str();
//string s2 = "TCM Outlier Image - Yc = " + convert2.str();
//imshow( s1 , img );
//imshow( s2 , img2 );
//waitKey(0);
//
//// free up memory
//delete [] data;
//delete [] cs;
//

//----//        for (uint j=0; j< table->codes.size(); j++) {
//
//            // If code is not used
//            if (table->codes[j] == 0xFFFFFFFF) {
//                continue;
//            }
//
//            // Represent the code in 8 bits
//            std::bitset<16> foo((int)table->codes[j]);
//
////            cout << "Codes with Length " << (table->codeLengths[j]) << " " << foo << endl;
//            cout << "Code no "<< j << " code " << foo  << " Length " << (table->codeLengths[j]) << endl;
//        }
//

//        for(int i = 0; i < 16; ++i)
//        {
//            for(int j = 0; j < codeLengths[i]; ++j)
//            {
//                if (table->codes[j] == 0xFFFFFFFF)
//                    continue;
//
//                std::bitset<8> foo((int)table->codes[j]);
//
//                if(j == 0) {
//                     cout << "Codes with Length " << (table->codeLengths[j]) << " are "  << " " << foo;
//                }
//                else {
//                    cout << " " << foo;
//                }
//
//
//            }
//
//            cout << "\n";
//        }

// ------------------------------

// August 10: Huffman print code
//        for (uint j=0; j< table->codes.size(); j++) {
//
//            // If code is not used
//            if (table->codes[j] == 0xFFFFFFFF) {
//                continue;
//            }
//
//            // Represent the code in 8 bits
//            std::bitset<16> foo((int)table->codes[j]);
//
////            cout << "Codes with Length " << (table->codeLengths[j]) << " " << foo << endl;
//            cout << "Code no "<< j << " code " << foo  << " Length " << (table->codeLengths[j]) << endl;
//        }
//

//        for(int i = 0; i < 16; ++i)
//        {
//            for(int j = 0; j < codeLengths[i]; ++j)
//            {
//                if (table->codes[j] == 0xFFFFFFFF)
//                    continue;
//
//                std::bitset<8> foo((int)table->codes[j]);
//
//                if(j == 0) {
//                     cout << "Codes with Length " << (table->codeLengths[j]) << " are "  << " " << foo;
//                }
//                else {
//                    cout << " " << foo;
//                }
//
//
//            }
//
//            cout << "\n";
//        }

// Main loop
//do {
//    // 3 bits is too small for a code
//    if (currentDataLength<3) continue;
//    
//    // Some stats
//    byteno++;
//#if DEBUGLEVEL>1
//    cout << "Byte "<<byteno<<" CDL "<<currentDataLength<<" value "<<hex<<data<<dec<<endl;
//#elif DEBUGLEVEL>0
//    if (byteno % 1000 == 0)
//    cout << "Byte "<<byteno<<" CDL "<<currentDataLength<<" value "<<hex<<data<<dec<<endl;
//#endif
//    
//    // Current Huffman table
//    HuffmanTable* htable = componentTablesDC[currentComponent];
//    if (ACDC == AC) htable = componentTablesAC[currentComponent];
//    
//    // Every one of 256 elements of the current Huffman table potentially has value, so we must go through all of them
//    for (int i=0; i<256; i++) {
//        // If code for i-th element is -1, then there is no Huffman code for i-th element
//        if (htable->codes[i] == 0xFFFFFFFF)
//            continue;
//        
//        // If current data length is greater or equal than n, compare first n bits (n - length of current Huffman code)
//        uint n = htable->codeLengths[i];
//        
//        if (currentDataLength >= n && htable->codes[i] == data >> (currentDataLength-n)) {
//#if DEBUGLEVEL>1
//            cout << "Found data "<<hex<<htable->codes[i]<<dec<<" len "<<n<<" at index "<<i<<endl;
//#endif
//            
//            // Remove first n bits from data;
//            currentDataLength -= n;
//            data = data - (htable->codes[i] << currentDataLength);
//            
//            // Reading of DC coefficients
//            if (ACDC == DC) {
//                unsigned char bitLength = i; // Next i bits represent DC coefficient value
//                
//                // Do we need to read more bits of data?
//                while (currentDataLength<bitLength) {
//                    if (!readMoreData(picture, data, currentDataLength)) {
//                        endOfFile = true;
//#if DEBUGLEVEL>0
//                        qDebug() << "End of file encountered inside a Huffman code!";
//#endif
//                        break;
//                    }
//                    byteno++;
//#if DEBUGLEVEL>1
//                    cout << "Byte(+) "<<byteno<<" CDL "<<currentDataLength<<" value "<<hex<<data<<dec<<endl;
//#endif
//                }
//                
//                // Read out DC coefficient
//                int DCCoeficient = data >> (currentDataLength-bitLength);
//                currentDataLength -= bitLength;
//                data = data - (DCCoeficient<<currentDataLength);
//                
//                // If MSB in DC coefficient starts with 0, then substract value of DC with 2^bitlength+1
//                //cout << "Before substract "<<DCCoeficient<<" BL "<<int(bitLength)<<endl;
//                if ( bitLength != 0 && (DCCoeficient>>(bitLength-1)) == 0 ) {
//                    DCCoeficient = DCCoeficient - (2 << (bitLength-1)) + 1;
//                }
//                //cout << "After substract "<<DCCoeficient<<" previousDC "<<previousDC[currentComponent]<<endl;
//                
//                dataBlock[m] = DCCoeficient + previousDC[currentComponent];
//#if DEBUGLEVEL>1
//                cout << "DC READ "<<dataBlock[m]<<" at index "<<m<<endl;
//#endif
//                m++;
//                
//                // No AC coefficients required?
//                if (ACcount == 0 || losslessFormat) return;
//                
//                // We generated our DC coefficient, next one is AC coefficient
//                ACDC = AC;
//                if (currentDataLength < 3) // If currentData length is < than 3, we need to read new byte, so leave this for loop
//                    break;
//                i =- 1; // CurrentDataLength is not zero, set i=0 to start from first element of array
//                htable = componentTablesAC[currentComponent];
//            }
//            
//            // Reading of AC coefficients
//            else {
//                unsigned char ACElement=i;
//                
//                /* Every AC component is composite of 4 bits (RRRRSSSS). R bits tells us relative position of
//                 non zero element from the previous non zero element (number of zeros between two non zero elements)
//                 SSSS bits tels us magnitude range of AC element
//                 Two special values:
//                 00 is END OF BLOCK (all AC elements are zeros)
//                 F0 is 16 zeroes */
//                
//                if (ACElement == 0x00)
//                    return;
//                
//                else if (ACElement == 0xF0) {
//                    for (int k=0;k<16;k++) {
//                        dataBlock[m] = 0;
//                        m++;
//                        if (m >= ACcount+1) {
//#if DEBUGLEVEL>0
//                            qDebug() << "Huffman error: 16 AC zeros requested, but only "<<k<<" left in block!";
//#endif
//                            return;
//                        }
//                    }
//                }
//                else {
//                    /* If AC element is 0xAB for example, then we have to separate it in two nibbles
//                     First nible is RRRR bits, second are SSSS bits
//                     RRRR bits told us how many zero elements are before this element
//                     SSSS bits told us how many binary digits our AC element has (if 1001 then we have to read next 9 elements from file) */
//                    
//                    // Let's separate byte to two nibles
//                    unsigned char Rbits = ACElement >> 4;
//                    unsigned char Sbits = ACElement & 0x0F;
//                    
//                    // Before our element there is Rbits zero elements
//                    for (int k=0; k<Rbits; k++) {
//                        if (m >= ACcount) {
//#if DEBUGLEVEL>0
//                            qDebug() << "Huffman error: "<<Rbits<<" preceeding AC zeros requested, but only "<<k<<" left in block!";
//#endif
//                            // in case of error, doing the other stuff will just do more errors so return here
//                            return;
//                        }
//                        dataBlock[m] = 0;
//                        m++;
//                    }
//                    
//                    // Do we need to read more bits of data?
//                    while (currentDataLength<Sbits) {
//                        if (!readMoreData(picture, data, currentDataLength)) {
//                            endOfFile = true;
//#if DEBUGLEVEL>0
//                            qDebug() << "End of file encountered inside a Huffman code!";
//#endif
//                            break;
//                        }
//                        byteno++;
//#if DEBUGLEVEL>1
//                        cout << "Byte(+) "<<byteno<<" CDL "<<currentDataLength<<" value "<<hex<<data<<dec<<" Sbits "<<int(Sbits)<<endl;
//#endif
//                    }
//                    
//                    // Read out AC coefficient
//                    int ACCoeficient = data >> (currentDataLength-Sbits);
//                    currentDataLength -= Sbits;
//                    data = data - (ACCoeficient<<currentDataLength);
//                    
//                    // If MSB in AC coefficient starts with 0, then substract value of AC with 2^bitLength+1
//                    if ( Sbits != 0 && (ACCoeficient>>(Sbits-1)) == 0 ) {
//                        ACCoeficient = ACCoeficient - (2 << (Sbits-1)) + 1;
//                    }
//                    dataBlock[m] = ACCoeficient;
//                    m++;
//                }
//                
//                // End of block
//                if (m >= ACcount+1)
//                    return;
//                
//                if (currentDataLength<3) // If currentData length is < 3, we need to read new byte, so leave this for loop
//                    break;
//                i =- 1; // currentDataLength is not zero, set i=0 to start from first element of array
//            }
//            
//        }
//    }
//    } while(readMoreData(picture, data, currentDataLength));



// not working is_exist_in_huffman_codes
//bool jpeg_decoder::is_exist_in_huffman_codes(int code, int codeLength, int currentComponent, int& decodedValue, bool is_dc) {
//
//    // pointer to a huffman table
//    HuffmanTable* hTable;
//
//    // is DC coefficient?
//    if(is_dc) {
//        hTable = componentTablesDC[currentComponent];
//    }
//    else {
//        hTable = componentTablesAC[currentComponent];
//    }
//
//    // Search for the input in your hashmap
////    auto search = hTable->huffData.find(huffKey(codeLength, code));
////    auto search = hTable->huffData.find(huffKey(code, codeLength));
//
//    std::map<huffKey, uint_8>::iterator iter;
//    for (iter = hTable->huffData.begin();iter != hTable->huffData.end(); ++iter) {
//
//        int codeInTable = iter->first.second;
//        int codeLengthInTable = iter->first.first;
//
//
//        // Print Code - Its Length : Equivalent letter
//        printf("    %04X at length %d = %02X or [%s] \n",
//               iter->first.second, iter->first.first,
//               iter->second, IntToBinary(iter->first.second, iter->first.first));
//    }
//
//    // not found
//    if(search != hTable->huffData.end())
//    {
//#if DEBUGLEVEL>2
//        cout << "Not found this code: " << code << ", with codeLength: " << codeLength << endl;
//#endif
//        return false;
//    }
//
//#if DEBUGLEVEL>2
//    printf("(Yaay) Found this code: %d, with codeLength: %d, and value = %02X\n",
//           code, codeLength, search->second);
//#endif
//
//    // set the value of decodedValue from huffman table
//    decodedValue = search->second;
//    return true;
//
//}

