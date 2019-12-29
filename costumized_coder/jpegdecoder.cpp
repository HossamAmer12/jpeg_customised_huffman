//
//  jpegdecoder.cpp
//  TCM
//
//  Created by Hossam Amer on 2018-08-09.
//  Author: Hossam Amer & Yanbing Jiang, University of Waterloo
//  Copyright © 2018 All rights reserved.
//
#pragma warning (disable : 4996)
#include "jpegdecoder.h"
#include "TypeDef.h"
#include <iostream>
#include <iomanip>
#include <bitset>
#include <cstring>

const double pi=3.1415926535897932384626433832795;
using namespace std;
// default constructor
jpeg_decoder::jpeg_decoder() : zigZagStart(0), zigZagEnd(63), losslessFormat(false) {
	// Debug messages used by parseSeg to tell us which segment we're at

	// intialize the rgb pointer to null so that you can create it later
	m_rgb = NULL;

	// initialize the m_Y and m_Cb and m_Cr pointers temp space to NULL
	m_Y = m_Cb = m_Cr = NULL;

	// init the Y, Cb, Cr picture buffer into null:
	m_YPicture_buffer = m_CbPicture_buffer = m_CrPicture_buffer = NULL;


#if PROFILE_JPEG_DECODE_PIC
	auto startTime = std::chrono::high_resolution_clock::now();
#endif

	// init the coefficients buffer into null:
	count_block_Y = count_block_Cb = count_block_Cr = 0;
}

// constructor
jpeg_decoder::jpeg_decoder(std::string filename) : zigZagStart(0), zigZagEnd(63), losslessFormat(false) {
    // Debug messages used by parseSeg to tell us which segment we're at
    segNames[0x00] = std::string("Baseline DCT; Huffman");
    segNames[0x01] = std::string("Extended sequential DCT; Huffman");
    segNames[0x02] = std::string("Progressive DCT; Huffman");
    segNames[0x03] = std::string("Spatial lossless; Huffman");
    segNames[0x04] = std::string("Huffman table");
    segNames[0x05] = std::string("Differential sequential DCT; Huffman");
    segNames[0x06] = std::string("Differential progressive DCT; Huffman");
    segNames[0x07] = std::string("Differential spatial; Huffman");
    segNames[0x08] = std::string("[Reserved: JPEG extension]");
    segNames[0x09] = std::string("Extended sequential DCT; Arithmetic");
    segNames[0x0A] = std::string("Progressive DCT; Arithmetic");
    segNames[0x0B] = std::string("Spatial lossless; Arithmetic");
    segNames[0x0C] = std::string("Arithmetic coding conditioning");
    segNames[0x0D] = std::string("Differential sequential DCT; Arithmetic");
    segNames[0x0E] = std::string("Differential progressive DCT; Arithmetic");
    segNames[0x0F] = std::string("Differential spatial; Arithmetic");
    segNames[0x10] = std::string("Restart");
    segNames[0x11] = std::string("Restart");
    segNames[0x12] = std::string("Restart");
    segNames[0x13] = std::string("Restart");
    segNames[0x14] = std::string("Restart");
    segNames[0x15] = std::string("Restart");
    segNames[0x16] = std::string("Restart");
    segNames[0x17] = std::string("Restart");
    segNames[0x18] = std::string("Start of image");
    segNames[0x19] = std::string("End of image");
    segNames[0x1A] = std::string("Start of scan");
    segNames[0x1B] = std::string("Quantization table");
    segNames[0x1C] = std::string("Number of lines");
    segNames[0x1D] = std::string("Restart interval");
    segNames[0x1E] = std::string("Hierarchical progression");
    segNames[0x1F] = std::string("Expand reference components");
    segNames[0x20] = std::string("JFIF header");
    segNames[0x21] = std::string("[Reserved: application extension]");
    segNames[0x22] = std::string("[Reserved: application extension]");
    segNames[0x23] = std::string("[Reserved: application extension]");
    segNames[0x24] = std::string("[Reserved: application extension]");
    segNames[0x25] = std::string("[Reserved: application extension]");
    segNames[0x26] = std::string("[Reserved: application extension]");
    segNames[0x27] = std::string("[Reserved: application extension]");
    segNames[0x28] = std::string("[Reserved: application extension]");
    segNames[0x29] = std::string("[Reserved: application extension]");
    segNames[0x2A] = std::string("[Reserved: application extension]");
    segNames[0x2B] = std::string("[Reserved: application extension]");
    segNames[0x2C] = std::string("[Reserved: application extension]");
    segNames[0x2D] = std::string("[Reserved: application extension]");
    segNames[0x2E] = std::string("[Reserved: application extension]");
    segNames[0x2F] = std::string("[Reserved: application extension]");
    segNames[0x30] = std::string("[Reserved: JPEG extension]");
    segNames[0x31] = std::string("[Reserved: JPEG extension]");
    segNames[0x32] = std::string("[Reserved: JPEG extension]");
    segNames[0x33] = std::string("[Reserved: JPEG extension]");
    segNames[0x34] = std::string("[Reserved: JPEG extension]");
    segNames[0x35] = std::string("[Reserved: JPEG extension]");
    segNames[0x36] = std::string("[Reserved: JPEG extension]");
    segNames[0x37] = std::string("[Reserved: JPEG extension]");
    segNames[0x38] = std::string("[Reserved: JPEG extension]");
    segNames[0x39] = std::string("[Reserved: JPEG extension]");
    segNames[0x3A] = std::string("[Reserved: JPEG extension]");
    segNames[0x3B] = std::string("[Reserved: JPEG extension]");
    segNames[0x3C] = std::string("[Reserved: JPEG extension]");
    segNames[0x3D] = std::string("[Reserved: JPEG extension]");
    segNames[0x3E] = std::string("Comment");
    segNames[0x3F] = std::string("[Invalid]");
    
    // set the fileName of the input with jpeg
    jpeg_filename = filename;
    
    // intialize the rgb pointer to null so that you can create it later
    m_rgb = NULL;
    
    // initialize the m_Y and m_Cb and m_Cr pointers temp space to NULL
    m_Y = m_Cb = m_Cr = NULL;
    
    // init the Y, Cb, Cr picture buffer into null:
    m_YPicture_buffer = m_CbPicture_buffer = m_CrPicture_buffer = NULL;
    
    
#if PROFILE_JPEG_DECODE_PIC
    auto startTime = std::chrono::high_resolution_clock::now();
#endif
    
    // init the coefficients buffer into null:
    count_block_Y = count_block_Cb = count_block_Cr = 0;
    readFile();

#if PROFILE_JPEG_DECODE_PIC
    auto endTime = std::chrono::high_resolution_clock::now();
    cout << "Decoding elapsed time: " <<  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " milliseconds" << endl;
#endif
    
}

// JPEG Decoder Destructor
jpeg_decoder::~jpeg_decoder() {
    deleteDCTComponentPointers();
    
    // TODO: you need to check that this deletion does not cause any erro
    // delete the raw pointer rgb
    delete [] m_rgb;
    
    // TODO: you need to check this
    deleteRawPictureBufferPointers();
    // delete the raw pointer for m_Y, m_Cb, m_Cr
    delete [] m_Y;
    delete [] m_Cb;
    delete [] m_Cr;
}

void jpeg_decoder::readFile()
{
    // Open the requested file, keep parsing blocks until we run
    // out of file, then close it.
    fp = fopen(jpeg_filename.c_str(), "rb");
    
    if (fp) {
        while(parseSeg() == JPEG_SEG_OK);
        // close the file
        fclose(fp);
    }
    else {
        perror("JPEG parse error");
        printf("Will pause for you to press Enter and exiting afterwards \n");
        getchar();
        exit(0);
    }
    
    
}

int jpeg_decoder::parseSeg()
{
    if (!fp) {
        printf("File failed to open.\n");
        return JPEG_SEG_ERR;
    }
    
    // file position
    uint_32 fpos = ftell(fp);
    
    // 16 bits word -- equivalent to a marker
    uint_16 id = READ_WORD(), size;
    
    if (id < 0xFFC0)
    {
        printf("Segment ID expected, not found.\n");
        return JPEG_SEG_ERR;
    }
    
#if PRINT_FIND_SEGMENTS_DECODER
    printf(
           "Found segment at file position %d: %s 0x%4x\n",
           fpos, segNames[id-0xFFC0].c_str(), id);
#endif
    
    switch (id) {
            // Application specific segement will be ignored (not useful for our purpose)
            // The SOI and EOI segments are the only ones not to have
            // a length, and are always a fixed two bytes long; do
            // nothing to advance the file position
            // Start of Image (SOI)
        case 0xFFD8:
            break;
            
            // End of Image (EOI)
        case 0xFFD9:
            // cout << "Size of TCOFF_Y = " << tCoeff_Y.size() << " x " << tCoeff_Y[0].size() << endl;
            // cout << components.size() << endl;
            
            // Progressive processing need after storing tcoeff: Dequantization, Inverse Zigzag, IDCT & Add block subsampling
            if (progressive_Huff_Format) {
                final_process_progressive();
            }
            // Convert the Y'CbCr image into RGB
            ycrcb_to_rgb24_image();
            
            return JPEG_SEG_EOF;
            
            // The DHT segment defines a Huffman table. The handler should
            // read exactly as many bytes from the file as are in the
            // segment; if not, something's gone wrong
        case 0xFFC4:
            size = READ_WORD() - 2;
            if (readHuffmanTables(size) != size) {
                cout<< "Unexpected end of DHT segment" << endl;
                return JPEG_SEG_ERR;
            }
            break;
            
            
        case 0xFFDB: // DQT = Define Quantization Tables
            size = READ_WORD() - 2;
            if (readQuantizationTables(size) != size) {
                cout<< "Unexpected end of DQT segment" << endl;
                return JPEG_SEG_ERR;
            }
            break;
            
        case 0xFFC2: // SOF2 = Progressive DCT; Huffman
            progressive_Huff_Format = true;
        case 0xFFC0: // SOF0 = Start of frame (Baseline DCT --
            // if you want to add more (lossless..etc), you need to send the maker, and probe other types
            
            size = READ_WORD();
            if (readFrameHeader(size) != size) {
                cout<< "Unexpected end of SOF segment" << endl;
                return JPEG_SEG_ERR;
            }
            // Initialize the necessary positions and buffers for storing the Y, Cb, Cr components:
            // if(!progressive_Huff_Format || counter_progressive == 0)
            initPositionsBuffersForPictureBuffer();
            break;
            
        case 0xFFDD: // DRI = Define Restart Interval
            // DRI is always 4 bytes long
            // We are currently just ignoring restart markers
            for (int i = 0; i < 4; ++i) {
                fgetc(fp);
            }
            break;
            
        case 0xFFDE: // DHP = Define Hierarchical Progression
        case 0xFFDF: // EXP = Expand Reference Components
            throw "Hierarchical JPEG not supported yet";
            
            // An SOS (Start of Scan) segment has a length determined only by the
            // length of the bitstream; for now, assume it's the rest
            // of the file less the two-byte EOI segment
        case 0xFFDA: // SOS = Start Of Scan
            // cout << ftell(fp) << endl;
            // counter_SOS++;
            size = READ_WORD();
            if (readScanHeader(size) != size) {
                cout<< "Unexpected end of SOS segment" << endl;
                return JPEG_SEG_ERR;
            }
            
            // Scan starts immediately after header
            // Try a new approach
            // Include both sequential and progressive
            readImageEntryPoint();
            //counter for counting the number of scans (i.e. number of 0xFFDA)
            //if (progressive_Huff_Format)
            break;
            // Any other segment has a length specified at its start,
            // so skip over that many bytes of file
            
            // Reserved for application segments
        case 0xFFE0:
        case 0xFFE1:
        case 0xFFE2:
        case 0xFFE3:
        case 0xFFE4:
        case 0xFFE5:
        case 0xFFE6:
        case 0xFFE7:
        case 0xFFE8:
        case 0xFFE9:
        case 0xFFEA:
        case 0xFFEB:
        case 0xFFEC:
        case 0xFFED:
        case 0xFFEE:
        case 0xFFEF:
            size = READ_WORD();
            application_size += size + 2; // Application length
            fseek(fp, size - 2, SEEK_CUR);
            break;
        
            
        default:
            size = READ_WORD();
            fseek(fp, size-2, SEEK_CUR);
            break;
    }
    
    return JPEG_SEG_OK;
}

uint_16 jpeg_decoder::readHuffmanTables(uint_16 tableLengthFromBitStream)
{
    
    // first two bytes are Huffman table length
    // third byte is separated in two nibbles. First nibble is Table Class,
    // second nibble is destination identifier or ID
    // Next 16 bytes are number of elements coded with 1-16 bits
    
    // First two bytes-table length
    tableLengthFromBitStream += 2;
    uint_8  tableID = 0; // Specifies one of component: 0 for luminance and 1 for chrominance
    uint_8  tableClass = 0; // Specifies is it DC element or AC element of table. 0-DC element 1-AC element
    uint_8  huffmanTableOptions = 0; // I will decompose this in two nibbles
    
    
    // counter of how many bytes have been read (2 bytes are read for the length of huffman table)
    int bytes_read = 2;
    
    
    while(bytes_read < tableLengthFromBitStream)
    {
        // read huffmantable options
        huffmanTableOptions = fgetc(fp);
        bytes_read++;
        
        tableID = huffmanTableOptions & 0x0F; // Right most 4 bits of the byte
        tableClass = huffmanTableOptions >> 4; // Left most 4 bits of the byte
        
#if PRINT_HUFFMAN_TABLE
        if(tableClass)
        {
            printf("Huffman table ID #%02X class 01 (AC Table):\n", tableID);
        }
        else {
            printf("Huffman table ID #%02X class 00 (DC Table):\n", tableID);
        }
        cout << "Table Length: " << std::dec << tableLengthFromBitStream << endl;
#endif
        
        // Looking for tableID in tables
        HuffmanTable* table = 0;
        /*for (uint i = 0; i< huffmanTables.size(); i++) {
         HuffmanTable* t = huffmanTables[i];
         if (t->tableID == tableID && t->tableClass == tableClass) {
         table = huffmanTables[i];
         break;
         }
         }*/
        
        // Not found, create a new table
        if (table == 0) {
            table = new HuffmanTable();
            table->tableID = tableID;
            table->tableClass = tableClass;
            table->tableSegmentLengthFromBitstream = tableLengthFromBitStream;
            huffmanTables.push_back(table);
        }
        
        // Next 16 bytes are number of elements coded with 1-16 bits
        uint_8 codeLengths[16];
        for(int j = 0; j < 16; ++j)
        {
            codeLengths[j] = fgetc(fp);
            bytes_read++;
            
            // save it:
            table->number_of_codes_for_each_1to16[j] = codeLengths[j];
        }
        
        uint_32 code = 0; // Huffman code which will be connected with element
        uint_8 element = 0; // Read element from file
        
        // Remaining bytes are the data values to be mapped
        // Build the Huffman map of (length, code) -> value
        for(int i = 0; i < 16; ++i)
        {
            for(int j = 0; j < codeLengths[i]; ++j)
            {
                // element is the actual unique element in the alphabet
                element = fgetc(fp);
                bytes_read++;
                
                table->codes[element] = code;
                table->codeLengths[element] = i + 1; // the length is at least 1 and at most 16
                
                // so far used for printing purposes
                // <Length, code> = element
                table->huffData[huffKey(i + 1, code)] = element;
                code++; //Elements on the same tree depth have code incremented by one
            }
            
            code <<= 1; // multiply by 2 for next iteration
        }
        
#if PRINT_HUFFMAN_TABLE
        // Once the map has been built, print it out
        std::map<huffKey, uint_8>::iterator iter;
        for (iter = table->huffData.begin();iter != table->huffData.end(); ++iter) {
            
            // Print Code - Its Length : Equivalent letter
            printf("    %04X at length %d = %02X\n",
                   iter->first.second, iter->first.first, iter->second);
        }
#endif
        
    } // end while loop
    
    // Compare the length of the table with bytes read minus the table length bytes (=2)
    return bytes_read-2;
} // end readHuffmanTables


uint_16 jpeg_decoder::readQuantizationTables(uint_16 tableLengthFromBitStream)
{
    // First two bytes from this point is Quantization Table length
    // Third byte is divided in two nibbles.
    // First nible is Quantization element table precision 8 or 16 bits
    // Second nibble is Quantization table destination identifier
    // Next 64 bytes are Quantization table element if the size of elements are 8 bits, or 128 bytes if the size of elements is two bytes
    // First two bytes-table length
    tableLengthFromBitStream += 2;
    uint_8  tableOptions = 0; // Specifies one of component: 0 for luminance and 1 for chrominance
    uint_8  tableIdentifier = 0; // Specifies is it DC element or AC element of table. 0-DC element 1-AC element
    uint_8  sizeOfElements = 0; // I will decompose this in two nibbles
    int tableData[64]; // 64 elements per table
    
    // counter of how many bytes have been read (2 bytes are read for the length of huffman table)
    int bytes_read = 2;
    
    while(bytes_read < tableLengthFromBitStream)
    {
        tableOptions = fgetc(fp);
        sizeOfElements  = tableOptions >> 4;     // Left most 4 bits of the byte
        tableIdentifier = tableOptions & 0x0F;   // Right most 4 bits of the byte
        
        QuantizationTable qTable(false, false);
        qTable.tableID = tableIdentifier;
        qTable.tableLength = tableLengthFromBitStream;
        
        
#if PRINT_QUANTIZATION_TABLE
        if(tableIdentifier)
        {
            printf("Quantization table ID #%02X class 01 (chrominance):\n", tableIdentifier);
        }
        else {
            printf("Quantization table ID #%02X class 00 (luminance):\n", tableIdentifier);
        }
        cout << "Table Length: " << std::dec << tableLengthFromBitStream << endl;
#endif
        
        // Elements are represented by 8 bits
        if (sizeOfElements == 0) {
            
#if PRINT_QUANTIZATION_TABLE
            cout << "Precision: 8 bits" << endl;
#endif
            for (int i = 0; i < 64; ++i) {
                tableData[i] = fgetc(fp); // Read a byte
            }
            
            // 1 byte and 64 bytes for the tableData
            bytes_read += 65;
        }
        // 16-bit elements
        else if(sizeOfElements == 1) {
            
#if PRINT_QUANTIZATION_TABLE
            cout << "Precision: 16 bits" << endl;
#endif
            
            for (int i = 0; i < 64; ++i) {
                
                tableData[i] = READ_WORD(); // Read 16 bits
            }
            
            // 1 byte and 64*2=128 bytes for the tableData
            bytes_read += 129;
            
        }
        
        // do inverse zigzag scanning
        inverseZigZagCoding (tableData, qTable.quantizationTableData);
        
        // push back
        quantizationTables.push_back(qTable);
        
#if PRINT_QUANTIZATION_TABLE
        cout << "Quantization Table: " << endl;
        for(int i = 0; i < 8; ++i)
        {
            for(int j = 0; j < 8; ++j)
            {
                if(qTable.quantizationTableData[i][j] < 10)
                {
                    cout << " " << qTable.quantizationTableData[i][j] << " ";
                }
                else
                {
                    cout  << qTable.quantizationTableData[i][j] << " ";
                }
                
            }
            cout << "\n";
        }
#endif
    } // end while loop
    // you compare the length of the table with bytes read minus the table length bytes (=2)
    return bytes_read-2;
} // end readQuantizationTables()


void jpeg_decoder::inverseZigZagCoding(int *array, int (*matrix)[8]) {
    int x_indices[64] = {0, 1, 0, 0, 1, 2, 3, 2,
        1, 0, 0, 1, 2, 3, 4, 5,
        4, 3, 2, 1, 0, 0, 1, 2,
        3, 4, 5, 6, 7, 6, 5, 4,
        3, 2, 1, 0, 1, 2, 3, 4,
        5, 6, 7, 7, 6, 5, 4, 3,
        2, 3, 4, 5, 6, 7, 7, 6,
        5, 4, 5, 6, 7, 7, 6, 7};
    
    int y_indices[64] = {0, 0, 1, 2, 1, 0, 0, 1,
        2, 3, 4, 3, 2, 1, 0, 0,
        1, 2, 3, 4, 5, 6, 5, 4,
        3, 2, 1, 0, 0, 1, 2, 3,
        4, 5, 6, 7, 7, 6, 5, 4,
        3, 2, 1, 2, 3, 4, 5, 6,
        7, 7, 6, 5, 4, 3, 4, 5,
        6, 7, 7, 6, 5, 6, 7, 7};
    
    for (int i = zigZagStart; i <= zigZagEnd; ++i)
    {
        matrix[y_indices[i]][x_indices[i]] = array[i-zigZagStart];
    }
    
} // end inverseZigZagCoding

void jpeg_decoder::inverseZigZagScanning(int out_matrix[8][8], const int input_matrix[8][8]) {
#if DEBUGLEVEL > 12
    cout << "Block:--------------" << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            cout << std::dec <<  input_matrix[i][j] << " ";
        }
        
        cout << "\n";
    }
#endif
    
    // traverse order
    int ZigZagArray[64] =
    {
        0,   1,   5,  6,   14,  15,  27,  28,
        2,   4,   7,  13,  16,  26,  29,  42,
        3,   8,  12,  17,  25,  30,  41,  43,
        9,   11, 18,  24,  31,  40,  44,  53,
        10,  19, 23,  32,  39,  45,  52,  54,
        20,  22, 33,  38,  46,  51,  55,  60,
        21,  34, 37,  47,  50,  56,  59,  61,
        35,  36, 48,  49,  57,  58,  62,  63,
    };
    
    
    // copy the matrix line by line over to temp
    int temp[64] = {0};
    for (int i = 0; i < 8; ++i) {
        for(int j = 0; j < 8; ++j)
        {
            int idx = j + i*8;
            temp[idx] = input_matrix[i][j];
        }
    }
    
    // do inverse zig zag scanning
    for (int i = 0; i < 8; ++i) {
        for(int j = 0; j < 8; ++j)
        {
            int idx = j + i*8;
            out_matrix[i][j] = temp[ZigZagArray[idx]];
        }
    }
    
    
#if DEBUGLEVEL > 12
    cout << "Block After:--------------" << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            cout << std::dec <<  out_matrix[i][j] << " ";
        }
        cout << "\n";
    }
#endif
} // end inverse zig zag scanning


uint_16 jpeg_decoder::readFrameHeader(uint_16 headerLengthFromBitStream)
{
    
    // bytes Read
    int bytes_read = 2; // 2 bytes for the length
    
    // Sample Precision:
    uint_8 sample_precision = fgetc(fp);
    
    // Taking picture height and width, which are stored in two bytes for each dimension
    uint_16 pictureHeight, pictureWidth;
    
    // Height in 2 bytes
    pictureHeight = READ_WORD();
    
    // Width in 2 bytes
    pictureWidth = READ_WORD();
    
    
    // Taking information for number of components and components data
    // Every component has 4 pillars: ComponentID, Horizontal sampling, Vertical Sampling, and Quantization table ID
    
    // Number of components:
    uint_8 numberOfComponents = fgetc(fp);
    
    // Increment bytesRead
    bytes_read += 6;
    
    if(headerLengthFromBitStream != numberOfComponents*3 + 8)
    {
        cout << "Bad SOF header length: expected " << (numberOfComponents*3+8) << "but given " << headerLengthFromBitStream << endl;
        return JPEG_SEG_ERR;
    }
    
    
    for(int i = 0; i < numberOfComponents; ++i)
    {
        uint_8 componentData[3];
        componentData[0] = fgetc(fp); // component id
        componentData[1] = fgetc(fp); // here you can get the horizontal (left most part)
        // and vertical factor (right most part)
        componentData[2] = fgetc(fp); // qTableID for the component
        
        // Increment bytesRead
        bytes_read += 3;
        
        // TODO: losless format check
        
        bool foundQuantizationTable = false;
        
        for(int j = 0; j < quantizationTables.size(); ++j)
        {
            if(quantizationTables[j].tableID == componentData[2])
            {
                Component newComponentMember(componentData[0], (componentData[1]>>4), (componentData[1] & 0x0F), componentData[2], quantizationTables[j]);
                this->components.push_back(newComponentMember);
                foundQuantizationTable = true;
            }
        } // end for loop
        
        if (!foundQuantizationTable) {
            cout << "Component with ID " << componentData[0] << " specifies unknown quantization table - skipping." << endl;
        }
    } // end for loop
    
    // Calculate scaling factors - much easier to work with
    int maxHFactor=1, maxVFactor=1;
    for (uint i = 0; i < components.size(); ++i) {
        if (components[i].HFactor > maxHFactor) maxHFactor = components[i].HFactor;
        if (components[i].VFactor > maxVFactor) maxVFactor = components[i].VFactor;
    }
    for (uint i = 0; i < components.size(); ++i) {
        components[i].HScale = maxHFactor / components[i].HFactor;
        components[i].VScale = maxVFactor / components[i].VFactor;
    }
    if (maxHFactor == 1 && maxVFactor == 1)
        this->hasSubSampling = false;
    else
        this->hasSubSampling = true;
    
    // Set the picture parameters
    jpegImageWidth = pictureWidth;
    jpegImageHeight = pictureHeight;
    jpegImageSamplePrecision = sample_precision;
    
    
    // set the mcu width and heights:
    mcu_width  = 8 * maxHFactor;
    mcu_height = 8 * maxVFactor;
    
    //set the mcu cols and rows
    mcu_rows = (jpegImageHeight + mcu_height - 1) / mcu_height ;
    mcu_cols = (jpegImageWidth  + mcu_width - 1) / mcu_width ;
    
    // Calculate the upscaled width and height to be N multiples of the MCU size:
    
    // Upscaled width and height
    upscale_width =  ceil(1.0* jpegImageWidth / mcu_width) * mcu_width;
    upscale_height = ceil(1.0* jpegImageHeight / mcu_height) * mcu_height;
    
    
    
    // TODO: unsed remove this
    total_block_Y = upscale_height * upscale_width / 64;
    total_block_C = ceil(1.0 * total_block_Y / (components.at(COMPONENT_Y).VFactor));
    
#if PRINT_FRAME_HEADER_SOF
    // Header Length
    cout << "============= MARKER: SOF" << endl;
    cout << "Header Length: " << headerLengthFromBitStream << endl;
    printf("Sample precision: %d\n", jpegImageSamplePrecision);
    cout << "Has subsampling: " << boolalpha << hasSubSampling << endl;
    cout << "Image Size (width*height): " << jpegImageWidth << "*" << jpegImageHeight << endl;
    cout << "Number of components: " << components.size() << endl;
    
    for(int i = 0; i < components.size(); ++i)
    {
    std:string component_type;
        switch(i)
        {
            case COMPONENT_Y:
                component_type = " (Lumin:  Y)";
                break;
            case COMPONENT_Cb:
                component_type = " (Chrom: Cb)";
                break;
            case COMPONENT_Cr:
                component_type = " (Chrom: Cr)";
                break;
        }
        cout << "Component(#" << i << "): ID = " << components.at(i).componentID << component_type << ", SubSamp "
        << components.at(i).HScale << "*" << components.at(i).VScale << ", Sel QT ID = " << components.at(i).componentQTableID << endl ;
    }
    
#endif
    
    return bytes_read;
}


uint_16 jpeg_decoder::readScanHeader(uint_16 headerLength) {
    
    /* --- Entropy-coded image data
     After FF DA marker starts header followed by data
     Header contains:
     - Header length (not data, just header)
     - Number of color components e.g. 3 for Y'CbCr, 4 for CMYK
     - Component IDs (2 bytes per component)
     - ZigZag definition (3 bytes)
     After header immediately comes data ending with marker FF D9
     */
    
    if (components.size() == 0) {
        throw "SOF header not present or specifies 0 components.";
    }
    
    // bytes Read
    int bytes_read = 2; // 2 bytes for the length
    
    // number of components //Ns in standard
    numberOfComponents = fgetc(fp);
    
    // catch the exception from CMYK
    if(numberOfComponents > MAX_NUMBER_SUPPORTED_COLOR_COMPONENTS)
    {
        // logCMYKErrorPictures();
    }

    
    // Increment bytes_read
    bytes_read += 3;
    
    if (headerLength != numberOfComponents*2 + 6){
        cout << "Bad SOF header length: expected " << (numberOfComponents*2 + 6) << ", but given " << headerLength << endl;
        return JPEG_SEG_ERR;
    }
    
    // For progressive: potential 1 component per scan...
    /*if (numberOfComponents != this->components.size()) {
     cout << "Number of components in SOS header is" << numberOfComponents << ", but in SOF header is " << this->components.size();
     cout << " -- We will trust SOF header, but this probably wont work :(" << endl;
     numberOfComponents = this->components.size();
     }*/
    
    if (numberOfComponents > ETF_FORMAT_MAX_COMPONENTS) {
        cout << "Specified " << numberOfComponents << " components, maximum is " << ETF_FORMAT_MAX_COMPONENTS;
        cout << " -- We will trust SOF header, but this probably wont work :(" << endl;
        numberOfComponents = ETF_FORMAT_MAX_COMPONENTS;
    }
    
    // Read components and their huffman tables
    for (int i = 0; i < numberOfComponents; ++i) {
        uint_8 currentID = fgetc(fp);
        componentID.push_back(currentID);
        uint_8 tableID = fgetc(fp);
        
        // Increment bytes_read
        bytes_read +=2;
        
        // Check component in components list
        if (!progressive_Huff_Format) {
            if (this->components.at(i).componentID != componentID[i]) {
                cout << "Component ID in SOS header is " << componentID[i] << ", but in SOF header is " << this->components[i].componentID;
                cout << " -- We will trust SOF header" << endl;
            }
        }
        
        
        // Find AC and DC Huffman table in tables list
        uint_8 tableDC = tableID >> 4; // left most part
        uint_8 tableAC = tableID & 0x0F; // right most part
        componentTablesDC[currentID - 1] = componentTablesAC[currentID - 1] = 0;
        
        // Find component huffman table
        for(int j = 0; j < huffmanTables.size(); ++j){
            HuffmanTable* tempHuffmanTable = huffmanTables[j];
            
            if(tempHuffmanTable->tableID == tableDC && !tempHuffmanTable->tableClass){
                componentTablesDC[currentID - 1] = tempHuffmanTable;
            }
            
            if(tempHuffmanTable->tableID == tableAC && tempHuffmanTable->tableClass){
                componentTablesAC[currentID - 1] = tempHuffmanTable;
            }
        }
        
        // DC Huffman tables not found
        if (componentTablesDC[currentID - 1] == 0) {
            // Get first usable table
            for (uint j = 0; j < huffmanTables.size(); ++j) {
                HuffmanTable* t = huffmanTables[j];
                if (t->tableClass == 0) {
                    componentTablesDC[currentID - 1] = t;
                    break;
                }
            }
            
            if (componentTablesDC[currentID - 1] == 0) {
                cout << "File contains no DC Huffman tables!" << endl;
                return JPEG_SEG_ERR;
            }
            
        } // DC not found
        
        // AC Huffman tables not found
        if (componentTablesAC[currentID - 1] == 0) {
            // Get first usable table
            for (uint j = 0; j < huffmanTables.size(); ++j) {
                HuffmanTable* t = huffmanTables[j];
                if (t->tableClass == 1) {
                    componentTablesAC[currentID - 1] = t;
                    break;
                }
            }
            if (componentTablesAC[currentID - 1] == 0) {
                // cout << "File contains no AC Huffman tables!" << endl; // Can happen in progressive JPEG!
                componentTablesAC[currentID - 1] = componentTablesDC[currentID - 1];
            }
            
            //cout << "SOS header specifies inexistant Huffman table " << tableAC << " - using " << componentTablesAC[i]->tableID << endl;
        } // AC not found
        
    } // end loop on components
    
    
    // Start and end point for zig-zag coding
    zigZagStart = fgetc(fp); // Ss in standard
    zigZagEnd   = fgetc(fp); // Se in standard
    
    //FOR PROGRESSIVE
    if (zigZagStart > 1) {
        EOB_run = 0;
    }
    //cout << zigZagStart << "  " << zigZagEnd << endl;
    // Increment bytes_read
    bytes_read += 2;
    
    // TODO: Bit approximation for progressive JPEG
    //unsigned char dummy;
    ssa = fgetc(fp);
    approximationH = ssa >> 4;    // Ah in standard
    approximationL = ssa & 0x0F;  // Al in standard
    // Increment bytes_read
    bytes_read ++;
    
#if PRINT_SOS
    cout << "============= MARKER: SOS" << endl;
    cout << "Header Length: " << headerLength << endl;
    cout << "Zigzag start: " << static_cast<int>(zigZagStart) << ", ZigZag End: " << static_cast<int>(zigZagEnd) << endl;
    cout << "Successive approximation bit position HIGH: " << static_cast<int>(approximationH) << ", Successive approximation LOW: " << static_cast<int>(approximationL) << endl;
    printf("Number of components in scan: %d\n", numberOfComponents);
    
    for(int i = 0; i < numberOfComponents; ++i){
        printf("Component[%d]: selector = %d, table = %d (DC), %d (AC) \n", i, components.at(i).componentID,
               componentTablesDC[componentID[i] - 1]->tableID, componentTablesAC[componentID[i] - 1]->tableID);
    }
    
#endif
    // you compare the length of the table with bytes read minus the table length bytes (=2)
    return bytes_read-2;
}

void jpeg_decoder::readProgressiveImageEntryPoint() {
}

void jpeg_decoder::readImageEntryPoint() {
    
    if (losslessFormat) {
        // TODO: lossless format is not supported here, but I added the code for it
        // Predictor requires to hold last line in memory
        //        for (uint i=0; i<components.size(); i++)
        //            scanLineCache[i] = new int[jpegImageWidth];
        cout << "-|- ##ERROR## LoslessFormat is not supported! " << endl;
        return;
    }
    
    if(!components.size()){
        cout << "-|- ##ERROR## Components vector is empty - you did not read any components! " << endl;
        return;
    }
    
    // TODO: Add a check here for images that are not multiples of 8
    
    // TODO: do you need to check the maximum here?
    // Y is the original hFactor and vFactor; others are subsampled
    int hFactor = components[COMPONENT_Y].HFactor;
    int vFactor = components[COMPONENT_Y].VFactor;
    
    
    // Set the maximum sample value to use it in clamping values
    // when you are calculating the YCbCr and rgb pixels
    maxSample = pow(2, jpegImageSamplePrecision)-1;
    
    // Create RGB24 buffer (8 bits for Red, 8 bits for Green, 8 bits for Blue):
    if (m_rgb == NULL) {
        int nRGBComponents = 3;
        int height = jpegImageHeight * nRGBComponents;
        int width  = jpegImageWidth  * nRGBComponents;
        
        // Create an RGB buffer with the same resolution even if they are not multiples of 8
        m_rgb = new unsigned char[width * height];
        memset(m_rgb, 0, width*height); // initialize all to zero
        
#if DEBUGLEVEL > 100
        cout << "Size of RGB Buffer is " << width << "x" << height << endl;
#endif
        
    }
    
    // Create m_Y, m_Cb, m_Cr temp storage space on the heap
    // you want this memory space until the end of the program
    // TODO: do you really need to allocate this memory on the heap?
    if (m_Y == NULL) {
        m_Y = new uint_8[64*4]; // at most 64 coefficients and 4 Y components (max 2x2 hFactor*vFactor, each block is 8x8 = 8x8x2x2 = 256 elements)
    }
    
    // Create the Cb if the number of components is 2 or more
    if(m_Cb == NULL && components.size()>=(COMPONENT_Cb+1)) {
        m_Cb = new uint_8[64];
    }
    
    // Create the Cr if the number of components is 3 or more
    if(m_Cr == NULL && components.size()>=(COMPONENT_Cr+1)) {
        m_Cr = new uint_8[64];
    }
    
    
    // Set all previous DC to zero
    for (uint iComponent = 0; iComponent < components.size(); ++iComponent) {
        previousDC[iComponent] = 0;
    }
    
    // Each block is 8x8 and you have 2x2 hFactor*vFactor
    // MCU: minimum coding unit
    int xstride_by_mcu = 8 * hFactor;
    int ystride_by_mcu = 8 * vFactor;
    
    // TODO: you can use this part to support one-channel images (Y only)
    // Don't forget to that block can be either 8 or 16 lines (components.size should normally be 3)
    uint_32 bytes_per_blocklines = static_cast<uint_32>(jpegImageWidth * components.size() * ystride_by_mcu);
    uint_32 bytes_per_mcu        = static_cast<uint_32>(components.size() * xstride_by_mcu);
    
#if DEBUGLEVEL > 20
    cout << "Horizontal factor " << hFactor << ", Vertical factor: " << vFactor << endl;
    cout << "xstride_by_mcu: " << xstride_by_mcu << ", ystride_by_mcu: " << ystride_by_mcu << endl;
    cout << "bytes_per_blocklines: " << bytes_per_blocklines << ", bytes_per_mcu: " << bytes_per_mcu << endl;
#endif
    
#if PRINT_BLOCK_PROGRESS
    int block_number = 0;
#endif
    
    // TODO: Support for images that are not multiples of 8x8, 420 (needs more work)
    //    int lim_x = ceil(1.0*jpegImageWidth/xstride_by_mcu);
    //    int lim_y = ceil(1.0*jpegImageHeight/ystride_by_mcu);
    //    for (int y = 0 ; y < lim_y; ++y) for (int x = 0; x < lim_x; ++x)
    
    // Just the decode the image by 'macroblock' (size is 8x8, 8x16, or 16x16)
    g_nbits_in_reservoir = 0;
    g_reservoir = 0;
    initCurrentPosition();
    
    for (int y = 0 ; y < jpegImageHeight; y+=ystride_by_mcu){
        for (int x = 0; x < jpegImageWidth; x+=xstride_by_mcu){
            //EOB_run = 0;
            // Decode the MCU plane
            if (!progressive_Huff_Format) {
                decode_mcu(hFactor, vFactor, x, y);
            }
            else {
                decode_mcu_progressive(hFactor, vFactor, x, y, componentID);
            }
            // cout << x << "---" << y << endl;
#if PRINT_BLOCK_PROGRESS
            cout << "Block " << (block_number++) << " is done!"  << endl;
#endif
        } // end inner
    }// end outer
    componentID.clear(); // Clear the ID vector for next use
    
#if DEBUGLEVEL > 50
    cout << "RGB Picture just before conversion: " << endl;
    int Height = jpegImageHeight;
    int Width = jpegImageWidth;
    const unsigned char* RGB = m_rgb;
    // Round up the width to the nearest DWORD boundary
    int iNumPaddedBytes = 4 - (Width * 3) % 4;
    iNumPaddedBytes = iNumPaddedBytes % 4;
    for (int y=Height-1; y>=0; y--)
    {
        for (int x=0; x<Width; x++)
        {
            int i = (x + (Width)*y) * 3;
            unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
            printf("%d ", rgbpix);
        }
        
        if (iNumPaddedBytes>0)
        {
            unsigned char pad = 0;
            printf("%d ", pad);
        }
    }
#endif
    
} // end readImageEntryPoint


//  Decoding for Y 2x2:
//  .-------.
//  | 1 | 2 |
//  |---+---|
//  | 3 | 4 |
//  '-------'

// Progressive Mode's Decode MCU
// Each scan may contain different component to scan, than need a buffer componentID to identify which component is scanned during current scan
void jpeg_decoder::decode_mcu_progressive(int componentWidth, int componentHeight, int currentX, int currentY, vector<uint_8> componentID) {
    for (int m = 0; m < componentID.size(); ++m) {
        counter_scan_blockidx = 0;
        // Y component
        if (componentID[m] == 1) {
            for (int y = 0; y < componentHeight; ++y)
            {
                for (int x = 0; x < componentWidth; ++x)
                {
                    int stride = componentHeight * 8;
                    int offset = x * 8 + y * 64 * componentWidth;
                    // Y component: (a data unit for huffman is an 8x8 block)
                    // huffman table process in progressive mode
                    process_huffmann_data_unit_progressive(COMPONENT_Y, currentX, currentY);
                    decode_single_block(offset, stride, COMPONENT_Y, &(m_Y[offset]), currentX, currentY);
                    counter_scan_blockidx++;
                } // end inner loop
            } // end outer loop
        }
        else if (componentID[m] == 2) {
            int stride = 8;
            int offset = 0;
            // Huffman table process in progressive mode
            process_huffmann_data_unit_progressive(COMPONENT_Cb, currentX, currentY);
            decode_single_block(offset, stride, COMPONENT_Cb, m_Cb, currentX, currentY);
        }
        else {// Cr component
            int stride = 8;
            int offset = 0;
            // Huffman table process in progressive mode
            process_huffmann_data_unit_progressive(COMPONENT_Cr, currentX, currentY);
            decode_single_block(offset, stride, COMPONENT_Cr, m_Cr, currentX, currentY);
        }
    }
} // end decode_mcu_progressive
//---------------------------------------------------------------------------------------------------

// Sequential Mode's Decode MCU
void jpeg_decoder::decode_mcu(int componentWidth, int componentHeight,int currentX, int currentY) {
    for (int y = 0; y < componentHeight; ++y )
    {
        for(int x = 0; x < componentWidth; ++x)
        {
            int stride = componentHeight * 8;
            int offset = x * 8 + y * 64 * componentWidth;
            
            // Y component: (a data unit for huffman is an 8x8 block)
            process_huffmann_data_unit(COMPONENT_Y, currentX, currentY);
            decode_single_block(offset, stride, COMPONENT_Y, &(m_Y[offset]), currentX, currentY);
            count_block_Y++;
            
        } // end inner loop
    } // end outer loop
    
    // The rest of the components if they exist:
    for (int iComponent = 1; iComponent < numberOfComponents; ++iComponent) //components.size()
    {
        int stride = 8;
        int offset = 0;
        
        // Cb and Cr:
        process_huffmann_data_unit(iComponent, currentX , currentY);
        
        if(iComponent == COMPONENT_Cb) {
            decode_single_block(offset, stride, iComponent, m_Cb, currentX, currentY);
        }
        else if(iComponent == COMPONENT_Cr) {
            decode_single_block(offset, stride, iComponent, m_Cr, currentX, currentY);
        }
        else {
            cout << "-|- ##ERROR## We only support three color components (Y'CbCr)." << endl;
        }
        
        count_block_Cb++;
        count_block_Cr++;
        
    }
} // end decode_mcu

// Process Huffman Data in Progressive Mode
void jpeg_decoder::process_huffmann_data_unit_progressive(int currentComponent , int currentX, int currentY) {
    
    // We memset it here, as later on we can just skip along, when we have lots
    // of leading zeros, for our AC run length encoding :)
    short DCT_tcoeff[64];
    memset(DCT_tcoeff, 0, sizeof(DCT_tcoeff)); //Initialize DCT_tcoeff
    
    // How many AC elements should we read?
    int ACcount;
    if (zigZagStart != 0) ACcount = zigZagEnd - zigZagStart + 1;
    else ACcount = 0;
    
    // found in huffman table
    bool found = false;
    
    // decodedValue from huffmanTable
    int decodedValue = 0;
    
    // DC/AC Refine
    bool progressive_DC_refine = false;
    bool progressive_AC_refine = false;
    
#if DEBUGLEVEL > 40
    printDCHuffmanCodes(currentComponent);
    printACHuffmanCodes(currentComponent);
#endif
    
    // DC Scans
    if (zigZagStart == 0) {
        if (approximationH != 0) progressive_DC_refine = true;
        if (!progressive_DC_refine) { // DC First Scans
            //    According to G.2 "In order to avoid repetition, detail flow diagrams
            //    of progressive decoder operation are not included. Decoder operation is
            //    defined by reversing the function of each stop described in the encoder
            //    flow charts, and performing the steps in reverse order."
            for (int codeLength = 1; codeLength <= 16; ++codeLength)
            {
                // Keep grabbing one bit at a time till we find one thats a huffman code
                int code = LookNBits(codeLength);
                // found in huffman table
                found = false;
                
                // Check if its one of our huffman codes
                // Current Huffman table (Initially DC)
                // We decode the first DC coeffient that same way as in a sequential
                // scan except for the point transform according to G.1.2.1
                if (is_exist_in_huffman_codes(code, codeLength, currentComponent, decodedValue)) {
                    // you looked ahead, skip those bits
                    SkipNBits(codeLength);
                    // set the found in huffman table to true
                    found = true;
                    // The decoded value is the number of bits we have to read in next
                    int numDataBits = decodedValue;
                    
#if DEBUGLEVEL > 50
                    cout << "The decoded value is the number of bits we have to read in next equals to " << numDataBits << endl;
#endif
                    // We know the next k bits are for the actual data
                    if (numDataBits == 0)
                    {
                        DCT_tcoeff[0] = previousDC[currentComponent];
                        
#if DEBUGLEVEL > 50
                        cout << "First DC coefficient is " << DCT_tcoeff[0] << endl;
#endif
                    }
                    else{
                        // residual DC:
                        short residualDC = GetNBits(numDataBits);
#if DEBUGLEVEL > 50
                        cout << "Delta of DC coefficient " << residualDC << endl;
#endif
                        residualDC = DetermineSign(residualDC, numDataBits);
#if DEBUGLEVEL > 50
                        cout << "Delta of DC coefficient with sign " << residualDC << endl;
#endif
                        DCT_tcoeff[0] = residualDC + previousDC[currentComponent];
                        previousDC[currentComponent] = DCT_tcoeff[0];
                        
#if DEBUGLEVEL > 50
                        //cout << "Component: " << currentComponent  << ", DC coefficient is: " << DCT_tcoeff[0] << " and prev DC coeff is updated" << endl;
                        //if (ftell(fp) >= 1930) {
                        cout << "DC Postion in Byte: " << ftell(fp) << endl;
                        cout << "Code: " << code << " -- Length: " << codeLength << endl;
                        printf("Code: %X \n", code);
                        cout << "Component:" << currentComponent << endl;
                        cout << "Decoded value in hex: " << std::hex << numDataBits << endl;
                        cout << std::dec << "DC coeff: " << DCT_tcoeff[0] << endl;
                        cout << "residual DC: " << residualDC << endl;
                        cout << "postion X:" << currentX << ", position Y:" << currentY << endl;
                        //}
#endif
                    } // end else
                    // Found so we can exit out
                    break;
                } // end if is_exist
            } // end for loop on codeLengths
            
            if (!found){
                cout << "-|- ##ERROR## (DC case) Code value not found in Huffman table: " << endl;
                // Log it
                // logErrorPictures();
                return;
            }
            
            // No AC coefficients required?
            if (ACcount == 0 || losslessFormat) { // i.e. DC output to m_DCT
                Component comp = components.at(currentComponent);
                comp.m_DCT[0] = (DCT_tcoeff[0] << approximationL);
                return;
            }
        }
        else {// DC Refine Scans
            //  This part saves the DC coefficient for a data unit in refining
            //  DC scans for the component.
            // Retrieve Data firstly
            switch (currentComponent) {
                case COMPONENT_Y:
                    DCT_tcoeff[0] = tCoeff_Y[currentY_dct[currentComponent]][currentX_dct[currentComponent]];
                    break;
                case COMPONENT_Cb:
                    DCT_tcoeff[0] = tCoeff_Cb[currentY_dct[currentComponent]][currentX_dct[currentComponent]];
                    break;
                case COMPONENT_Cr:
                    DCT_tcoeff[0] = tCoeff_Cr[currentY_dct[currentComponent]][currentX_dct[currentComponent]];
                    break;
                default:
                    DCT_tcoeff[0] = tCoeff_Y[currentY_dct[currentComponent]][currentX_dct[currentComponent]];
                    break;
            } // end switch
            Component comp = components.at(currentComponent);
            if (GetNBits(1) != 0) {
                DCT_tcoeff[0] |= (1 << approximationL);
                //DCT_tcoeff[0] = 128;
            }
            comp.m_DCT[0] = DCT_tcoeff[0];
        }
    }
    else {
        // NOTE: rrrr = count_0  ssss = size_val
        // Zigzag != 0 --> AC Scans
        int coeff_counter;
        if (zigZagStart != 0 && approximationH != 0) progressive_AC_refine = true;
        bool progressive_done = false;
        bool EOB_found = false;
        
        if (!progressive_AC_refine) { // AC First Scans
            if (EOB_run > 0) {
                // If a previous call created a nonzero EOB run then we decrement the
                // counter and return.
                --EOB_run;
            }
            else {
                for (coeff_counter = zigZagStart; coeff_counter <= zigZagEnd;) {
                    int codeLength = 0;
                    for (codeLength = 1; codeLength <= 16; ++codeLength) {
                        // found in huffman table
                        found = false;
                        
                        // Keep grabbing one bit at a time till we find one thats a huffman code
                        int code = LookNBits(codeLength);
                        
                        // Check if its one of our huffman codes
                        // Current Huffman table (AC table)
                        if (is_exist_in_huffman_codes(code, codeLength, currentComponent, decodedValue, false))
                        {
                            // Skip over k bits, since we found the huffman value
                            // and looked ahead before
                            SkipNBits(codeLength);
                            // set the found in huffman table to true
                            found = true;
                            
                            // Let's separate byte to two nibbles
                            int valCode = decodedValue;
                            uint_8 count_0 = valCode >> 4;	// Number RunLengthZeros
                            uint_8 size_val = valCode & 0xF;	// Number of bits for our data
                            
                            if (size_val == 0) {
                                // A zero value ssss with rrrr == 15 means to skip
                                // 16 zero coefficients.
                                if (count_0 == 0xF) {
                                    coeff_counter += 16;  // skip 16 zeros
                                }
                                else {
                                    // A zero value ssss with rrrr != 15 means to create
                                    // End of Band run.
                                    
                                    // The EOB run includes the current block. This is why we
                                    // do no processing for rrrr = 0 and substract one when
                                    // rrrr != 0.
                                    if (count_0 != 0) {
                                        int length_temp = count_0;
                                        int bits = GetNBits(length_temp);
                                        EOB_run = (1 << count_0) + bits - 1;
                                    }
                                    EOB_found = true;
                                    break;
                                }
                                
                            } // end if (size_val == 0)
                            else
                            {
                                // When ssss != 0, rrrr gives the number of zero elements to skip
                                // before the next non-zero coefficient.
                                coeff_counter += count_0; //skip count_0 zeroes
                                // if the coeffs counter is greater than normally 63 coeffs
                                if (coeff_counter > 63) {
                                    cout << ftell(fp) << endl;
                                    cout << "-|- ##ERROR## Coefficients counter = " << coeff_counter << " is greater than ACcount " << ACcount << endl;
                                    // in case of error, doing the other stuff will just do more errors so return here
                                    
                                    // Log it
                                    // logErrorPictures();
                                    return;
                                } // end if
                                
                                short ac_coeff = GetNBits(size_val);
                                ac_coeff = DetermineSign(ac_coeff, size_val);
                                DCT_tcoeff[coeff_counter++] = (ac_coeff << (int) approximationL);
                                
#if DEBUGLEVEL > 50
                                //cout << "Decoded value, run_len, size_val: " << (int) valCode << "--" << (int) count_0 << "--" << (int)size_val << endl;
                                //cout  << "AC coeff: " << ac_coeff << endl;
                                //cout  << "Code: " << (int)code << " -- Length: " << (int)codeLength << endl;
                                if (ftell(fp) >= 32936) {
                                    cout << std::dec << "Postion in Byte: " << ftell(fp) << endl;
                                    cout << "Code: " << code << " -- Length: " << codeLength << endl;
                                    printf("Code: %X \n", code);
                                    cout << "Component:" << currentComponent << endl;
                                    cout << "Decoded value, run_len, size_val: " << (int)valCode << "--" << (int)count_0 << "--" << (int)size_val << endl;
                                    cout << "Decoded value in hex: " << std::hex << valCode << endl;
                                    cout << std::dec << "AC coeff: " << ac_coeff << endl;
                                    cout << "postion X:" << currentX << ", position Y:" << currentY << endl;
                                    cout << std::hex << "reservoir @ end: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
                                }
#endif
                            } // end else
                            // Found so we can exit out
                            break;
                        }// end if is_exist_in_huffman_codes
                    } // end for codeLengths
                    
                    if (!found) {
                        cout << ftell(fp) << endl;
                        cout << currentX << "  " << currentY;
                        cout << "-|- ##ERROR## (AC case) Code value not found in Huffman table: " << endl;
                        // Log it
                        // logErrorPictures();
                        return;
                    }
                    
                    if (codeLength > 16) {
                        coeff_counter++;
                    }
                    if (EOB_found) break;
                }
            }
        }
        else { //AC Refine
            //    Section G.1.2.3 defines how to encode refining scans for AC
            //    coefficients. Unfortunately this section is vague and
            //    undecipherable. Reversing an undecipherable process results
            //    in something unimaginable. This is a "best-guess" interpretation
            //    that seems to work.
            //
            //    The basic process at work is that zero counts do not include nonzero
            //    values. Whenever we skip a value due to zero count or End of Band runs
            //    we have to read one bit to refine each non-zero value we skip. The
            //    process is ugly and it means that data is encoding out of order.
            // Retrieve Saved DCT value for current block
            for (int i = 0; i < 64; ++i)
            {
                switch (currentComponent) {
                    case COMPONENT_Y:
                        DCT_tcoeff[i] = tCoeff_Y [currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                        break;
                    case COMPONENT_Cb:
                        DCT_tcoeff[i] = tCoeff_Cb[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                        break;
                    case COMPONENT_Cr:
                        DCT_tcoeff[i] = tCoeff_Cr[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                        break;
                    default:
                        DCT_tcoeff[i] = tCoeff_Y [currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                        break;
                } // end switch
                
            }
            
            for (coeff_counter = zigZagStart; coeff_counter <= zigZagEnd;) {
                if (EOB_run != 0) {
                    // An EOB run has caused us to skip entire data units. We need
                    // to refine any previously non-zero coefficients.
                    // Notice that we do not initialize kk here. We could be using
                    // an EOB run to skip all the remaining coefficients in the current one.
                    for (; coeff_counter <= zigZagEnd; ++coeff_counter) {
                        // REFINE AC COEFFICIENT
                        if (DCT_tcoeff[coeff_counter] != 0) {
                            if (DCT_tcoeff[coeff_counter] > 0) {
                                if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (1 << approximationL);
                            }
                            else if (DCT_tcoeff[coeff_counter] < 0) {
                                if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (-1 << approximationL);
                            }
                        }
                    }
                    --EOB_run;
                }
                else {
                    int codeLength = 0;
                    for (codeLength = 1; codeLength <= 16; ++codeLength) {
                        // found in huffman table
                        found = false;
                        
                        // Keep grabbing one bit at a time till we find one thats a huffman code
                        int code = LookNBits(codeLength);
                        
                        // Check if its one of our huffman codes
                        // Current Huffman table (AC table)
                        if (is_exist_in_huffman_codes(code, codeLength, currentComponent, decodedValue, false))
                        {
                            // Skip over k bits, since we found the huffman value
                            // and looked ahead before
                            /*if (ftell(fp) >= 1930) {
                             cout << std::hex << "reservoir: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
                             }*/
                            SkipNBits(codeLength);
                            // set the found in huffman table to true
                            found = true;
                            
                            // Let's separate byte to two nibbles
                            int valCode = decodedValue;
                            uint_8 count_0 = valCode >> 4;	// Number RunLengthZeros
                            uint_8 size_val = valCode & 0xF;	// Number of bits for our data
                            
                            if (size_val == 0) {
                                // ssss == 0 means that we either have an EOB run or we need to
                                // 16 non-zero coefficients.
                                if (count_0 == 0xF) {
                                    //coeff_counter += 16;  // skip 16 zeros
                                    for (unsigned int ii = 0; coeff_counter <= zigZagEnd && ii < 16; ++coeff_counter) {
                                        if (coeff_counter > zigZagEnd) cout << " - | -##ERROR## Coefficient Counter out of range @ AC Refine" << endl;
                                        // REFINE AC COEFFICIENT
                                        if (DCT_tcoeff[coeff_counter] != 0) {
                                            if (DCT_tcoeff[coeff_counter] > 0) {
                                                if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (1 << approximationL);
                                            }
                                            else if (DCT_tcoeff[coeff_counter] < 0) {
                                                if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (-1 << approximationL);
                                            }
                                        }
                                        else {
                                            ++ii;
                                        }
                                    }
                                }
                                else {
                                    if (count_0 == 0) {
                                        EOB_run = 1;
                                    }
                                    else {
                                        int length_temp = count_0;
                                        int bits = GetNBits(length_temp);
                                        EOB_run = (1 << count_0) + bits;
                                    }
                                }
                            } // end if (size_val == 0)
                            else if (size_val == 1)
                            {
                                // ssss == 1 means that we are creating a new non-zero
                                // coefficient. rrrr gives the number of zero coefficients to
                                // skip before we reach this one.
                                // Save the value for the new coefficient. Unfortunately the data
                                // is stored out of order.
                                int length_temp = 1;
                                int newvalue = GetNBits(length_temp);
                                
                                for (unsigned int zerocount = 0; coeff_counter <= 63 && (zerocount < count_0 || DCT_tcoeff[coeff_counter] != 0); ++coeff_counter) {
                                    
                                    if (coeff_counter > zigZagEnd) {
                                        cout << "-|- ##ERROR## IN AC REFINE 1: Coefficients counter = " << coeff_counter << " is greater than ACcount " << ACcount << endl;
                                        return;
                                    } // end if
                                    
                                    // REFINE AC COEFFICIENT
                                    if (DCT_tcoeff[coeff_counter] != 0) {
                                        if (DCT_tcoeff[coeff_counter] > 0) {
                                            if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (1 << approximationL);
                                        }
                                        else if (DCT_tcoeff[coeff_counter] < 0) {
                                            if (GetNBits(1) != 0) DCT_tcoeff[coeff_counter] += (-1 << approximationL);
                                        }
                                    }
                                    else {
                                        ++zerocount;
                                    }
                                }
                                
                                if (coeff_counter > zigZagEnd) {
                                    cout << ftell(fp) << endl;
                                    cout << "-|- ##ERROR## IN AC REFINE 2: Coefficients counter = " << coeff_counter << " is greater than ACcount " << ACcount << endl;
                                    return;
                                } // end if
                                
                                if (newvalue) {
                                    DCT_tcoeff[coeff_counter++] = (1 << (int)approximationL);
                                }
                                else {
                                    DCT_tcoeff[coeff_counter++] = (-1 << (int)approximationL);
                                }
                                //short ac_coeff = GetNBits(size_val);
                                //ac_coeff = DetermineSign(ac_coeff, size_val);
                                //DCT_tcoeff[coeff_counter++] = (ac_coeff << approximationL);
                                
#if DEBUGLEVEL > 50
                                //cout << "Decoded value, run_len, size_val: " << (int) valCode << "--" << (int) count_0 << "--" << (int)size_val << endl;
                                //cout  << "AC coeff: " << ac_coeff << endl;
                                //cout  << "Code: " << (int)code << " -- Length: " << (int)codeLength << endl;
                                if (ftell(fp) >= 32936) {
                                    cout << std::dec << "Postion in Byte: " << ftell(fp) << endl;
                                    cout << "Code: " << code << " -- Length: " << codeLength << endl;
                                    printf("Code: %X \n", code);
                                    cout << "Component:" << currentComponent << endl;
                                    cout << "Decoded value, run_len, size_val: " << (int)valCode << "--" << (int)count_0 << "--" << (int)size_val << endl;
                                    cout << "Decoded value in hex: " << std::hex << valCode << endl;
                                    cout << std::dec << "AC coeff: " << ac_coeff << endl;
                                    cout << "postion X:" << currentX << ", position Y:" << currentY << endl;
                                    cout << std::hex << "reservoir @ end: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
                                }
#endif
                            } // end else
                            // Found so we can exit out
                            break;
                        }// end if is_exist_in_huffman_codes
                    } // end for codeLengths
                    
                    if (!found) {
                        cout << ftell(fp) << endl;
                        cout << currentX << "  " << currentY;
                        cout << "-|- ##ERROR## (AC case) Code value not found in Huffman table: " << endl;
                        // Log it
                        // logErrorPictures();
                        return;
                    }
                    
                    if (codeLength > 16) {
                        coeff_counter++;
                    }
                    if (EOB_found) break;
                }
            }
        }
        
        // We've decoded a block of data, so copy it across to our buffer
        Component comp = components.at(currentComponent);
        for (int j = zigZagStart; j <= zigZagEnd; ++j) {
            comp.m_DCT[j] = DCT_tcoeff[j];
        }
        
#if DEBUGLEVEL > 30
        if (currentX == 224 && currentY == 0) {
            dumpDCTValues(comp.m_DCT);
        }
        
#endif
    }
} // end process_huffmann_data_unit_progressive

// Process Huffman Data in Sequential Mode
void jpeg_decoder::process_huffmann_data_unit(int currentComponent, int currentX, int currentY) {
    // We memset it here, as later on we can just skip along, when we have lots
    // of leading zeros, for our AC run length encoding :)
    short DCT_tcoeff[64];
    memset(DCT_tcoeff, 0, sizeof(DCT_tcoeff)); //Initialize DCT_tcoeff
    // How many AC elements should we read?
    int ACcount = zigZagEnd - zigZagStart;
    
    // found in huffman table
    bool found = false;
    
    // decodedValue from huffmanTable
    int decodedValue = 0;
    
#if DEBUGLEVEL > 40
    printDCHuffmanCodes(currentComponent);
    printACHuffmanCodes(currentComponent);
#endif
    
    // First thing is get the 1 DC coefficient at the start of our 64 element
    // block (the length of the code word is maximum 16 bits or 2 bytes)
    for (int codeLength = 1; codeLength <= 16; ++codeLength)
    {
        // Keep grabbing one bit at a time till we find one thats a huffman code
        int code = LookNBits(codeLength);
        
        // found in huffman table
        found = false;
        
        // Check if its one of our huffman codes
        // Current Huffman table (Initially DC)
        if (is_exist_in_huffman_codes(code, codeLength, currentComponent, decodedValue))
        {
            /*if (ftell(fp) >= 1930) {
             cout << std::hex << "reservoir: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
             }*/
            // you looked ahead, skip those bits
            SkipNBits(codeLength);
            
            // set the found in huffman table to true
            found = true;
            
            // The decoded value is the number of bits we have to read in next
            int numDataBits = decodedValue;
            
#if DEBUGLEVEL > 50
            cout << "The decoded value is the number of bits we have to read in next equals to " << numDataBits << endl;
#endif
            // We know the next k bits are for the actual data
            if (numDataBits == 0)
            {
                DCT_tcoeff[0] = previousDC[currentComponent];
                
#if DEBUGLEVEL > 50
                cout << "First DC coefficient is " << DCT_tcoeff[0] << endl;
#endif
            }
            else
            {
                // residual DC:
                short residualDC = GetNBits(numDataBits);
#if DEBUGLEVEL > 50
                cout << "Delta of DC coefficient " << residualDC << endl;
#endif
                
                residualDC = DetermineSign(residualDC, numDataBits);
                
#if DEBUGLEVEL > 50
                cout << "Delta of DC coefficient with sign " << residualDC << endl;
#endif
                
                DCT_tcoeff[0] = residualDC + previousDC[currentComponent];
                previousDC[currentComponent] = DCT_tcoeff[0];
                
#if DEBUGLEVEL > 50
                //cout << "Component: " << currentComponent  << ", DC coefficient is: " << DCT_tcoeff[0] << " and prev DC coeff is updated" << endl;
                if (ftell(fp) >= 1930) {
                    cout << "DC Postion in Byte: " << ftell(fp) << endl;
                    cout << "Code: " << code << " -- Length: " << codeLength << endl;
                    printf("Code: %X \n", code);
                    cout << "Component:" << currentComponent << endl;
                    cout << "Decoded value in hex: " << std::hex << numDataBits << endl;
                    cout << std::dec << "DC coeff: " << DCT_tcoeff[0] << endl;
                    cout << "residual DC: " << residualDC << endl;
                    cout << "postion X:" << currentX << ", position Y:" << currentY << endl;
                }
#endif
            } // end else
            // Found so we can exit out
            break;
        } // end if is_exist
    } // end for loop on codeLengths
    
    if (!found){
        cout << "-|- ##ERROR## (DC case) Code value not found in Huffman table: " << endl;
        // Log it
        // logErrorPictures();
        return;
    }
    
    // No AC coefficients required?
    if (ACcount == 0 || losslessFormat) {
        cout << "-|- ##ERROR## (ACcount == 0 || losslessFormat) No AC coefficients reading is required! " << endl;
        return;
    }
    
    // Second, the 63 AC coefficient
    int coeff_counter = 1;
    bool EOB_found = false;
    do {
        
        int codeLength = 0;
        for (codeLength = 1; codeLength <= 16; ++codeLength) {
            
            // found in huffman table
            found = false;
            
            // Keep grabbing one bit at a time till we find one thats a huffman code
            int code = LookNBits(codeLength);
            
            // Check if its one of our huffman codes
            // Current Huffman table (AC table)
            if (is_exist_in_huffman_codes(code, codeLength, currentComponent, decodedValue, false))
            {
                // Skip over k bits, since we found the huffman value
                // and looked ahead before
                /*if (ftell(fp) >= 1930) {
                 cout << std::hex << "reservoir: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
                 }*/
                SkipNBits(codeLength);
                // set the found in huffman table to true
                found = true;
                
                /* If AC element is 0xAB for example, then we have to separate it in two nibbles
                 * First nible (repeating RLE) is RRRR bits, second is SSSS bits
                 * RRRR bits defines the #zero elements are before this element
                 * SSSS bits defines the #binary digits our AC element has (if 1001
                 * then we have to read next 9 elements from file)
                 */
                
                // Let's separate byte to two nibbles
                int valCode = decodedValue;
                uint_8 count_0 = valCode >> 4;	// Number RunLengthZeros
                uint_8 size_val = valCode & 0xF;	// Number of bits for our data
                
                /*cout << "Code: " << code << " -- Length: " << codeLength << endl;
                 printf("Code: %X \n", code);*/
                
                if (size_val == 0){
                    // RLE
                    if (count_0 == 0) {
                        EOB_found = true;	// EOB found, go out
                    }
                    else if (count_0 == 0xF) {
                        coeff_counter += 16;  // skip 16 zeros
                    }
                    
                } // end if (size_val == 0)
                else
                {
                    
                    coeff_counter += count_0; //skip count_0 zeroes
                    
                    // if the coeffs counter is greater than normally 63 coeffs
                    if (coeff_counter > ACcount)
                    {
                        cout << "-|- ##ERROR## Coefficients counter = " << coeff_counter << " is greater than ACcount " << ACcount << endl;
                        // in case of error, doing the other stuff will just do more errors so return here
                        
                        // Log it
                        // logErrorPictures();
                        
                        return;
                    } // end if
                    
                    short ac_coeff = GetNBits(size_val);
                    ac_coeff = DetermineSign(ac_coeff, size_val);
                    DCT_tcoeff[coeff_counter++] = ac_coeff;
#if DEBUGLEVEL > 50
                    //cout << "Decoded value, run_len, size_val: " << (int) valCode << "--" << (int) count_0 << "--" << (int)size_val << endl;
                    //cout  << "AC coeff: " << ac_coeff << endl;
                    //cout  << "Code: " << (int)code << " -- Length: " << (int)codeLength << endl;
                    if (ftell(fp) >= 1930) {
                        cout << std::dec << "Postion in Byte: " << ftell(fp) << endl;
                        cout << "Code: " << code << " -- Length: " << codeLength << endl;
                        printf("Code: %X \n", code);
                        cout << "Component:" << currentComponent << endl;
                        cout << "Decoded value, run_len, size_val: " << (int)valCode << "--" << (int)count_0 << "--" << (int)size_val << endl;
                        cout << "Decoded value in hex: " << std::hex << valCode << endl;
                        cout << std::dec << "AC coeff: " << ac_coeff << endl;
                        cout << "postion X:" << currentX << ", position Y:" << currentY << endl;
                        cout << std::hex << "reservoir @ end: " << g_reservoir << "--" << g_nbits_in_reservoir << endl;
                    }
#endif
                } // end else
                
                // Found so we can exit out
                break;
            }// end if is_exist_in_huffman_codes
        } // end for codeLengths
        
        if (!found){
            cout << "-|- ##ERROR## (AC case) Code value not found in Huffman table: " << endl;
            // Log it
            // logErrorPictures();
            return;
        }
        
        if (codeLength > 16)
        {
            coeff_counter++;
        }
        
    } while ((coeff_counter <= 63) && (!EOB_found));
    
    // We've decoded a block of data, so copy it across to our buffer
    Component comp = components.at(currentComponent);
    for (int j = 0; j < 64; ++j)
    {
        comp.m_DCT[j] = DCT_tcoeff[j];
    }
    
#if DEBUGLEVEL > 30
    if (currentX == 224 && currentY == 0) {
        dumpDCTValues(comp.m_DCT);
    }
    
#endif
    
} // end process_huffmann_data_unit

void jpeg_decoder::decode_single_block(int offset, int stride, int currentComponent, uint_8 *outputBuf, int currentX, int currentY) {
    
    // fetch the DCT coefficients for the component
    short* input_dct_coeffs = components.at(currentComponent).m_DCT;
    
    // Create a temp 8x8, i.e. 64 array for the data
    int data[64] = { 0 };
    
    // Progressive Only: Retrieve stored DCT value from the t_Coeff
    if (progressive_Huff_Format) {
        for (int i = 0; i < 64; ++i)
        {
            switch (currentComponent) {
                case COMPONENT_Y:
                    data[i] = tCoeff_Y[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                    break;
                case COMPONENT_Cb:
                    data[i] = tCoeff_Cb[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                    break;
                case COMPONENT_Cr:
                    data[i] = tCoeff_Cr[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                    break;
                default:
                    data[i] = tCoeff_Y[currentY_dct[currentComponent] + floor(i / 8)][currentX_dct[currentComponent] + (i % 8)];
                    break;
            } // end switch
        }
    }
    
    // Copy our data into the temp array
    for (int i = zigZagStart; i <= zigZagEnd; ++i)
    {
        data[i] = input_dct_coeffs[i];
 
	}
    
#if DEBUGLEVEL > 14
    cout << "Block init 1: \n" << endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            
            int c = j + 8 * i;
            cout << data[c] << " ";
        }
        
        cout << "\n";
    }
#endif
    
    // Create a reshaped matrix of size 8x8
    int dataReshapedInto8x8[8][8] = { 0 };
    reshapeArrayInto8x8(dataReshapedInto8x8, data);
	
#if DEBUGLEVEL > 20
    if (currentX == 1008 && currentY == 0 && currentComponent == COMPONENT_Cb) {
        cout << "after reshapeArrayInto8x8 " << endl;
        //    dumpDecodedBlock(dataReshapedInto8x8);
        dumpDecodedBlockInDecimal(dataReshapedInto8x8);
    }
#endif
    
#if DEBUGLEVEL > 50
    //    // TODO: remove FANKOOSH
    static int fankoosh = 0;
    if(currentComponent == COMPONENT_Y) {
        fankoosh++;
        ofstream myfile;
        std::string path_to_files = "C:/Users/y77jiang/OneDrive - University of Waterloo/5e. TCM-Inception C++/jpeg_tcm/JPEG/";
        std::string output_csv_name = path_to_files + "goose_Y_dec.csv";
        myfile.open (output_csv_name, std::ofstream::out | std::ofstream::app);
        
        std::stringstream oss;
        std::size_t found = output_csv_name.find_last_of(".");
        std::string path_with_name = output_csv_name.substr(0, found);
        found = output_csv_name.find_last_of("/\\");
        std::string name_file_only = path_with_name.substr(found+1);
        
        myfile << currentX << "-" << currentY << "\n";
        
        if (currentX == 1008 && currentY == 0) {
            cout << currentX << "-" << currentY << "\n";
        }
        
        for(int i = 0; i < 8; ++i) {
            for(int j = 0; j < 8; ++j) {
                myfile << dataReshapedInto8x8[i][j] << ",";
                
                if (currentX == 1008 && currentY == 0) {
                    cout << dataReshapedInto8x8[i][j] << ",";
                }
            }
            
            if (currentX >= 1008 && currentY == 0) {
                cout << "\n";
            }
            
            if( (i + 1) < jpegImageHeight){
                myfile << "\n";
            }
        }
        myfile.close();
    }
#endif
    
#if DEBUGLEVEL > 50
    //    // TODO: remove FANKOOSH
    static int fankoosh = 0;
    if (currentComponent == COMPONENT_Y) {
        fankoosh++;
        ofstream myfile;
        std::string path_to_files = "C:/Users/y77jiang/OneDrive - University of Waterloo/5e. TCM-Inception C++/jpeg_tcm/dataset/";
        std::string output_csv_name = path_to_files + "goose_Y_dec.csv";
        myfile.open(output_csv_name, std::ofstream::out | std::ofstream::app);
        
        std::stringstream oss;
        std::size_t found = output_csv_name.find_last_of(".");
        std::string path_with_name = output_csv_name.substr(0, found);
        found = output_csv_name.find_last_of("/\\");
        std::string name_file_only = path_with_name.substr(found + 1);
        
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                if(!(i==0 && j ==0))
                    myfile << dataReshapedInto8x8[j][i] << ",";
            }
        }
        
        myfile << "\n";
        myfile.close();
    }
#endif
    
    // If sequential, process firsly; if progressive, save and process later
    if (!progressive_Huff_Format) {
        // Perform dequantization
        multiplyWithQuantizationTable(dataReshapedInto8x8, currentComponent);
        
        // Perform inverse zig zag scanning
        int dataReshapedInto8x8Zig[8][8] = { 0 };
        inverseZigZagScanning(dataReshapedInto8x8Zig, dataReshapedInto8x8);
		
#if DEBUGLEVEL > 20
        if (currentX == 1008 && currentY == 0 && currentComponent == COMPONENT_Cb) {
            cout << "Before IDCT and after quantization in decimal representation: " << endl;
            dumpDecodedBlockInDecimal(dataReshapedInto8x8Zig);
        }
#endif
        // Add the raw DCT coefficients after dequantization after inverse zigzag scanning:
        addtCoeffBlock(dataReshapedInto8x8Zig, currentComponent, currentX, currentY);
        
        // Perform IDCT
        int output_value_dct[8][8] = { 0 };
        perform_idct(output_value_dct, dataReshapedInto8x8Zig);
 
#if DEBUGLEVEL > 20
        cout << "After IDCT and after quantization in decimal representation: " << endl;
        dumpDecodedBlockInDecimal(output_value_dct);
#endif
        
        // Note: it does not matter the component has subsampling or not
        // addBlockSubsampling will take care of it - the method is general
        addBlockSubsampling(output_value_dct, currentComponent);
		

        // Perform transpose (your block is transposed compared to the reference
        // code)
        int output_value[8][8] = { 0 };
        transpose(output_value, output_value_dct);
        
        // Level Shift each element (i.e. add 128), and copy to our
        // output
        unsigned char *outptr = outputBuf;
        for (int y = 0; y < 8; ++y){
            for (int x = 0; x < 8; ++x)
            {
                // TODO: use switch instead of if
                // Only 8 and 12-bit is supported by DCT per JPEG standard
                if (jpegImageSamplePrecision == 8) {
                    output_value[x][y] += 128;
                }
                else if (jpegImageSamplePrecision == 12) {
                    output_value[x][y] += 2048;
                }
                
                outptr[x] = Clamp(output_value[x][y]);
            } // end inner loop
            outptr += stride;
        } // end outer loop
    }
    else {
        // Progressive: ADD tcoeff block in zigzaged form, and quantization/idct done later
        addtCoeffBlock(dataReshapedInto8x8, currentComponent, currentX, currentY);
    }
#if DEBUGLEVEL > 60
    dumpDecodedBlockInDecimal(output_value);
#endif
    
#if DEBUGLEVEL > 40
    printf("# Decoded 8x8 Block#\n");
    outptr = outputBuf;
    for(int y = 0; y < 8; ++y)
    {
        for(int x = 0; x < 8; ++x)
        {
            
            printf("%d ", outptr[x]);
        }
        
        outptr += stride;
        cout << "\n";
    }
#endif
} // end decode_single_block

float jpeg_decoder::C(int u)
{
    if (u == 0)
        return (1.0f/sqrtf(2));
    else
        return 1.0f;
} // end C

int jpeg_decoder::func(int x, int y, const int block[8][8])
{
    const float PI = 3.14f;
    static bool isFirstTime = true;
    if(isFirstTime)
    {
        const double inv16 = 1.0 / 16.0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                cosine_idct[j][i] = cosf( (2.0 * i + 1) * j * PI * inv16 );
            }
        }
        isFirstTime = false;
    }
    
  
    float sum = 0;
    for( int u = 0; u < 8; ++u)
    {
        for(int v = 0; v < 8 ; ++v)
        {
            sum += ( C(u) * C(v) ) * block[u][v] * cosine_idct[u][x]  * cosine_idct[v][y];
        } // end inner loop
    } // end outer loop
    return (int) ((1.0/4.0) * sum);
} // end func

void jpeg_decoder::perform_idct(int outBlock[8][8], const int inBlock[8][8])
{
    for(int y = 0; y < 8; ++y)
    {
        for(int x = 0; x < 8; ++x)
        {
            outBlock[x][y]  =  func( x, y, inBlock);
        } // end inner loop
    } // end outer loop
} // end perform_idct


void jpeg_decoder::transpose(int outBlock[8][8], const int inBlock[8][8])
{
    int temp [8][8];
    
    // copy the inBlock over to temp
    for(int y=0; y<8; y++)
    {
        for(int x=0; x<8; x++)
        {
            temp[x][y] = inBlock[x][y];
        }
    }
    // outBlock transpose
    for(int y=0; y<8; y++)
    {
        for(int x=0; x<8; x++)
        {
            outBlock[x][y]  =  temp[y][x];
        } // end inner loop
    } // end outer loop
} // end transpose


void jpeg_decoder::reshapeArrayInto8x8(int outArray[8][8], const int inArray[64]) {
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            int c = j + 8*i;
            outArray[i][j] = inArray[c];
        } // end inner loop
    } // end outer loop
} // end reshapeArrayInto8x8

// check exist in huffman table
bool jpeg_decoder::is_exist_in_huffman_codes(int code, int codeLength, int currentComponent, int& decodedValue, bool is_dc) {
    
    // pointer to a huffman table
    HuffmanTable* hTable;
    
    // is DC coefficient?
    if(is_dc) {
        hTable = componentTablesDC[currentComponent];
    }
    else {
        hTable = componentTablesAC[currentComponent];
    }
    
    std::map<huffKey, uint_8>::iterator iter;
    for (iter = hTable->huffData.begin();iter != hTable->huffData.end(); ++iter) {
        
        int codeInTable = iter->first.second;
        int codeLengthInTable = iter->first.first;
        
        if(code == codeInTable && codeLength == codeLengthInTable)
        {
            decodedValue = iter->second;
#if DEBUGLEVEL>40
            printf("(Yaay) Found this code: %d, with codeLength: %d, and value = %02X\n",
                   code, codeLength, iter->second);
#endif
            return true;
        }
        
#if DEBUGLEVEL>20
        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X or [%s] \n",
               iter->first.second, iter->first.first,
               iter->second, IntToBinary(iter->first.second, iter->first.first));
#endif
    }
    
    // not found
#if DEBUGLEVEL>20
    cout << "Not found this code: " << code << ", with codeLength: " << codeLength << endl;
#endif
    return false;
    
} // end is_exist_in_huffman_codes

// Color conversion functions:
void jpeg_decoder::convert_yccb_to_rgb24(int y, int cb, int cr, int& r, int& g, int& b) {
    float red, green, blue;
    
    red   = y + 1.402f * (cb-128) ;
    green = y - 0.34414f * (cr - 128) - 0.71414f *(cb - 128);
    blue  = y + 1.772f *( cr - 128);
    
    r = static_cast<int>(Clamp(static_cast<int>(red)));
    g = static_cast<int>(Clamp(static_cast<int>(green)));
    b = static_cast<int>(Clamp(static_cast<int>(blue)));
    
} // end convert_yccb_to_rgb24

void jpeg_decoder::ycrcb_to_rgb24_image() {
    
    // <red, green, blue> pixel
    int r, g, b;
    constexpr int nRGBComponents = 3;
    
    for (int y = 0; y < upscale_height; ++y) {
        for (int x = 0; x < upscale_width; ++x) {
            int pixel_offset = x * nRGBComponents + jpegImageWidth * nRGBComponents * y;
            
            // Fetch the Y'CbrCr pixel
            int yc = m_YPicture_buffer[y][x]; // y is the column, x is the row
            
            // TODO: Support for grayscale components with various sample precision
            // Note: Cb, Cr should be 128 if there is only one component
            // and it is a BMP format
            int cb = 128;
            int cr = 128;
            if(components.size() >= (COMPONENT_Cb+1)) {
                cb = m_CbPicture_buffer[y][x];
            }
            
            if(components.size() >= (COMPONENT_Cr+1)) {
                cr = m_CrPicture_buffer[y][x];
            }
            
            // Convert the Y'CbrCr pixel into RGB pixel
            convert_yccb_to_rgb24(yc, cr, cb, r, g, b);
            
            // Set the pixel to its RGB <pattern is: RGB RGB RGB...etc>
            m_rgb[pixel_offset]     = Clamp(r);
            m_rgb[pixel_offset + 1] = Clamp(g);
            m_rgb[pixel_offset + 2] = Clamp(b);
        } // end inner
    } // end outer
    
} // end ycrcb_to_rgb24_image

void jpeg_decoder::multiplyWithQuantizationTable(int dataBlock[8][8], int currentComponent) {
    
    
#if DEBUGLEVEL > 399
    dumpQuantizationTable(currentComponent);
#endif
    
    // Note: quantization table needs to be dezig-zagged in order
    // to be able to do proper dequantization
    QuantizationTable* table = this->components[currentComponent].componentQuantizationTable;
    
    int x_indices[64] = {0, 1, 0, 0, 1, 2, 3, 2,
        1, 0, 0, 1, 2, 3, 4, 5,
        4, 3, 2, 1, 0, 0, 1, 2,
        3, 4, 5, 6, 7, 6, 5, 4,
        3, 2, 1, 0, 1, 2, 3, 4,
        5, 6, 7, 7, 6, 5, 4, 3,
        2, 3, 4, 5, 6, 7, 7, 6,
        5, 4, 5, 6, 7, 7, 6, 7};
    
    int y_indices[64] = {0, 0, 1, 2, 1, 0, 0, 1,
        2, 3, 4, 3, 2, 1, 0, 0,
        1, 2, 3, 4, 5, 6, 5, 4,
        3, 2, 1, 0, 0, 1, 2, 3,
        4, 5, 6, 7, 7, 6, 5, 4,
        3, 2, 1, 2, 3, 4, 5, 6,
        7, 7, 6, 5, 4, 3, 4, 5,
        6, 7, 7, 6, 5, 6, 7, 7};
    
    
#if DEBUGLEVEL > 30
    cout << "Block:--------------" << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            cout << std::dec <<  dataBlock[i][j] << " ";
        }
        
        cout << "\n";
    }
    
    cout << "QTable:--------------" << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            
            cout <<  table->quantizationTableData[i][j] << " ";
        }
        
        cout << "\n";
    }
#endif
    
    
#if DEBUGLEVEL > 30
    cout << "Print Quantization Operation:--------------" << endl;
    cout << "A) Block: " << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            int val = dataBlock[i][j];
            cout << "[ " << std::dec <<  val << " ] ";
        }
        cout << "\n";
    }
    
    cout << "B) Operation: " << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            int idx = j + i*8;
            int val = dataBlock[i][j];
            int q = table->quantizationTableData[y_indices[idx]][x_indices[idx]];
            cout << "[ " << std::dec <<  val << " * " << q << " ] ";
        }
        cout << "\n";
    }
#endif
    
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            int idx = j + i*8;
            dataBlock[i][j] = dataBlock[i][j] * table->quantizationTableData[y_indices[idx]][x_indices[idx]];
        } // end inner
    } // end outer
    
    
    
#if DEBUGLEVEL > 30
    cout << "Block After dequant:--------------" << endl;
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            cout << std::dec << dataBlock[i][j] << " ";
        }
        
        cout << "\n";
    }
#endif
    
} // end multiplyWithQuantizationTable

/* This method will be called repeatedly to place 8x8 blocks on the luminance picture, from
 left to right and then onto next line.
 */
/* This method is much more complex, having to deal with subsampling.
 Also somewhat slower due to all the counters. */
void jpeg_decoder::addBlockSubsampling(int dataBlock[8][8], int comp /* Current component */ ) {
    
    // Avoid function calls... these functions are inline but the library is external...
    int width = jpegImageWidth, height = jpegImageHeight;
    // pixel jump is used for storing YCbYCr YCbCr...etc; will leave it, but I will only store Y
    int pixelJump = static_cast<int>(components.size());
    // TODO: adding the code for 12 bits, but it is not supported
    if (jpegImageSamplePrecision == 12) {
        pixelJump *= 2;
    }
    
    // How many Y's in the horizontal and vertical direction (2x2 is the usual case)
    int HFactor = components[comp].HFactor, VFactor = components[comp].VFactor;
    
    // These Y's should be scaled/repeated how many times horizonatally and vertically
    int HScale = components[comp].HScale,   VScale = components[comp].VScale;
    
    // raster block sanning/traverse using x, y positions:
    for (int y = 0; y < 8; ++y) {
        
        if (currentX[comp] >= width) break;
        if (currentY[comp] + y >= height) break;
        
        // Repeat each line VScale times
        for (int vfy = 0; vfy < VScale; ++vfy) {
            
            int picture_y = currentY[comp] + y * VScale + vfy;
            if (picture_y >= height) break;
            
            for (int x = 0; x < 8 ; ++x) {
                
                if (currentX[comp] + x * HScale >= width) break;
                
                // x and y are reverted compared to image
                int value = dataBlock[y][x];
                
                // save the initial raw value
                int raw_value = value;
                
                
                // level-shift value because it is not initially level shifted
                if (jpegImageSamplePrecision == 8) {
                    value += 128;
                }
                else if (jpegImageSamplePrecision == 12) {
                    value += 2048;
                }
                
                // Clamp the value:
                value = Clamp(value);
                
                // Repeat each pixel HScale times
                int realx = currentX[comp] + x * HScale;
                for (int i = 0; i < HScale; ++i) {
                    
                    // Set the pixel <Xcor: realX, Ycor: picture_y>
                    
                    switch(comp) {
                        case COMPONENT_Y:
							
                            m_YPicture_buffer[picture_y][realx] = value;
                            break;
                        case COMPONENT_Cb:
                            m_CbPicture_buffer[picture_y][realx] = value;
                            break;
                        case COMPONENT_Cr:
                            m_CrPicture_buffer[picture_y][realx] = value;
                            break;
                        default:
                            m_YPicture_buffer[picture_y][realx] = value;
                            cout << "-|- ##ERROR## You are sending the wrong component - will assume it's luminance. Most probably will fail!" << endl;
                            
                            break;
                    } // end switch
                    
                    if (++realx >= width) break;
                } // end for i
            } // end for x
        } // end for vfy
    } // end for y
    
    
    // currentX and currentY are now 8 pixels below values when function is called
    // they should be 8 pixels to the right multiplied by HScale
    // unless this is image edge, in which case they should be in the next line
    // Update starting X and Y for next block, taking into account subsampling
    currentX[comp] += 8 * HScale;
    currentBlockHFactor[comp]++;
    
    // you made a line of blocks
    if(currentBlockHFactor[comp] >= HFactor) {
        
        // restore the current X to its initial position and reset the counters
        currentX[comp] -= 8 * HScale * HFactor;
        currentBlockHFactor[comp] = 0;
        
        // go to next line
        currentY[comp] += 8 * VScale;
        currentBlockVFactor[comp]++;
        
        // you made a column of blocks
        if (currentBlockVFactor[comp] >= VFactor) {
            
            // restore the current Y to its initial position and reset the counters
            currentY[comp] -= 8  * VScale * VFactor;
            currentBlockVFactor[comp] = 0;
            
            currentX[comp] += 8 * HScale * HFactor;
            if (currentX[comp] >= width) {
                currentX[comp] = 0;
                currentY[comp] += 8  * VScale * VFactor;
                
                // TODO: Progressive JPEG not supported
                // Force end here because Progressive JPEG file has more stuff after this
                if (currentY[comp] >= height && comp == components.size()-1) {
                    endOfFile = true;
                } // end if (currentX[comp] >= width)
            } // end if (currentBlockVFactor[comp] >= VFactor)
            
        } // end  if (currentBlockVFactor[comp] >= VFactor)
    } // end if(currentBlockHFactor[comp] >= HFactor)
    
} // end addBlockSubSampling



void jpeg_decoder::addtCoeffBlock(int dataBlock[8][8], int comp /* Current component */, int currentX, int currentY ) {
    
    //    cout << comp << ") Adding dct block at X: " << currentX_dct[comp]  << ", Y: " << currentY_dct[comp] << endl;
    
    // How many Y's in the horizontal and vertical direction (2x2 is the usual case)
    int HFactor = components[comp].HFactor, VFactor = components[comp].VFactor;
    
    // These Y's should be scaled/repeated how many times horizonatally and vertically
    int HScale = components[comp].HScale,   VScale = components[comp].VScale;
    
    // Note: you should use the upscaled width and upscaled height
    int comp_height = ceil(1.0 * upscale_height   /  VScale);
    int comp_width  = ceil(1.0 * upscale_width    / HScale);
    int width = comp_width, height = comp_height;
    
    
#if DEBUGLEVEL > 20
    if (currentX == 1008 && currentY == 0 && comp == COMPONENT_Cb) {
        cout << comp << ") Reporting information about the image and the component... " << endl;
        cout << "Width: " << width << " Height: " << height << endl;
        cout << "HScale: " << HScale << ", VScale: " << VScale << endl;
        cout << "Sent currentX , currentY " << currentX << ", " << currentY << endl;
        cout << "Calculated currentX, currentY  " << currentX_dct[comp]
        << ", " << currentY_dct[comp] << endl;
        
        
        cout << "Addtcoeff block init decoder input block: " << endl;
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                cout << dataBlock[i][j] << " ";
            }
            
            cout << "\n";
        }
        
        cout << "Bye " << endl;
    }
    
    
#endif
    int block_idx = (currentX_dct[comp] / 8) + (currentY_dct[comp] / 8) * (width / 8);
    
    // raster block sanning/traverse using x, y positions:
    for (int y = 0; y < 8; ++y) {
        if (currentX_dct[comp] >= width) break;
        if (currentY_dct[comp] + y >= height) break;
        
        int picture_y = currentY_dct[comp] + y;
        
        // NEW to TCM: store the DCT coefficient in terms of AC indexes
        
        for (int x = 0; x < 8; ++x) {
            
            if (currentX_dct[comp] + x >= width) break;
            
            int raw_dct_value = dataBlock[y][x];
            int realx = currentX_dct[comp] + x;
            int dct_coeff_idx = x + y * 8;
            // Set the pixel <Xcor: realX, Ycor: picture_y>
            switch (comp) {
                case COMPONENT_Y:
                    
                    tCoeff_Y[picture_y][realx] = raw_dct_value;
                    
                    // NEW to TCM: store the DCT coefficient in terms of AC indexes
                    tCoeff_Y_AC[dct_coeff_idx][block_idx] = raw_dct_value;
                    
                    //cout << "X: " << realx << ", Y: " << picture_y << ", val: " <<  tCoeff_Y[picture_y][realx] << endl;
                    break;
                case COMPONENT_Cb:
                    
                    
                    tCoeff_Cb[picture_y][realx] = raw_dct_value;
                    // cout << "X Cb: " << realx << ", Y Cb: " << picture_y << ", val: " <<  tCoeff_Cb[picture_y][realx] << endl;
                    break;
                case COMPONENT_Cr:
                    tCoeff_Cr[picture_y][realx] = raw_dct_value;
                    //                          cout << "X Cr: " << realx << ", Y Cr: " << picture_y << ", val: " <<  tCoeff_Cb[picture_y][realx] << endl;
                    break;
                default:
                    tCoeff_Y[picture_y][realx] = raw_dct_value;
                    cout << "-|- ##ERROR## You are sending the wrong component - will assume it's luminance. Most probably will fail!" << endl;
                    
                    break;
            } // end switch
            
            if (++realx >= width) break;
        } // end for x
    } // end for y
    
#if DEBUGLEVEL > 20
    
    if (currentX == 1008 && currentY == 0 && comp == COMPONENT_Cb) {
        
        int debug[8][8] = { 0 };
        
        
        cout << "Reverse operation in decoder: " << endl;
        cout << "A) Block: " << endl;
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                // debug[i][j] = dataBlock[i][j];
                
                int realx = currentX_dct[comp] + j;
                int realy = currentY_dct[comp] + i;
                debug[i][j] = tCoeff_Y[realy][realx];
                cout << debug[i][j] << " ";
                // cout << "X Cb: " << realx << ", Y Cb: " << realy << ", val: " << tCoeff_Cb[realy][realx] << endl;
            }
            
            cout << "\n";
        }
        
        
        cout << "B) Operation: " << endl;
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                int idx = j + i * 8;
                int q = components.at(comp).componentQuantizationTable->quantizationTableData[i][j];
                debug[i][j] = dataBlock[i][j] / q;
                
                cout << debug[i][j] << " ";
            }
            
            cout << "\n";
        }
        
        cout << "Yes I found it :P " << endl;
    }
    
#endif
    // currentX_dct and currentY_dct are now 8 pixels below values when function is called
    // they should be 8 pixels to the right multiplied by HScale
    // unless this is image edge, in which case they should be in the next line
    // Update starting X and Y for next block, taking into account subsampling
    
    if (!progressive_Huff_Format || (progressive_Huff_Format && zigZagStart == 0 && approximationL != 0)) {
        HScale = 1; // note: you will store the dct coefficients without any scaling
        VScale = 1; // note: algorithm below works, but you have to set the scales into 1x1
        currentX_dct[comp] += 8 * HScale;
        currentBlockHFactor_dct[comp]++;
        
        // you made a line of blocks
        if (currentBlockHFactor_dct[comp] >= HFactor) {
            
            // restore the current X to its initial position and reset the counters
            currentX_dct[comp] -= 8 * HScale * HFactor;
            currentBlockHFactor_dct[comp] = 0;
            
            // go to next line
            currentY_dct[comp] += 8 * VScale;
            currentBlockVFactor_dct[comp]++;
            
            // you made a column of blocks
            if (currentBlockVFactor_dct[comp] >= VFactor) {
                
                // restore the current Y to its initial position and reset the counters
                currentY_dct[comp] -= 8 * VScale * VFactor;
                currentBlockVFactor_dct[comp] = 0;
                
                currentX_dct[comp] += 8 * HScale * HFactor;
                if (currentX_dct[comp] >= width) {
                    currentX_dct[comp] = 0;
                    currentY_dct[comp] += 8 * VScale * VFactor;
                    
                    // TODO: Progressive JPEG not supported
                    // Force end here because Progressive JPEG file has more stuff after this
                    if (currentY_dct[comp] >= height && comp == components.size() - 1) {
                        endOfFile = true;
                    } // end if (currentX_dct[comp] >= width)
                } // end if (currentBlockVFactor_dct[comp] >= VFactor)
                
            } // end  if (currentBlockVFactor_dct[comp] >= VFactor)
        } // end if(currentBlockHFactor_dct[comp] >= HFactor)}
    }
    else if ((progressive_Huff_Format && zigZagStart != 0) || (progressive_Huff_Format && zigZagStart == 0 && approximationL == 0)){
        currentX_dct[comp] += 8;
        if (currentX_dct[comp] + 8 > width) {
            currentX_dct[comp] = 0;
            currentY_dct[comp] += 8;
        }
    }
} // end addtCoeffBlock

void jpeg_decoder::final_process_progressive() {
    for (int currentComponent = COMPONENT_Y; currentComponent < components.size(); ++currentComponent) {
        int HFactor = components[currentComponent].HFactor, VFactor = components[currentComponent].VFactor;
        // These Y's should be scaled/repeated how many times horizonatally and vertically
        int HScale = components[currentComponent].HScale, VScale = components[currentComponent].VScale;
        
        // Note: use the upscaled width and height to store the DCT coefficientss
        int comp_height = ceil(1.0* upscale_height / VScale);
        int comp_width = ceil(1.0* upscale_width / HScale);
        int width = comp_width, height = comp_height;
        
        // Initialize the positions
        int currentX = 0;
        int currentY = 0;
        int currentBlockHFactor = 0;
        int currentBlockVFactor = 0;
        
        // loop on the entire picture
        for (uint i = 0; i < comp_height; i += 8) {
            for (uint j = 0; j < comp_width; j += 8) {
                int dataReshapedInto8x8[8][8] = { 0 };
                // fetch the block
                uint y;
                for (y = 0; y < 8; ++y) {
                    //vector<vector<double> >dataReshapedInto8x8(8, vector<double>(8));
                    
                    if (currentX >= width) break;
                    if (currentY + y >= height) break;
                    
                    int picture_y = currentY + y;
                    
                    uint x;
                    for (x = 0; x < 8; ++x) {
                        
                        if (currentX + x >= width) break;
                        int realx = currentX + x;
                        
                        if (currentComponent == COMPONENT_Y)
                            dataReshapedInto8x8[y][x] = tCoeff_Y[picture_y][realx];
                        else if (currentComponent == COMPONENT_Cb)
                            dataReshapedInto8x8[y][x] = tCoeff_Cb[picture_y][realx];
                        else
                            dataReshapedInto8x8[y][x] = tCoeff_Cr[picture_y][realx];
                    }
                }
                
                // Perform dequantization
                multiplyWithQuantizationTable(dataReshapedInto8x8, currentComponent);
                
                // Perform inverse zig zag scanning
                int dataReshapedInto8x8Zig[8][8] = { 0 };
                inverseZigZagScanning(dataReshapedInto8x8Zig, dataReshapedInto8x8);
                
                //store performed data (inverse ZigZag / Dequantization) back to tCoeff
                for (y = 0; y < 8; ++y) {
                    int picture_y = currentY + y;
                    uint x;
                    for (x = 0; x < 8; ++x) {
                        int realx = currentX + x;
                        if (currentComponent == COMPONENT_Y)
                            tCoeff_Y[picture_y][realx] = dataReshapedInto8x8Zig[y][x];
                        else if (currentComponent == COMPONENT_Cb)
                            tCoeff_Cb[picture_y][realx] = dataReshapedInto8x8Zig[y][x];
                        else
                            tCoeff_Cr[picture_y][realx] = dataReshapedInto8x8Zig[y][x];
                    }
                }
                
                // Perform IDCT
                int output_value_dct[8][8] = { 0 };
                perform_idct(output_value_dct, dataReshapedInto8x8Zig);
                addBlockSubsampling(output_value_dct, currentComponent);
 

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
    }
}

void jpeg_decoder::IDCT (int block[8][8], int transformedBlock[8][8]) {
    int sum = 0;
    int counter = 0;
    
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            sum = 0;
            for(int u = 0; u < 8; ++u){
                for(int v = 0 ;v < 8; ++v){
                    sum = sum + block[u][v] * coefficients[ counter++ ];
                } // end v loop
            } // end u loop
            
            // All coefficients are multiplied by 1024 since they are int
            sum = sum >> 10;
            
            // Only 8 and 12-bit is supported by DCT per JPEG standard
            if (precision == 8) {
                transformedBlock[i][j] = sum + 128;
            }
            else if (precision == 12) {
                transformedBlock[i][j] = sum + 2048;
            }
        } // end j loop
    } // end i loop
} // end IDCT


void jpeg_decoder::write_bmp24(const std::string output_bmp_filename, int Width, int Height, unsigned char *RGB) {
    
#pragma pack(1)
    struct stBMFH // BitmapFileHeader & BitmapInfoHeader
    {
        // BitmapFileHeader
        char         bmtype[2];     // 2 bytes - 'B' 'M'
        unsigned int iFileSize;     // 4 bytes
        short int    reserved1;     // 2 bytes
        short int    reserved2;     // 2 bytes
        unsigned int iOffsetBits;   // 4 bytes
        // End of stBMFH structure - size of 14 bytes
        // BitmapInfoHeader
        unsigned int iSizeHeader;    // 4 bytes - 40
        unsigned int iWidth;         // 4 bytes
        unsigned int iHeight;        // 4 bytes
        short int    iPlanes;        // 2 bytes
        short int    iBitCount;      // 2 bytes
        unsigned int Compression;    // 4 bytes
        unsigned int iSizeImage;     // 4 bytes
        unsigned int iXPelsPerMeter; // 4 bytes
        unsigned int iYPelsPerMeter; // 4 bytes
        unsigned int iClrUsed;       // 4 bytes
        unsigned int iClrImportant;  // 4 bytes
        // End of stBMIF structure - size 40 bytes
        // Total size - 54 bytes
    };
#pragma pack()
    
    // Round up the width to the nearest DWORD boundary
    int iNumPaddedBytes = 4 - (Width * 3) % 4;
    iNumPaddedBytes = iNumPaddedBytes % 4;
    
    stBMFH bh;
    memset(&bh, 0, sizeof(bh));
    bh.bmtype[0]='B';
    bh.bmtype[1]='M';
    bh.iFileSize = (Width*Height*3) + (Height*iNumPaddedBytes) + sizeof(bh);
    bh.iOffsetBits = sizeof(stBMFH);
    bh.iSizeHeader = 40;
    bh.iPlanes = 1;
    bh.iWidth = Width;
    bh.iHeight = Height;
    bh.iBitCount = 24;
    FILE* fp = fopen(output_bmp_filename.c_str(), "wb");
    
    
    fwrite(&bh, sizeof(bh), 1, fp);
    for (int y = Height-1; y>=0; y--)
    {
        for (int x = 0; x < Width; x++)
        {
            int i = (x + (Width)*y) * 3;
            unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
            fwrite(&rgbpix, 3, 1, fp);
        }
        if (iNumPaddedBytes>0)
        {
            unsigned char pad = 0;
            fwrite(&pad, iNumPaddedBytes, 1, fp);
        }
        
        
    }
    fclose(fp);
    
    
#if DEBUGLEVEL > 50
    cout << "RGB Picture exit is here: " << endl;
    for (int y=Height-1; y>=0; y--)
    {
        for (int x=0; x<Width; x++)
        {
            int i = (x + (Width)*y) * 3;
            unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
            printf("%d ", rgbpix);
        }
        
        if (iNumPaddedBytes>0)
        {
            unsigned char pad = 0;
            printf("%d ", pad);
        }
        
        printf("\n");
    }
#endif
    
    std::size_t found = output_bmp_filename.find_last_of("/\\");
    std::string name_file_only = output_bmp_filename.substr(found+1);
    cout << "Success: " << name_file_only << " has been created!" << endl;
} // end write_bmp24


#if OPEN_CV_ENABLED
int jpeg_decoder::display_jpg_yuv(std::string window_name, uint_32 comp) const {
    
    if (m_Y==NULL)
    {
        printf("Y channel buffer is null");
        return 0;
    }
    
    if(comp >= components.size())
    {
        cout << "Number of components in the bitstream is " << components.size() << ", but input comp is: " << comp << ". Will set comp to Y'" << endl;
        comp = COMPONENT_Y;
    }
    
    uint_8** ptr;
    std::string compType = " Y'";
    switch(comp) {
        case COMPONENT_Y:
            ptr = m_YPicture_buffer;
            compType = " Y'";
            break;
        case COMPONENT_Cb:
            ptr = m_CbPicture_buffer;
            compType = " Cb";
            break;
        case COMPONENT_Cr:
            ptr = m_CrPicture_buffer;
            compType = " Cr";
            break;
        default:
            ptr = m_YPicture_buffer;
            compType = " Y'";
            break;
    } // end switch
    
    window_name += " " + compType;
    
    uint8_t* temp = new uint8_t[jpegImageWidth*jpegImageHeight];
    
    for(int i = 0; i < jpegImageHeight; ++i) {
        for(int j = 0; j < jpegImageWidth; ++j) {
            int idx = j + i * jpegImageWidth;
            temp[idx] = ptr[i][j];
        }
    }
    
    cv::Mat img(jpegImageHeight, jpegImageWidth, CV_8U, temp, cv::Mat::AUTO_STEP);
    imshow(window_name, img);
    cv::waitKey(0);
    cout << "Success: " << compType << " has been displayed!" << endl;
    delete [] temp;
    return 1;
    
} // end display_jpg_yuv


int jpeg_decoder::display_jpg_bmp(const std::string output_bmp_filename) const {
    
    std::stringstream oss;
    std::size_t found = output_bmp_filename.find_last_of(".");
    std::string path_with_name = output_bmp_filename.substr(0, found);
    found = output_bmp_filename.find_last_of("/\\");
    std::string name_file_only = output_bmp_filename.substr(found+1);
    
    
    oss << path_with_name << "_" << jpegImageWidth << "x" << jpegImageHeight << ".yuv";
    cv::Mat img = cv::imread(output_bmp_filename, CV_LOAD_IMAGE_COLOR);
    imshow(name_file_only, img);
    cv::waitKey(0);
    
    
    oss.clear(); oss.str("");
    oss << "BMP_File" << "_" << jpegImageWidth << "x" << jpegImageHeight << ".bmp";
    cout << "Success: " << oss.str() << " has been displayed!" << endl;
    return 1;
    
} // end display_jpg_bmp

#endif // end open_cv_enabled

// Note: file path has to be path_to_files + filename.yuv (filename should not contain a .)
int jpeg_decoder::write_yuv_from_jpg_in_yuv(const std::string output_yuv_filename) {
    
    std::stringstream oss;
    std::size_t found = output_yuv_filename.find_last_of(".");
    std::string path_with_name = output_yuv_filename.substr(0, found);
    found = output_yuv_filename.find_last_of("/\\");
    std::string name_file_only = path_with_name.substr(found+1);
    
    oss << path_with_name << "_" << jpegImageWidth << "x" << jpegImageHeight << ".yuv";
    std::string output = oss.str();
    FILE* pFile = fopen(output.c_str(), "wb");
    uint_8** ptr = m_YPicture_buffer;
    
    for(int i = 0; i < jpegImageHeight; ++i) {
        for(int j = 0; j < jpegImageWidth; ++j) {
            uint_8 uc = ptr[i][j];
            fwrite( &uc, sizeof(uint_8), 1, pFile );
        }
    }
    fclose(pFile);
    
    
    oss.clear(); oss.str("");
    oss << name_file_only << "_" << jpegImageWidth << "x" << jpegImageHeight << ".yuv";
    cout << "Success: " << oss.str() << " has been created!" << endl;
    return 1;
    
} // end write_yuv_from_jpg_in_yuv


int jpeg_decoder::write_yuv_from_jpg_in_csv(const std::string output_csv_name) {
    
    ofstream myfile;
    myfile.open (output_csv_name);
    uint_8** ptr = m_YPicture_buffer;
    std::stringstream oss;
    std::size_t found = output_csv_name.find_last_of(".");
    std::string path_with_name = output_csv_name.substr(0, found);
    found = output_csv_name.find_last_of("/\\");
    std::string name_file_only = path_with_name.substr(found+1);
    
    for(int i = 0; i < jpegImageHeight; ++i) {
        for(int j = 0; j < jpegImageWidth; ++j) {
            uint_8 uc = ptr[i][j];
            myfile << static_cast<int>(uc) << ",";
        }
        
        if( (i + 1) < jpegImageHeight){
            myfile << "\n";
        }
    }
    
    myfile.close();
    
    oss.clear(); oss.str("");
    oss << name_file_only << "_" << jpegImageWidth << "x" << jpegImageHeight << ".csv";
    cout << "Success: " << oss.str() << " has been created!" << endl;
    return 1;
    
} // end write_yuv_from_jpg_in_yuv


int jpeg_decoder::convert_jpg_to_bmp(const std::string output_bmp_filename) {
    
    if (m_rgb==NULL)
    {
        printf("Failed to decode jpg and RGB buffer is null\n");
        return 0;
    }
    
    // save the output bmp filename
    jpeg_output_bmp_filename = output_bmp_filename;
    
    // you have already read the JPEG image - write it directly to BMP
    write_bmp24(output_bmp_filename, jpegImageWidth, jpegImageHeight, m_rgb);
    
    
#if OPEN_CV_ENABLED && DISPLAY_BMP_IMAGE
    display_jpg_bmp(output_bmp_filename);
#endif
    
    // TODO: you should think about how you will release the rgb memory here
    return 1;
    
} // end convert_jpg_to_bmp


int jpeg_decoder::write_tcoeff_y_from_jpg_in_csv(const std::string output_csv_name, int comp) {
    ofstream myfile;
    myfile.open (output_csv_name, std::ofstream::out);
    
    std::stringstream oss;
    std::size_t found = output_csv_name.find_last_of(".");
    std::string path_with_name = output_csv_name.substr(0, found);
    found = output_csv_name.find_last_of("/\\");
    std::string name_file_only = path_with_name.substr(found+1);
    
    // These Y's should be scaled/repeated how many times horizonatally and vertically
    int HScale = components[comp].HScale,   VScale = components[comp].VScale;
    int comp_height = ceil(jpegImageHeight / VScale);
    int comp_width  = ceil(jpegImageWidth  /  HScale);
    
    
    cout << "comp: " << comp << endl;
    cout << "comp_width " << comp_width << ", " << comp_height << endl;
    for(int i = 0; i < comp_height; ++i) {
        for(int j = 0; j < comp_width; ++j) {
            int uc = tCoeff_Y[i][j];
            
            if(comp == COMPONENT_Cb) {
                uc = tCoeff_Cb[i][j];
            }
            else if(comp == COMPONENT_Cr) {
                uc = tCoeff_Cr[i][j];
            }
            
            myfile << static_cast<int>(uc) << ",";
        }
        
        if( (i + 1) < comp_height){
            myfile << "\n";
        }
    }
    myfile.close();
    
    oss.clear(); oss.str("");
    oss << name_file_only << "_" << comp_width << "x" << comp_height << ".csv";
    cout << "Success: " << oss.str() << " has been created!" << endl;
    return 1;
} // end write_tcoeff_y_from_jpg_in_csv


void jpeg_decoder::initPositionsBuffersForPictureBuffer() {   
    // [Buffers for DCT coefficients read from bitstream: ]
    // Create the luminance 2D buffer for the entire picture:
    // ***************** Change to upscaled size to fit the multiple of 8
    int comp_height = ceil(1.0* upscale_height / components.at(COMPONENT_Y).VScale);
    int comp_width = ceil(1.0* upscale_width / components.at(COMPONENT_Y).HScale);
    
    // create the transformed coeffiecients:
    tCoeff_Y.resize(comp_height, vector<int> (comp_width, 0));
    
    // NEW for TCM: create the tCoeff_Y_AC to store the DCT coefficient in terms of AC indexes
    tCoeff_Y_AC.resize(64, vector<int>(total_block_Y, 0));
    
    
    if(components.size() >= COMPONENT_Cb + 1) {
        // jpegImageHeight * jpegImageWidth
        comp_height = ceil(1.0* upscale_height / components.at(COMPONENT_Cb).VScale);
        comp_width  = ceil(1.0* upscale_width / components.at(COMPONENT_Cb).HScale);
        
        // create the transformed coeffiecients:
        tCoeff_Cb.resize(comp_height, vector<int>(comp_width, 0));
        
    }
    
    if(components.size() >= COMPONENT_Cr + 1) {
        comp_height = ceil(1.0* upscale_height / (double) components.at(COMPONENT_Cr).VScale);
        comp_width = ceil(1.0* upscale_width / (double) components.at(COMPONENT_Cr).HScale);
        
        // create the transformed coeffiecients:
        tCoeff_Cr.resize(comp_height, vector<int>(comp_width, 0));
    }
    
    // Create the luminance 2D buffer for the entire picture:
    m_YPicture_buffer = new uint_8*[upscale_height]; // height = rows //NEW: Changed to upscaled size for smoothing purpose
    for(int i = 0; i <  upscale_height; ++i){
        m_YPicture_buffer[i] = new uint_8[upscale_width]; // width = columns  //NEW: Changed to upscaled size for smoothing purpose
    }
    
    // Create the Cb 2D buffer for the entire picture:
    m_CbPicture_buffer = new uint_8*[upscale_height];
    for(int i = 0; i < upscale_height; ++i){
        m_CbPicture_buffer[i] = new uint_8[upscale_width];
    }
    
    // Create the Cr 2D buffer for the entire picture:
    m_CrPicture_buffer = new uint_8*[upscale_height];
    for(int i = 0; i < upscale_height; ++i){
        m_CrPicture_buffer[i] = new uint_8[upscale_width];
    }
    
    initCurrentPosition();
    
} // end initPositionsBuffersForPictureBuffer

void jpeg_decoder::initCurrentPosition() {
    // Initialize the positions (where you are currently at in the picture):
    for (uint i = 0; i < components.size(); ++i) {
        // Set the co-ordinates of the next block (x, y)
        currentX[i] = currentY[i] = 0;
        
        // intitialize BlockHFactor and BlockVFactor (used for subsampling)
        // counts how many times you scaled your block
        currentBlockHFactor[i] = currentBlockVFactor[i] = 0;
        
        // Set the co-ordinates of the next DCT block (x, y)
        currentX_dct[i] = currentY_dct[i] = 0;
        
        // intitialize BlockHFactor and BlockVFactor (used for dct)
        // counts how many times you scaled your block
        currentBlockHFactor_dct[i] = currentBlockVFactor_dct[i] = 0;
    }
}

void jpeg_decoder::deleteDCTComponentPointers() {
    
    // Delete the buffers for holding the Y'CbCr data
    for(int iComponent = 0; iComponent < components.size(); ++iComponent){
        if(components.at(iComponent).m_DCT != nullptr)
        {
            delete [] components.at(iComponent).m_DCT;
            components.at(iComponent).m_DCT = nullptr;
        }
    }
} // end deleteDCTComponentPointers

void jpeg_decoder::deleteRawPictureBufferPointers() {
    
    // Delete the buffers for holding the Y'CbCr data
    for(int i = 0; i < upscale_height; ++i){
        // delete each sub-array:
        if(m_YPicture_buffer[i] != nullptr) {
            delete [] m_YPicture_buffer[i];
            m_YPicture_buffer[i] = nullptr;
        }
        
        if(m_CbPicture_buffer[i] != nullptr) {
            delete [] m_CbPicture_buffer[i];
            m_CbPicture_buffer[i] = nullptr;
        }
        
        if(m_CrPicture_buffer[i] != nullptr) {
            delete [] m_CrPicture_buffer[i];
            m_CrPicture_buffer[i] = nullptr;
        }
        
    } // end loop
    
    if(m_YPicture_buffer != nullptr) {
        delete [] m_YPicture_buffer;
        m_YPicture_buffer = nullptr;
    }
    
    if(m_CbPicture_buffer != nullptr) {
        delete [] m_CbPicture_buffer;
        m_CbPicture_buffer = nullptr;
    }
    
    if(m_CrPicture_buffer != nullptr) {
        delete [] m_CrPicture_buffer;
        m_CrPicture_buffer = nullptr;
    }
    
} // end deleteRawPictureBufferPointers

// decoding functions from the other source
/***************************************************************************/

void DeZigZag(int outBlock[64], const int inBlock[64])
{
    constexpr static int ZigZagArray[64] =
    {
        0,   1,   5,  6,   14,  15,  27,  28,
        2,   4,   7,  13,  16,  26,  29,  42,
        3,   8,  12,  17,  25,  30,  41,  43,
        9,   11, 18,  24,  31,  40,  44,  53,
        10,  19, 23,  32,  39,  45,  52,  54,
        20,  22, 33,  38,  46,  51,  55,  60,
        21,  34, 37,  47,  50,  56,  59,  61,
        35,  36, 48,  49,  57,  58,  62,  63,
    };
    for(int i=0; i<64; i++){
        outBlock[ i ] = inBlock[ZigZagArray[i]];
    }
}


/***************************************************************************/

// debugging functions:

void jpeg_decoder::dumpBlock(int block[8][8]) {
    cout << "Dump block:";
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
            cout << setw(5) << block[i][j] << " ";
        cout << endl;
    }
} // end dumpBlock

void jpeg_decoder::printDCHuffmanCodes(int currentComponent) {
    HuffmanTable * hTable = componentTablesDC[currentComponent];
    
    printf("Huffman table ID #%02X class 00 (DC Table):\n", hTable->tableID);
    // Once the map has been built, print it out
    std::map<huffKey, uint_8>::iterator iter;
    for (iter = hTable->huffData.begin();iter != hTable->huffData.end(); ++iter) {
        
        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X or [%s] \n",
               iter->first.second, iter->first.first,
               iter->second, IntToBinary(iter->first.second, iter->first.first));
    } // end for loop
} // end printDCHuffmanCodes


void jpeg_decoder::printACHuffmanCodes(int currentComponent) {
    HuffmanTable * hTable = componentTablesAC[currentComponent];
    
    printf("Huffman table ID #%02X class 01 (AC Table):\n", hTable->tableID);
    // Once the map has been built, print it out
    std::map<huffKey, uint_8>::iterator iter;
    for (iter = hTable->huffData.begin();iter != hTable->huffData.end(); ++iter) {
        // Print Code - Its Length : Equivalent letter
        printf("    %04X at length %d = %02X or [%s] \n",
               iter->first.second, iter->first.first,
               iter->second, IntToBinary(iter->first.second, iter->first.first));
    } // end for loop
} // end printACHuffmanCodes


void jpeg_decoder::dumpDCTValues(short dct[64]) {
    printf("\n#Extracted DCT values from SOS#\n");
    int c = 0;
    for (int i=0; i<64; i++)
    {
        printf("% 4d  ", dct[c++]);
        
        if ( (c>0) && (c%8==0) ) printf("\n");
    }
    printf("\n");
} // end dumpDCTValues

void jpeg_decoder::dumpDecodedBlock(int val[8][8]) {
    printf("# Decoded 8x8 Block#\n");
    for( int y=0; y<8; y++){
        for( int x=0; x<8; x++)
        {
            printf("%2x ", val[x][y]);
        }
        printf("\n");
    }
} // end dumpDecodedBlock


void jpeg_decoder::dumpDecodedBlockInDecimal(int val[8][8]) {
    printf("# Decoded 8x8 Block#\n");
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            cout << std::dec <<  val[i][j] << " ";
        }
        cout << "\n";
    }
} // end dumpDecodedBlock

void jpeg_decoder::dumpQuantizationTable(int currentComponent) {
    
    if(currentComponent){
        printf("Quantization table ID #%02X class 01 (chrominance):\n", currentComponent);
    }else {
        printf("Quantization table ID #%02X class 00 (luminance):\n", currentComponent);
    }
    
    QuantizationTable* qTable = this->components[currentComponent].componentQuantizationTable;
    
    cout << "Quantization Table: " << endl;
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 8; ++j){
            if(qTable->quantizationTableData[i][j] < 10){
                cout << " " << qTable->quantizationTableData[i][j] << " ";
            }else{
                cout  << qTable->quantizationTableData[i][j] << " ";
            }
            
        }
        cout << "\n";
    }   
} // end dumpQuantizationTable
