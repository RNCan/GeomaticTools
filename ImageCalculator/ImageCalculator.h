//***********************************************************************
#pragma once

#include <deque>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include <SDKDDKVer.h>


#include "Geomatic/UtilGDAL.h"
#include "Geomatic/ImageParser.h"


namespace WBSF
{


class CImageCalculatorOption : public CBaseOptions
{
public:

    //enum TFilePath { INPUT_FILE_PATH, OUTPUT_FILE_PATH, NB_FILE_PATH };


    CImageCalculatorOption();
    virtual ERMsg ParseOption(int argc, char* argv[]);
    virtual ERMsg ProcessOption(int& i, int argc, char* argv[]);


    std::vector<std::string> m_equations;
};



typedef std::deque < std::vector<float>> OutputData;

//typedef std::pair<double, size_t> NBRPair;
typedef std::vector<CGDALDatasetEx> CGDALDatasetExVector;

class CImageCalculator
{
public:

    ERMsg Execute();

    std::string GetDescription()
    {
        return  std::string("ImageCalculatorImages version ") + VERSION + " (" + __DATE__ + ")";
    }

    ERMsg OpenAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS);
    void ReadBlock(CGDALDatasetExVector& inputDS, int xBlock, int yBlock, std::deque<CRasterWindow>& windows);
    void ProcessBlock(int xBlock, int yBlock, const std::deque<CRasterWindow>& windows, OutputData& outputData);
    void WriteBlock(int xBlock, int yBlock, CGDALDatasetEx& outputDS, OutputData& outputData);
    ERMsg CloseAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS);

    //CVariableToData GetVariableToData(std::vector<std::pair<std::string, std::string > > equations);



    CImageCalculatorOption m_options;

    std::deque < std::deque< CImageParser > > m_parser;
    std::deque < boost::dynamic_bitset<> > m_bands_to_read;


    static const char* VERSION;
    static const size_t NB_THREAD_PROCESS;
};
}
