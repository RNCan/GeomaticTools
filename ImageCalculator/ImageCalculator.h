//***********************************************************************
#pragma once

#include <boost/dynamic_bitset.hpp>
#include <deque>
#include "geomatic/UtilGDAL.h"
#include "geomatic/LandsatDataset2.h"


namespace WBSF
{

class CImageCalculatorOption : public CBaseOptions
{
public:

    enum TFilePath		{ INPUT_FILE_PATH, OUTPUT_FILE_PATH, NB_FILE_PATH };


    CImageCalculatorOption();
    virtual ERMsg ParseOption(int argc, char* argv[]);
    virtual ERMsg ProcessOption(int& i, int argc, char* argv[]);


    std::vector<std::pair<std::string, std::string > > m_equations;
};

typedef std::deque < std::vector<double>> OutputData;

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
    void ReadBlock(CGDALDatasetExVector& inputDS, int xBlock, int yBlock, CRasterWindow& bandHolder);
    void ProcessBlock(int xBlock, int yBlock, const CRasterWindow& bandHolder, OutputData& outputData);
    void WriteBlock(int xBlock, int yBlock, CGDALDatasetEx& outputDS, OutputData& outputData);
    void CloseAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS);

    CImageCalculatorOption m_options;

    static const char* VERSION;
    static const size_t NB_THREAD_PROCESS;
};
}
