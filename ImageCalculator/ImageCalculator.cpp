//***********************************************************************
// program to merge Landsat image image over a period
//
//***********************************************************************
// version
// 4.0.0	02/09/2025	Rémi Saint-Amant	cross-platform version. Replace MTParser by muparser.
// 3.2.0	13/09/2023	Rémi Saint-Amant	Compile with GDAL 3.7.1
// 3.1.0	20/12/2021	Rémi Saint-Amant	Compile with VS 2019 and GDAL 3.0.3
// 3.0.1	10/10/2018	Rémi Saint-Amant	Compile with VS 2017. Add -BLOCK_THREADS
// 3.0.0	03/11/2017	Rémi Saint-Amant	Compile with GDAL 2.02
// 2.0.1	13/06/2015	Rémi Saint-Amant	Add -hist option
// 2.0.0    11/03/2015	Rémi Saint-Amant	New Equation form "E=mc²"
// 1.7.5	05/02/2015	Rémi Saint-Amant	Bug correction in UnionExtent
// 1.7.4	30/01/2015	Rémi Saint-Amant	don't modify input VRT file. Bug correction in IntersectRect
// 1.7.3    27/01/2015	Rémi Saint-Amant	Compile with GDAL 1.11.1
// 1.7.2    25/10/2014	Rémi Saint-Amant	bug correction in help
// 1.7.1    24/07/2014	Rémi Saint-Amant	Bug corrections
// 1.7.0	26/06/2014	Rémi Saint-Amant	GDAL 1.11, UNICODE, VC 2013
// 1.6.4	22/07/2013	Rémi Saint-Amant	Add progress for -Overview and -Stats
// 1.6.3	22/07/2013	Rémi Saint-Amant	Recompilation of cleaned code. Add -stats options
//											Better management of noData NAN and infinite number
//											Limit of no data to -1.0e38 and -1.0e308
// 1.6.2	18/04/2013	Rémi Saint-Amant	new compilation, bug correction in mask
// 1.6.1	12/04/2013  Rémi Saint-Amant	new compilation without debug info. bug correction in write time.
// 1.6.0	27/01/2013  Rémi Saint-Amant	Read by block instead of by row
// 1.5.0	06/11/2012	Rémi Saint-Amant	Improvements
//											ToDo : include new frameworks for image
// 1.4.0	05/04/2012  Rémi Saint-Amant	correction of bug in ManageNoData
// 1.3.0	03/04/2012  Rémi Saint-Amant	correction of bug in omp
// 1.2.0	15/03/2012	Rémi Saint-Amant	Little changes
// 1.1.0	14/02/2012	Rémi Saint-Amant	Bug correction in the srcNoData
// 1.0.0	30/01/2012	Rémi Saint-Amant	Initial version


//D:\Travaux\Landsat\Landsat(2000-2018)
//Input\Median_2000-2018.vrt Output\Median_2000_NBR.vrt -e "NBR=(i1b4-i1b7)/(i1b4+i1b7)" -of VRT -overwrite -co "COMPRESS=LZW" -te 1022538.9 6663106.0 1040929.5 6676670.7 -multi



#include <cmath>
#include <array>
#include <utility>
#include <iostream>


#include <boost/filesystem.hpp>


#include "basic/OpenMP.h"
#include "Basic/muparser/muParser.h"
#include "geomatic/UtilGDAL.h"


#include "ImageCalculator.h"




using namespace std;


namespace WBSF
{
const char* CImageCalculator::VERSION = "4.0.0";
const size_t CImageCalculator::NB_THREAD_PROCESS = 2;


//*********************************************************************************************************************

CImageCalculatorOption::CImageCalculatorOption()
{

    m_appDescription = "This software compute bands from many inputs images into destination image using user define math equations. This software use muParser to evaluate equations.";

    static const COptionDef OPTIONS[] =
    {
        {"-Equation",1,"\"formulas\"",true,"Equation to evaluate. One equation per output band. Variables in equation must have the form \"Name=i#b#\" where Name is the name of the output variable (use when output is VRT), # is 1 to the number of images (i) or bands (b)."},
        { "-e", 1, "\"formulas\"", true, "Equilvalent to -Equation." },
        {"srcfile",0,"",true, "Input images file path (no limits)."},
        {"dstfile",0,"",false, "Output image file path."}
    };


    for (size_t i = 0; i < sizeof(OPTIONS) / sizeof(COptionDef); i++)
        AddOption(OPTIONS[i]);



    static const CIOFileInfoDef IO_FILE_INFO[] =
    {
        {"Input image","srcfile","1","*","","Multi-layer input"},
        {"Output Image", "dstfile","1","*","","Multi-layer output. One layer by equation"},
    };

    for (size_t i = 0; i < sizeof(IO_FILE_INFO) / sizeof(CIOFileInfoDef); i++)
        AddIOFileInfo(IO_FILE_INFO[i]);
}

ERMsg CImageCalculatorOption::ParseOption(int argc, char* argv[])
{
    ERMsg msg = CBaseOptions::ParseOption(argc, argv);

    if (msg)
    {
        if (m_filesPath.size() < 2)
        {
            msg.ajoute("ERROR: Invalid argument line. At least 2 images is needed: input and output image.\n");
            msg.ajoute("Argument found: ");
            for (size_t i = 0; i < m_filesPath.size(); i++)
                msg.ajoute("   " + to_string(i + 1) + "- " + m_filesPath[i]);
        }

        if (m_equations.size() == 0)
            msg.ajoute("ERROR: Invalid argument line. At least 1 equation must be define.\n");
    }

    if (m_outputType == GDT_Unknown)
        m_outputType = GDT_Float32;

    if (m_dstNodata == MISSING_NO_DATA)
        m_dstNodata = WBSF::GetDefaultNoData(GDT_Int16);//use Int16 missing value



    return msg;
}





ERMsg CImageCalculatorOption::ProcessOption(int& i, int argc, char* argv[])
{
    ERMsg msg;
    if (IsEqual(argv[i], "-Equation") || IsEqual(argv[i], "-e"))
    {
        m_equations.push_back(argv[++i]);
    }

    /*else if (IsEqual(argv[i], "-ManageSrcNoData"))
    {
    	m_bManageSrcNoData = true;
    }*/
    else
    {
        //it's a base option
        msg = CBaseOptions::ProcessOption(i, argc, argv);
    }

    return msg;
}



ERMsg CImageCalculator::Execute()
{
    ERMsg msg;


    size_t nbImages = m_options.m_filesPath.size();

    if (!m_options.m_bQuiet)
    {
        cout << "Output: " << m_options.m_filesPath[nbImages - 1] << endl;
        cout << "From:   " << endl;
        for (size_t i = 0; i < nbImages - 1; i++)
            cout << "   " << m_options.m_filesPath[i] << endl;
        cout << "Equations:" << endl;
        for (size_t i = 0; i < m_options.m_equations.size(); i++)
            cout << "   " << m_options.m_equations[i] << endl;

        cout << endl;
    }



    GDALAllRegister();

    CGDALDatasetExVector inputDSvector;
    //CGDALDatasetEx maskDS;
    CGDALDatasetEx outputDS;


    msg = OpenAll(inputDSvector, outputDS);
    if (!msg)
        return msg;


    assert(m_options.m_CPU > 0);
    assert(!m_options.m_equations.empty());
    m_parser.resize(m_options.m_CPU);
    for (size_t t = 0; t < m_parser.size() && msg; t++)
    {
        m_parser[t].resize(m_options.m_equations.size());

        for (size_t e = 0; e < m_parser[t].size(); e++)
            msg += m_parser[t][e].Compile(m_options.m_equations[e], inputDSvector);
    }


    if (msg)
    {
        assert(!m_parser.empty());

        m_bands_to_read.resize(inputDSvector.size());
        for (size_t e = 0; e < m_parser[0].size(); e++)
        {
            const std::deque<CVarDef>& vars = m_parser[0][e].GetParsedVars();
            for (size_t v = 0; v < vars.size(); v++)
            {
                size_t i = vars[v].m_ib.m_image;
                size_t b = vars[v].m_ib.m_band;
                if (m_bands_to_read[i].empty())
                    m_bands_to_read[i].resize(inputDSvector[i].GetRasterCount());

                m_bands_to_read[i].set(b);
            }

        }

        if (!m_options.m_bQuiet && m_options.m_bCreateImage)
            cout << "Create output images (" << outputDS.GetRasterXSize() << " C x " << outputDS.GetRasterYSize() << " R x " << outputDS.GetRasterCount() << " B) with " << m_options.m_CPU << " threads..." << endl;

        CGeoExtents extents = m_options.m_extents;
        m_options.ResetBar((size_t)extents.m_xSize * extents.m_ySize);

        vector<pair<int, int>> XYindex = extents.GetBlockList();
        map<int, bool> treadNo;


        //omp_set_nested(1);//for inner parallel llop
        //#pragma omp parallel for schedule(static, 1) num_threads( NB_THREAD_PROCESS ) if (m_options.m_bMulti)
        for (int b = 0; b < (int)XYindex.size(); b++)
        {
            int xBlock = XYindex[b].first;
            int yBlock = XYindex[b].second;

            deque<CRasterWindow> inputData;
            OutputData outputData;

            ReadBlock(inputDSvector, xBlock, yBlock, inputData);
            ProcessBlock(xBlock, yBlock, inputData, outputData);
            WriteBlock(xBlock, yBlock, outputDS, outputData);

        }//for all blocks
    }

    msg += CloseAll(inputDSvector, outputDS);

    return msg;
}



ERMsg CImageCalculator::OpenAll(CGDALDatasetExVector& inputDSvector, CGDALDatasetEx& outputDS)
{
    ERMsg msg;

    if (!m_options.m_bQuiet)
        cout << endl << "Open input image..." << endl;

    //msg = inputDS.OpenInputImage(m_options.m_filesPath[CImageCalculatorOption::INPUT_FILE_PATH], m_options);
    size_t nbImages = m_options.m_filesPath.size();
    //Get input file path in
    vector < string > filePathIn(nbImages - 1);
    for (size_t i = 0; i < nbImages - 1; i++)
        filePathIn[i] = m_options.m_filesPath[i];

    string filePathOut = m_options.m_filesPath[nbImages - 1];


    array<CGeoExtents, 2> extents;
    inputDSvector.resize(filePathIn.size());
    for (size_t i = 0; i < filePathIn.size() && msg; i++)
    {
        msg += inputDSvector[i].OpenInputImage(filePathIn[i], m_options);

        if (msg)
        {
            //compute minimum and maximum extents
            CImageCalculatorOption op;
            inputDSvector[i].UpdateOption(op);
            extents[0].IntersectExtents(op.m_extents, CGeoExtents::RES_MAX);
            extents[1].UnionExtents(op.m_extents, CGeoExtents::RES_MAX);
        }
    }


    if (msg)
    {
        bool bInitextent = m_options.m_extents.IsRectNull();
        inputDSvector[0].UpdateOption(m_options);

        if (bInitextent)
            m_options.m_extents = extents[1];//used maximum extent. We can add option for that.


        if (!m_options.m_bQuiet)
        {
            CGeoExtents extents = m_options.m_extents;
            cout << "    Extents        = X:{" << to_string(extents.m_xMin) << ", " << to_string(extents.m_xMax) << "}  Y:{" << to_string(extents.m_yMin) << ", " << to_string(extents.m_yMax) << "}" << endl;
            cout << "    Block Size     = " << to_string(extents.m_xBlockSize) << " X " << to_string(extents.m_xBlockSize) << endl << endl;
        }
    }


    if (msg && m_options.m_bCreateImage)
    {
        CImageCalculatorOption options(m_options);
        options.m_scenes_def.clear();
        options.m_nbBands = m_options.m_equations.size();

        if (!m_options.m_bQuiet)
        {
            cout << "Open output images..." << endl;
            cout << "    Size           = " << options.m_extents.m_xSize << " cols x " << options.m_extents.m_ySize << " rows x " << options.m_nbBands << " bands" << endl;
            cout << "    Extents        = X:{" << to_string(options.m_extents.m_xMin) << ", " << to_string(options.m_extents.m_xMax) << "}  Y:{" << to_string(options.m_extents.m_yMin) << ", " << to_string(options.m_extents.m_yMax) << "}" << endl << endl;
        }


        string filePath = options.m_filesPath.back();

        //replace the common part by the new name
        for (size_t zz = 0; zz < m_options.m_equations.size(); zz++)
        {
            string subName = "Eq" + to_string(zz + 1);
            string::size_type pos = m_options.m_equations[zz].find('=');
            if (pos != string::npos)
                subName = m_options.m_equations[zz].substr(0, pos);

//            size_t z = m_options.m_scene_extents[0] + zz;
            options.m_VRTBandsName += GetFileTitle(filePath) + "_" + subName + ".tif|";
        }

        msg += outputDS.CreateImage(filePath, options);
    }

    return msg;
}

void CImageCalculator::ReadBlock(CGDALDatasetExVector& inputDS, int xBlock, int yBlock, deque<CRasterWindow>& windows)
{
    #pragma omp critical(BlockIO)
    {

        m_options.m_timerRead.resume();

        CGeoExtents extents = m_options.m_extents.GetBlockExtents(xBlock, yBlock);


        windows.resize(inputDS.size());
        for (size_t i = 0; i < inputDS.size(); i++)
        {
            inputDS[i].ReadBlock(extents, windows[i], 1, NOT_INIT, NOT_INIT, m_bands_to_read[i]);
        }

        m_options.m_timerRead.stop();

    }
}



void CImageCalculator::ProcessBlock(int xBlock, int yBlock, const deque<CRasterWindow>& windows, OutputData& outputData)
{
    CGeoExtents extents = m_options.GetExtents();
    CGeoSize blockSize = extents.GetBlockSize(xBlock, yBlock);

    if (windows.empty())
    {
        int nbCells = blockSize.m_x * blockSize.m_y;

        #pragma omp atomic
        m_options.m_xx += nbCells;
        m_options.UpdateBar();

        return;
    }



    //init memory
    if (m_options.m_bCreateImage)
    {
        outputData.resize(m_options.m_equations.size());
        for (size_t s = 0; s < outputData.size(); s++)
            outputData[s].insert(outputData[s].begin(), blockSize.m_x * blockSize.m_y, float(m_options.m_dstNodata));
    }




    #pragma omp critical(ProcessBlock)
    {
        m_options.m_timerProcess.resume();

        #pragma omp parallel for num_threads( m_options.m_CPU ) if (m_options.m_bMulti )
        for (int y = 0; y < blockSize.m_y; y++)
        {
            int t = omp_get_thread_num();
            for (int x = 0; x < blockSize.m_x; x++)
            {
                int xy = y * blockSize.m_x + x;

                for (size_t e = 0; e < m_parser[t].size(); e++)
                {
                    const std::deque<CVarDef>& vars = m_parser[t][e].GetParsedVars();

                    vector<double> var_values(vars.size());
                    for (size_t v = 0; v < vars.size(); v++)
                    {
                        size_t i = vars[v].m_ib.m_image;
                        size_t b = vars[v].m_ib.m_band;
                        var_values[v] = windows[i][b].at(x, y);
                    }

                    double result = m_parser[t][e].Evaluate(var_values);
                    outputData[e][xy] = (float)max(GetTypeLimit(m_options.m_outputType, true), min(GetTypeLimit(m_options.m_outputType, false), result));
                }

                #pragma omp atomic
                m_options.m_xx++;

            }//for x
            m_options.UpdateBar();
        }//for y

        m_options.m_timerProcess.stop();

    }//critical process



}


void CImageCalculator::WriteBlock(int xBlock, int yBlock, CGDALDatasetEx& outputDS, OutputData& outputData)
{
    #pragma omp critical(BlockIO)
    {
        m_options.m_timerWrite.resume();

        if (outputDS.IsOpen())
        {
            CGeoExtents extents = outputDS.GetExtents();
            CGeoRectIndex outputRect = extents.GetBlockRect(xBlock, yBlock);

            assert(outputRect.m_x >= 0 && outputRect.m_x < outputDS.GetRasterXSize());
            assert(outputRect.m_y >= 0 && outputRect.m_y < outputDS.GetRasterYSize());
            assert(outputRect.m_xSize > 0 && outputRect.m_xSize <= outputDS.GetRasterXSize());
            assert(outputRect.m_ySize > 0 && outputRect.m_ySize <= outputDS.GetRasterYSize());

            for (size_t z = 0; z < outputData.size(); z++)
            {
                GDALRasterBand* pBand = outputDS.GetRasterBand(z);
                if (!outputData.empty())
                {
                    assert(outputData.size() == outputDS.GetRasterCount());
                    pBand->RasterIO(GF_Write, outputRect.m_x, outputRect.m_y, outputRect.Width(), outputRect.Height(), &(outputData[z][0]), outputRect.Width(), outputRect.Height(), GDT_Float32, 0, 0);
                }
                else
                {
                    double noData = outputDS.GetNoData(z);
                    pBand->RasterIO(GF_Write, outputRect.m_x, outputRect.m_y, outputRect.Width(), outputRect.Height(), &noData, 1, 1, GDT_Float32, 0, 0);
                }
            }
        }

        m_options.m_timerWrite.stop();
    }
}

ERMsg CImageCalculator::CloseAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS)
{
    ERMsg msg;

    for (CGDALDatasetExVector::iterator it = inputDS.begin(); it != inputDS.end(); it++)
        it->Close();



    m_options.m_timerWrite.start();

    msg += outputDS.Close(m_options);

    m_options.m_timerWrite.stop();


    return msg;
}

}
