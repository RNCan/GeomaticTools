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


//"D:\Travaux\Landsat\Landsat(2000-2018)\Input\Landsat_2000-2018(2).vrt" "D:\Travaux\Landsat\Landsat(2000-2018)\Output\test3.vrt" -of VRT -overwrite -co "COMPRESS=LZW"   -te 1022538.9 6663106.0 1040929.5 6676670.7 -multi -SpikeThreshold 0.75


//#define BOOST_NO_CXX11_SCOPED_ENUMS

//#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include <cmath>
#include <array>
#include <utility>
#include <iostream>


#include <boost/filesystem.hpp>

//#include <boost/filesystem.hpp>


#include "basic/OpenMP.h"
#include "geomatic/UtilGDAL.h"


#include "ImageCalculator.h"
//#include "geomatic/LandTrendUtil.h"
//#include "geomatic/LandTrendCore.h"





using namespace std;
//using namespace WBSF::Landsat2;
//using namespace LTR;



namespace WBSF
{
	const char* CImageCalculator::VERSION = "1.0.0";
	const size_t CImageCalculator::NB_THREAD_PROCESS = 2;


	//*********************************************************************************************************************

	CImageCalculatorOption::CImageCalculatorOption()
	{
			
		m_appDescription = "This software compute bands from many inputs images into destination image using user define math equation. This Software use MTParser (Mathieu Jacques) to evaluate equation.";
		std::string indicesName = Landsat2::GetIndiceNames();

		//AddOption("-RGB");
		static const COptionDef OPTIONS[] =
		{
			{"-Equation",1,"\"formulas\"",true,"Equation to evaluate. One equation per output band. Variables in equation must have the form \"Name=i#b#\" where Name is the nameof the output variable (use when output is VRT), # is 1 to the number of images or bands."},
			{ "-e", 1, "\"formulas\"", true, "Equilvalent to -Equation." },
			{"srcfile",0,"",true, "Input images file path (no limits)."},
			{"dstfile",0,"",false, "Output image file path."}
		};

//		AddOption("-ty");
		for (size_t i = 0; i < sizeof(OPTIONS) / sizeof(COptionDef); i++)
			AddOption(OPTIONS[i]);


//        boost::filesystem::path p1 = "c:/test/a.txt";
  //      boost::filesystem::path p2 = "c:/test/a.txt";
    //    boost::filesystem::absolute(p1, p2);
		//Pour les trigger Bande 1 c’est - 125 quand on fait  ex.b1 1994 – b1 1995 ou b1 1996 – b1 1995.
		//Pour le tassel Cap brightness c’est + 750  ex.tcb1994 – tcb 1995 ou tcb 1996 – tcb 1995


		static const CIOFileInfoDef IO_FILE_INFO[] =
		{
			{"Input image","srcfile","1","*","","Multi-layer input"},
			{"Output Image", "dstfile","1","*","","Same as input file"},
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
				msg.ajoute("ERROR: Invalid argument line. At least 2 images is needed: input image and destination image.\n");
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
			string name;
			string equation = argv[++i];
			string::size_type pos = equation.find('=');
			if (pos != string::npos)
			{
				name = equation.substr(0, pos);
				equation = equation.substr(pos + 1);
			}
			else
			{
				name = "Equation" + to_string(m_equations.size() + 1);
			}

			m_equations.push_back({ name,equation });
			//m_name.push_back(name);
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




		if (!m_options.m_bQuiet)
		{
			cout << "Output: " << m_options.m_filesPath[CImageCalculatorOption::OUTPUT_FILE_PATH] << endl;
			cout << "From:   " << m_options.m_filesPath[CImageCalculatorOption::INPUT_FILE_PATH] << endl;

			if (!m_options.m_maskName.empty())
				cout << "Mask:   " << m_options.m_maskName << endl;

		}


		/*m_options.m_modifier = -1;
		if (m_options.m_indice == Landsat2::I_B4)
			m_options.m_modifier = 1;*/


		GDALAllRegister();

		CGDALDatasetExVector inputDS;
		//CGDALDatasetEx maskDS;
		CGDALDatasetEx outputDS;
		

		msg = OpenAll(inputDS, outputDS);
		if (!msg)
			return msg;


		if (!m_options.m_bQuiet && m_options.m_bCreateImage)
			cout << "Create output images (" << outputDS.GetRasterXSize() << " C x " << outputDS.GetRasterYSize() << " R x " << outputDS.GetRasterCount() << " B) with " << m_options.m_CPU << " threads..." << endl;

		CGeoExtents extents = m_options.m_extents;
		m_options.ResetBar((size_t)extents.m_xSize * extents.m_ySize);

		vector<pair<int, int>> XYindex = extents.GetBlockList();
		map<int, bool> treadNo;

		//omp_set_nested(1);//for IOCPU
		//#pragma omp parallel for schedule(static, 1) num_threads( NB_THREAD_PROCESS ) if (m_options.m_bMulti)
		for (int b = 0; b < (int)XYindex.size(); b++)
		{
			int xBlock = XYindex[b].first;
			int yBlock = XYindex[b].second;

			Landsat2::CLandsatWindow inputData;
			OutputData outputData;

			ReadBlock(inputDS, xBlock, yBlock, inputData);
			ProcessBlock(xBlock, yBlock, inputData, outputData);
			WriteBlock(xBlock, yBlock, outputDS, outputData);
		}//for all blocks

		CloseAll(inputDS, outputDS);

		return msg;
	}



	ERMsg CImageCalculator::OpenAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS)
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


		inputDS.resize(filePathIn.size());
		for (size_t i = 0; i < filePathIn.size() && msg; i++)
			msg += inputDS[i].OpenInputImage(filePathIn[i], m_options);


		if (msg)
		{
			inputDS[0].UpdateOption(m_options);

			if (!m_options.m_bQuiet)
			{
				CGeoExtents extents = inputDS[0].GetExtents();
				//            CProjectionPtr pPrj = inputDS.GetPrj();
				//            string prjName = pPrj ? pPrj->GetName() : "Unknown";

				cout << "    Size           = " << inputDS->GetRasterXSize() << " cols x " << inputDS->GetRasterYSize() << " rows x " << inputDS.GetRasterCount() << " bands" << endl;
				cout << "    Extents        = X:{" << to_string(extents.m_xMin) << ", " << to_string(extents.m_xMax) << "}  Y:{" << to_string(extents.m_yMin) << ", " << to_string(extents.m_yMax) << "}" << endl;
				//          cout << "    Projection     = " << prjName << endl;
				cout << "    NbBands        = " << inputDS.GetRasterCount() << endl;
			}
		}


		if (msg && m_options.m_bCreateImage)
		{
			CImageCalculatorOption options(m_options);
			options.m_scenes_def.clear();
			options.m_nbBands = m_options.m_equations.size();

			if (!m_options.m_bQuiet)
			{
				cout << endl;
				cout << "Open output images..." << endl;
				cout << "    Size           = " << options.m_extents.m_xSize << " cols x " << options.m_extents.m_ySize << " rows x " << options.m_nbBands << " bands" << endl;
				cout << "    Extents        = X:{" << to_string(options.m_extents.m_xMin) << ", " << to_string(options.m_extents.m_xMax) << "}  Y:{" << to_string(options.m_extents.m_yMin) << ", " << to_string(options.m_extents.m_yMax) << "}" << endl;
				//cout << "    NbBands        = " << options.m_nbBands << endl;
				cout << "    Nb. Scenes     = " << nb_scenes << endl;
			}

			//std::string indices_name = Landsat2::GetIndiceName(m_options.m_indice);
			string filePath = options.m_filesPath[CImageCalculatorOption::OUTPUT_FILE_PATH];

			//replace the common part by the new name
			for (size_t zz = 0; zz < m_options.m_equations.size(); zz++)
			{
				size_t z = m_options.m_scene_extents[0] + zz;
				string subName = m_options.m_equations[zz].first;
				options.m_VRTBandsName += GetFileTitle(filePath) + "_" + subName + ".tif|";
			}

			msg += outputDS.CreateImage(filePath, options);
		}

		return msg;
	}

	void CImageCalculator::ReadBlock(CGDALDatasetExVector& inputDS, int xBlock, int yBlock, CRasterWindow& block_data)
	{
#pragma omp critical(BlockIO)
		{
			m_options.m_timerRead.start();

			CGeoExtents extents = m_options.m_extents.GetBlockExtents(xBlock, yBlock);
			inputDS.ReadBlock(extents, block_data, 1, m_options.m_IOCPU);

			m_options.m_timerRead.stop();



			const vector<string>& formulas = m_options.m_equations;

			//create virtual band
			unordered_set<string> varList;

			m_equations.clear();
			m_equations.resize(formulas.size());
			for (size_t i = 0; i < formulas.size() && msg; i++)
			{
				m_equations[i].m_bManageSrcNoData = option.m_bManageSrcNoData;
				m_equations[i].m_formula = formulas[i];
				msg += m_equations[i].Compile(inputDSVector, option.m_CPU);

				if (msg)
				{
					const CMTParserVariableVector& vars = m_equations[i].GetVars();
					for (CMTParserVariableVector::const_iterator it = vars.begin(); it != vars.end(); it++)
						varList.insert(it->m_name);
				}
			}

			if (!msg)
				return msg;

			//now create the link variable bandHolder
			m_imageBandToData.clear();
			m_variableToData.clear();

			m_variableToData.resize(m_equations.size());
			for (size_t i = 0; i < m_equations.size(); i++)
			{
				const CMTParserVariableVector& var = m_equations[i].GetVars();

				m_variableToData[i].resize(var.size());
				for (size_t j = 0; j < m_variableToData[i].size() && msg; j++)
				{
					CImageBandPos IBpos = GetImageBand(var[j]);

					if (msg)
					{
						vector< pair<size_t, size_t> >::const_iterator item = std::search_n(m_imageBandToData.begin(), m_imageBandToData.end(), 1, IBpos);

						if (item == m_imageBandToData.end())
						{
							//not already in the list
							m_variableToData[i][j] = (int)m_imageBandToData.size();
							m_imageBandToData.push_back(IBpos);
						}
						else
						{
							//already in the list
							m_variableToData[i][j] = int(item - m_imageBandToData.begin());
						}
					}
				}
			}

			//now, add band to band holder
			for (CImageBandPosVector::const_iterator it = m_imageBandToData.begin(); it != m_imageBandToData.end(); it++)
			{
				CSingleBandHolderPtr pBandHolder = inputDSVector[it->first].GetSingleBandHolder(it->second);
				AddBand(pBandHolder);
			}

			m_entireExtents = option.GetExtents();

			for (int i = 0; i < (int)m_bandHolder.size(); i++)
				msg += m_bandHolder[i]->Load(this);


			if (!option.m_bQuiet)
			{
				cout << "Memory block: " << m_entireExtents.XNbBlocks() << " Blocks per line x " << m_entireExtents.YNbBlocks() << " Lines" << endl;
				cout << "Memory block size: " << GetBlockSizeX() << " Cols x " << GetBlockSizeY() << " Rows x " << GetBlockSizeZ() << " Bands (max)" << endl;
			}

		}
	}



	void CImageCalculator::ProcessBlock(int xBlock, int yBlock, const CRasterWindow& window, OutputData& outputData)
	{
		CGeoExtents extents = m_options.GetExtents();
		CGeoSize blockSize = extents.GetBlockSize(xBlock, yBlock);

		if (window.empty())
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
			outputData.resize(window.GetNbScenes());
			for (size_t s = 0; s < outputData.size(); s++)
				outputData[s].insert(outputData[s].begin(), blockSize.m_x * blockSize.m_y, m_options.m_dstNodata);
		}



#pragma omp critical(ProcessBlock)
		{
			m_options.m_timerProcess.start();


			//#pragma omp parallel for num_threads( m_options.m_CPU ) if (m_options.m_bMulti )
			for (int y = 0; y < blockSize.m_y; y++)
			{
				for (int x = 0; x < blockSize.m_x; x++)
				{
					int xy = y * blockSize.m_x + x;

					//Get pixel
					const CImageBandPosVector& imageBandToData = bandHolder.GetImageBandToData();
					const CVariableToData& variableToData = bandHolder.GetVariableToData();


					int nbVariables = (int)imageBandToData.size();
					int nbFormulas = (int)bandHolder.GetVirtualBandVector().size();


					for (size_t z = 0; z < window.size(); z++)
					{
						size_t zz = z;
						data[z] = window.getpix.GetPixelIndice(zz, m_options.m_indice, x, y, m_options.m_rings);
					}

					//for all equation (band)
					for (size_t z = 0; z < variableToData.size(); z++)
					{
						vector<float> vars(variableToData[z].size());
						for (size_t v = 0; v < variableToData[z].size(); v++)
						{
							size_t pos = variableToData[z][v];
							vars[v] = input[pos]->at(x, y);
						}

						double value = bandHolder.Evaluate((int)z, vars);
						output[z][y][x] = (float)outputDS.PostTreatment(value);
					}

					assert(all_y.size() == data3.size());
					double val = max(GetTypeLimit(m_options.m_outputType, true), min(GetTypeLimit(m_options.m_outputType, false), output_corr_factore[1]));
					outputData[z][xy] = val;

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
			m_options.m_timerWrite.start();

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
						pBand->RasterIO(GF_Write, outputRect.m_x, outputRect.m_y, outputRect.Width(), outputRect.Height(), &(outputData[z][0]), outputRect.Width(), outputRect.Height(), GDT_Float64, 0, 0);
					}
					else
					{
						double noData = outputDS.GetNoData(z);
						pBand->RasterIO(GF_Write, outputRect.m_x, outputRect.m_y, outputRect.Width(), outputRect.Height(), &noData, 1, 1, GDT_Float64, 0, 0);
					}
				}
			}

			m_options.m_timerWrite.stop();
		}
	}

	void CImageCalculator::CloseAll(CGDALDatasetExVector& inputDS, CGDALDatasetEx& outputDS)
	{
		for (CGDALDatasetExVector::iterator it = inputDS.begin(); it != inputDS.end(); it++)
			it->Close();
		
		

		m_options.m_timerWrite.start();

		outputDS.Close(m_options);

		m_options.m_timerWrite.stop();

	}

}
