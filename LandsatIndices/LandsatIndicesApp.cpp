//***********************************************************************
// program to merge Landsat image image over a period
//									 
//***********************************************************************


#include <boost/timer/timer.hpp>
#include <stdio.h>
#include <iostream>
#include "LandsatIndices.h"



using namespace std;
using namespace WBSF;

//***********************************************************************
//									 
//	Main                                                             
//						 	 		 
//***********************************************************************
int main(int argc, char* argv[])
{
    boost::timer::cpu_timer timer;
    timer.start();

    //Create a mergeImages object
    CLandsatIndices app;

    ERMsg msg = app.m_options.ParseOption(argc, argv);

    if (!app.m_options.m_bQuiet)
        cout << app.GetDescription() << endl;


    if (msg)
        msg = app.Execute();

    if (!msg)
    {
        PrintMessage(msg);
        return -1;
    }


    timer.stop();

    if (!app.m_options.m_bQuiet)
    {
        app.m_options.PrintTime();
        cout << endl << "Total time = " << SecondToDHMS(timer.elapsed().wall / 1e9) << endl;
    }



    return 0;
}
