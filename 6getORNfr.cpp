#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

const double vth=-50;
const double vreset=-70;
const double taum=20;
const double ve=50;
const double vi=-75;
const double vr=-70;
const double gl=1;
const double excitnoise=0.28;
const double inhibnoise=0.5;
const double inputscale=2.5;
const double trefra=2;
const int maxinput=1;
//const double stepsize = 0.0025;
//const int noofstep=maxinput/stepsize;
const double tauadapt=60;
const int noofbisect=9;
const double maxadaptcurrent=40;
const double noofgroupoffile=4;
const int noofodourant = 160;
const int noofodour = 16;
const int noofconc=7;
const int minconc=-6;
const int nooffile = noofconc*2;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double condeff;
	double effinput;
	double efftaum;
	double step = 0;
	double v;
	double adaptcurrent;
	double t;	
	double tupper;
	double tlower;
	double output;
	double input;
	double adaptscale;
	string filename;
	stringstream filenameindex[nooffile];	

	int fileindexcounter=0;
	for (m=0;m<noofconc;m++)
	{
		ifstream infile;
		filenameindex[fileindexcounter]<<"receptorrstar"<<minconc+m<<".txt";
		filename=filenameindex[fileindexcounter].str();
		fileindexcounter=fileindexcounter+1;
		infile.open(filename.c_str());
		ofstream outfile;
		filenameindex[fileindexcounter]<<"ORNfr"<<minconc+m<<".txt";
		filename=filenameindex[fileindexcounter].str();
		fileindexcounter=fileindexcounter+1;
		outfile.open(filename.c_str());	
		for (j=0;j<noofodour;j++)
		{
			for (k=0;k<noofodourant;k++)
			{
				infile>>input;
				condeff=1/(1+excitnoise+inhibnoise+input*inputscale);
				effinput=(ve*(input*inputscale+excitnoise)+vi*inhibnoise+vr*gl)*condeff;
				efftaum=taum*condeff;
				v=vreset;
				adaptcurrent=maxadaptcurrent*sqrt(input);
				t=0;
				if (effinput<=vth)
					output=0;
				else
				{
					adaptscale = tauadapt*adaptcurrent/(tauadapt-efftaum);
					while (v<vth)
					{
						t=t+1;
						v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
					}
					tupper=t;
					tlower=t-1;
					t=(tupper+tlower)*0.5;
					for(i=0;i<noofbisect;i++)
					{
						v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
						if (v>vth)
							tupper=t;
						else
							tlower=t;
						t=(tupper+tlower)*0.5;								
					}
					output=1000/(t+trefra);				
				}
				outfile<<output<<"	";
			}
			outfile<<endl;
		}
		infile.close();
		outfile.close();
	}
}
