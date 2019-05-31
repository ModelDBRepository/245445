#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <sstream>
using namespace std;
#include "RandNum.h"
CRandNum randNum;

#define RAN 	randNum.GenRandReal_10()
#define RAN00 	randNum.GenRandReal_00()
#define RANDINT randNum.GenRandInt32()
#define Norm	randNum.Gaussian_Noise()

const double vth=-50;
const double vreset=-70;
const double taum=20;
const double ve=50;
const double vi=-75;
const double vr=-70;
const double gl=1;
const double PNconstexcitconmean = 0.24;
const double PNconstinhibconmean = 0.15;
const double PNconstexcitconvar = 0.02;
const double PNconstinhibconvar = 0.01;
const double inputscalePNconst=0.045;
const double inputscaleLNconst=0.013;
const double inputscaleinterPNconst=0.04;
const double inputscaleinterLNconst=0.004;
const double otherinhibconst=1;
const double trefra=2;
const int noofodourant=160;
const int noofodour=16;
const double adaptcurrentmaxPNmean=4.5;
const double adaptcurrentmaxLNmean=1.8;
const double adaptcurrentmaxPNvar=0.4;
const double adaptcurrentmaxLNvar=0.2;
const double tauadapt=25;
const int noofbisect=9;
const int noofiteration=100;
const double corrconnratio = 0.01;
const double noiseconnratio = 0.001;
const double baseconn = 0.006;
const double baseconnnoise = 0.002;
const double noisePN = 0.01;
const double noiseLN = 0.003;
const double noiseinterPN = 0.01;
const double noiseinterLN = 0.001;
const int noofconc=7;
const int minconc=-6;
const int nooffile = noofconc*2;

int main()
{
	int i;
	int j;
	int m;
	int k;
	int n;
	int fileindexcounter;
	double condeff;
	double temp;
	double tempone;
	double temptwo;
	double tempthree;
	double tempfour;
	double tempfive;
	double effinput;
	double efftaum;
	double rantempone;
	double rantemptwo;
	double adaptcurrentmaxPN[noofodourant];
	double adaptcurrentmaxLN[noofodourant];
	double ORNmaxresponse[noofodour][noofodourant];
	double interconnectivity[noofodourant][noofodourant];
	double rawinput[noofodourant];
	double outputPN[noofodourant];
	double outputLN[noofodourant];
	double outputLNold[noofodourant];
	double inputscalePN[noofodourant];
	double inputscaleLN[noofodourant];
	double inputscaleinterPN[noofodourant];
	double inputscaleinterLN[noofodourant];
	double PNconstexcitcon[noofodourant];
	double PNconstinhibcon[noofodourant];
	double otherinhibinput;
	double v;
	double t;	
	double tupper;
	double tlower;
	double adaptscale;
	string filename;
	stringstream filenameindex[nooffile];

	srand(time(NULL));
	
	for (i=0;i<noofodourant;i++)
	{
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		inputscalePN[i] = inputscalePNconst+noisePN*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		inputscaleLN[i] = inputscaleLNconst+noiseLN*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		inputscaleinterPN[i] = inputscaleinterPNconst+noiseinterPN*rantempone;		
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		inputscaleinterLN[i] = inputscaleinterLNconst+noiseinterLN*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		PNconstexcitcon[i] = PNconstexcitconmean+PNconstexcitconvar*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		PNconstinhibcon[i] = PNconstinhibconmean+PNconstinhibconvar*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		adaptcurrentmaxPN[i] = adaptcurrentmaxPNmean+adaptcurrentmaxPNvar*rantempone;
		do	{rantempone = Norm;}
		while (abs(rantempone)>2);
		adaptcurrentmaxLN[i] = adaptcurrentmaxLNmean+adaptcurrentmaxLNvar*rantempone;				
	}
	
	//find connectivity between LN and PN+other LN
	for (i=0;i<noofodourant;i++)
		interconnectivity[i][i]=0; 

	ifstream inputone;
	inputone.open("ORNfr0.txt");
	
	for (k=0;k<noofodour;k++)
		for (j=0;j<noofodourant;j++)
			inputone>>ORNmaxresponse[k][j];
	inputone.close();		

	//calculate correlation matrix for ORN-PN connectivity
	ofstream outtest;
	outtest.open("corrmatrix1.txt");
	for (j=1;j<noofodourant;j++)
	{
		for (i=0;i<j;i++)
		{	
			tempone=0;
			temptwo=0;
			tempthree=0;
			tempfour=0;
			tempfive=0;
			for (k=0;k<noofodour;k++)
			{
				tempone=tempone+ORNmaxresponse[k][i];
				temptwo=temptwo+ORNmaxresponse[k][j];
				tempthree=tempthree+ORNmaxresponse[k][i]*ORNmaxresponse[k][i];
				tempfour=tempfour+ORNmaxresponse[k][j]*ORNmaxresponse[k][j];
				tempfive=tempfive+ORNmaxresponse[k][i]*ORNmaxresponse[k][j];						
			}
			tempone=tempone/(double)noofodour;
			temptwo=temptwo/(double)noofodour;
			tempthree=tempthree/(double)noofodour;
			tempfour=tempfour/(double)noofodour;
			tempfive=tempfive/(double)noofodour;
			temp=(tempfive-tempone*temptwo)/(sqrt((tempthree-tempone*tempone)*(tempfour-temptwo*temptwo)));	
			outtest<<temp<<"	";
			if (temp<0)
			{
				do	{rantempone = Norm;}
				while (abs(rantempone)>2);
				interconnectivity[i][j]=baseconn+baseconnnoise*rantempone;
				do	{rantempone = Norm;}
				while (abs(rantempone)>2);
				interconnectivity[j][i]=baseconn+baseconnnoise*rantempone;					
			}
			else
			{
				do	{rantempone = Norm;}
				while (abs(rantempone)>2);
				do	{rantemptwo = Norm;}
				while (abs(rantemptwo)>2);
				interconnectivity[i][j]=baseconn+baseconnnoise*rantempone+temp*(corrconnratio+noiseconnratio*rantemptwo);
				do	{rantempone = Norm;}
				while (abs(rantempone)>2);
				do	{rantemptwo = Norm;}
				while (abs(rantemptwo)>2);
				interconnectivity[j][i]=baseconn+baseconnnoise*rantempone+temp*(corrconnratio+noiseconnratio*rantemptwo);					
			}
		}							
	outtest<<endl;
	}
	outtest.close();
	fileindexcounter=0;
	for (m=0;m<noofconc;m++)
	{
		ifstream infile;
		filenameindex[fileindexcounter]<<"ORNfr"<<minconc+m<<".txt";
		filename=filenameindex[fileindexcounter].str();
		fileindexcounter=fileindexcounter+1;
		infile.open(filename.c_str());
		ofstream outfile;
		filenameindex[fileindexcounter]<<"PNfr"<<minconc+m<<".txt";
		filename=filenameindex[fileindexcounter].str();
		fileindexcounter=fileindexcounter+1;
		outfile.open(filename.c_str());				
		for (j=0;j<noofodour;j++)
		{
			for (i=0;i<noofodourant;i++)
				infile>>rawinput[i];
		
			for (i=0;i<noofodourant;i++)
				outputLNold[i]=0;
		
		//iterate until steady state
			for (n=0;n<noofiteration;n++)
			{		
				for (i=0;i<noofodourant;i++)
				{
				//find total (raw) inhibitory input received 		
				
					otherinhibinput=0;
					for	(k=0;k<noofodourant;k++)
						otherinhibinput=otherinhibinput+interconnectivity[i][k]*outputLNold[k];	
					otherinhibinput=otherinhibinput*otherinhibconst;

				//find response of PN		
					condeff=1/(1+(rawinput[i]*inputscalePN[i]+PNconstexcitcon[i]+PNconstinhibcon[i]+inputscaleinterPN[i]*otherinhibinput));
					effinput=(ve*(rawinput[i]*inputscalePN[i]+PNconstexcitcon[i])+vi*(inputscaleinterPN[i]*otherinhibinput+PNconstinhibcon[i])+vr*gl)*condeff;
					efftaum=taum*condeff;
					v=vreset;
					t=0;
					if (effinput<=vth)
						outputPN[i]=0;
					else
					{
						adaptscale = tauadapt*adaptcurrentmaxPN[i]*(sqrt(outputPN[i]))/(tauadapt-efftaum);
						while (v<vth)
						{
							t=t+1;
							v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
						}
						tupper=t;
						tlower=t-1;
						t=(tupper+tlower)*0.5;
						for(k=0;k<noofbisect;k++)
						{
							v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
							if (v>vth)
								tupper=t;
							else
							tlower=t;
							t=(tupper+tlower)*0.5;								
						}
						outputPN[i]=1000/(t+trefra);				
					}
		
				//find response of LN				
					condeff=1/(1+(rawinput[i]*inputscaleLN[i]+inputscaleinterLN[i]*otherinhibinput));
					effinput=(ve*rawinput[i]*inputscaleLN[i]+vi*inputscaleinterLN[i]*otherinhibinput+vr*gl)*condeff;
					efftaum=taum*condeff;
					v=vreset;
					t=0;
					if (effinput<=vth)
						outputLN[i]=0;
					else
					{
						adaptscale = tauadapt*adaptcurrentmaxLN[i]*(sqrt(outputLN[i]))/(tauadapt-efftaum);
						while (v<vth)
						{
							t=t+1;
							v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
						}
						tupper=t;
						tlower=t-1;
						t=(tupper+tlower)*0.5;
						for(k=0;k<noofbisect;k++)
						{
							v=vreset*exp(-t/efftaum)+effinput*(1-exp(-t/efftaum))-adaptscale*(exp(-t/tauadapt)-exp(-t/efftaum));
							if (v>vth)
								tupper=t;
							else
								tlower=t;
							t=(tupper+tlower)*0.5;								
						}
						outputLN[i]=1000/(t+trefra);						
					}
				}
				for (i=0;i<noofodourant;i++)
					outputLNold[i]=outputLN[i];			
			}
	
			for (i=0;i<noofodourant;i++)
				outfile<<outputPN[i]<<"	";
			outfile<<endl;	
		}
		infile.close();		
		outfile.close(); 	
	}
}
