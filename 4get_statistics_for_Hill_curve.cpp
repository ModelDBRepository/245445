#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;
#include "RandNum.h"
CRandNum randNum;

#define RAN 	randNum.GenRandReal_10()
#define RAN00 	randNum.GenRandReal_00()
#define RANDINT randNum.GenRandInt32()
#define Norm	randNum.Gaussian_Noise()

const int noofodour = 16;
const int noofodourant = 160;
const double chalfmax= -0.4;
const double chalfmin = -4.2;
const double chalfmean = -2.8;
const double chalfsd = 1;
const double slopemin = 0.6;
const double slopemax = 2.8;
const double logslopemean = (log(slopemin)+log(slopemax))*0.5;
const double logslopesd = 0.25;
//const double corrnfmax = 0.2;
const double assumedataconc = 0;
const int minconc=-6;
const int noofconc = 7;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double chalf[noofodourant][noofodour];
	double fmax[noofodourant][noofodour];
	double fmaxmean = 0;
	double maxfmax=-1000;
	double fdata[noofodourant][noofodour];
	double slope[noofodourant];
	string filename;
	stringstream filenameindex[noofconc];
	srand(time(NULL));
				
	ifstream inputfile;
	inputfile.open("chemsimresponse.txt");
	for (j=0;j<noofodour;j++)
		for (i=0;i<noofodourant;i++)
			inputfile>>fdata[i][j];
	inputfile.close();	
	
	//find mean of response
	for (i=0;i<noofodourant;i++)
		for (j=0;j<noofodour;j++)
			fmaxmean = fmaxmean+fdata[i][j];
	fmaxmean=fmaxmean/(noofodourant*noofodour);
		
	for (i=0;i<noofodourant;i++)
	{
		for (j=0;j<noofodour;j++)
		{
		//find chalf slope and temporary value for fmax 
			do
				chalf[i][j]=chalfmean+chalfsd*Norm;		
			while ((chalf[i][j]>chalfmax) or (chalf[i][j]<chalfmin));
			fmax[i][j]= fdata[i][j];
		}
		do
			slope[i]=exp(logslopemean+logslopesd*Norm);
		while ((slope[i]>slopemax) or (slope[i]<slopemin));		
	}
		
	//infer new fmax from data/older value and find new mean for fmax
	for (i=0;i<noofodourant;i++)
		for (j=0;j<noofodour;j++)
		{
			fmax[i][j]=fmax[i][j]*(1+exp(-slope[i]*(assumedataconc-chalf[i][j])));
			if (fmax[i][j]>maxfmax)
				maxfmax=fmax[i][j];				
		}
	
	ofstream outone;
	outone.open("fmax.txt");
	for (j=0;j<noofodour;j++)
	{
		for (k=0;k<noofodourant;k++)
			outone<<fmax[k][j]<<"	";
		outone<<endl;
	}	
	outone.close(); 
	
	ofstream outtwo;
	outtwo.open("nprime.txt");

	for (k=0;k<noofodourant;k++)
		outtwo<<slope[k]<<"	";
	outtwo.close(); 
	
	ofstream outthree;
	outthree.open("chalf.txt");
	for (j=0;j<noofodour;j++)
	{
		for (k=0;k<noofodourant;k++)
			outthree<<chalf[k][j]<<"	";
		outthree<<endl;
	}	
	outthree.close(); 
	

	for (m=0;m<noofconc;m++)
	{
		ofstream outfour;
		filenameindex[m]<<"receptorrstar"<<minconc+m<<".txt";
		filename = filenameindex[m].str(); 		
		outfour.open(filename.c_str());
		for (j=0;j<noofodour;j++)
		{
			for (k=0;k<noofodourant;k++)
				outfour<<fmax[k][j]/(1+exp(-slope[k]*(double(minconc+m)-chalf[k][j])))<<"	";
			outfour<<endl;
		}	
		outfour.close(); 
	}
	
	ofstream outsix;
	outsix.open("maxfmax.txt");
	outsix<<maxfmax<<endl;
	outsix.close();	
return 1;	
}
