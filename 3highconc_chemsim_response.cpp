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
const int knownodourant = 26;
const double whochangemagicnoone = 0.2;
const double whochangemagicnotwo = 0.8;
const double whochangemagicnothree = 0.5;
const double whochangemagicnofour = 0.25;
const double scaleresmagicno = 1;
const double ratioofsimchange = 1;
const double smallchange = 0.05;
const double smallres = 0.02;
const int totalnoit = 100;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double defchemsim[noofodour][noofodour-1];
	double responsesim[noofodour][noofodour-1];
	double response[noofodourant][noofodour];
	double scalefactor;
	double chemsimdiff;
	double maxchemsim;
	double perturbone;
	double perturbtwo;
	double tempone;
	double temptwo;
			
	ifstream inputfile;
	inputfile.open("ORNdatachemsim.txt");
	maxchemsim=0;
	for (j=1;j<noofodour;j++)
		for (m=0;m<j;m++)
			{
				inputfile>>defchemsim[j][m];
				if (defchemsim[j][m]>maxchemsim)
					maxchemsim=	defchemsim[j][m];
			}

	inputfile.close();	
	
//	ifstream inputresponsesim;
//	inputfile.open("alldatachemsim.txt");
//	for (j=1;j<noofodour;j++)
//		for (m=0;m<j;m++)
//		inputfile>>responsesim[j][m];
//	inputfile.close();	
	
	ifstream inputresponse;	
	inputresponse.open("responseinter.txt");
	for (j=0;j<noofodour;j++)
		for (i=0;i<noofodourant;i++)
			inputresponse>>response[i][j];
	inputresponse.close();	
	
	for (k=1;k<=totalnoit;k++)
	{
		scalefactor=(double)(1/(double)(k));            //less changes made in later iteration
		for (j=1;j<noofodour;j++)  
			for (m=0;m<j;m++)
			{
				responsesim[j][m]=0;                   //finding similarity of odor in the current responses
				for (i=knownodourant;i<noofodourant;i++)
					responsesim[j][m]=responsesim[j][m]+(response[i][j]-response[i][m])*(response[i][j]-response[i][m]);
				responsesim[j][m]=sqrt(responsesim[j][m]/noofodourant);
			}		
		
		for (j=1;j<noofodour;j++)
		{
			for (m=0;m<j;m++)
			{
				chemsimdiff=defchemsim[j][m]-responsesim[j][m];
				for (i=knownodourant;i<noofodourant;i++)
				{
					tempone=RAN;
					temptwo=RAN;
					
					if (chemsimdiff<0)
					{
						perturbone = scaleresmagicno*scalefactor*ratioofsimchange*temptwo*(-sqrt(-chemsimdiff))*(maxchemsim-defchemsim[j][m]);
						perturbtwo = scaleresmagicno*scalefactor*ratioofsimchange*temptwo*(-sqrt(-chemsimdiff))*(maxchemsim-defchemsim[j][m]);	  //determine size of changes				
					}
					else
					{
						perturbone = scaleresmagicno*scalefactor*temptwo*(sqrt(chemsimdiff)*defchemsim[j][m]);
						perturbtwo = scaleresmagicno*scalefactor*temptwo*(sqrt(chemsimdiff)*defchemsim[j][m]);							
					}
					if (response[i][j]>=response[i][m])
					{
						if (tempone<whochangemagicnoone)
							response[i][j]=response[i][j]+perturbone+perturbtwo;
						else if (tempone>=whochangemagicnotwo)
							response[i][m]=response[i][m]-perturbone-perturbtwo;
						else
						{
							response[i][j]=response[i][j]+perturbone;
							response[i][m]=response[i][m]-perturbtwo;
						}
					}
					if (response[i][m]>response[i][j])
					{
						if (tempone<whochangemagicnoone)
							response[i][j]=response[i][j]-perturbone-perturbtwo;
						else if (tempone>=whochangemagicnotwo)
							response[i][m]=response[i][m]+perturbone+perturbtwo;
						else
						{
							response[i][j]=response[i][j]-perturbone;
							response[i][m]=response[i][m]+perturbtwo;
						}
					}						
			
				}				
			}		
		}
	}

	for (i=knownodourant;i<noofodourant;i++)
		for (j=0;j<noofodour;j++)
			if (response[i][j]<smallres)
				response[i][j]=smallres*RAN;           //reset negative responses to a very small positive value

	ofstream outres;
	outres.open("chemsimresponse.txt");
	for (j=0;j<noofodour;j++)
	{
		for (i=0;i<noofodourant;i++)
			outres<<response[i][j]<<"	";
		outres<<endl;
	}	
	outres.close(); 
	
	for (j=1;j<noofodour;j++)  
		for (m=0;m<j;m++)
		{
			responsesim[j][m]=0;
			for (i=knownodourant;i<noofodourant;i++)
				responsesim[j][m]=responsesim[j][m]+(response[i][j]-response[i][m])*(response[i][j]-response[i][m]);
			responsesim[j][m]=sqrt(responsesim[j][m]/noofodourant);
		}		
		
	ofstream outchemsim;
	outchemsim.open("newchemsim.txt");
	for (j=1;j<noofodour;j++)
		{
		for (m=0;m<j;m++)
			outchemsim<<responsesim[j][m]<<"	";
		outchemsim<<endl;
		}
		outchemsim.close(); 				
			
return 1;	
}
