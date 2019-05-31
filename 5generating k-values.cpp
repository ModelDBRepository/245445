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

const double ktwomean = 0.1;
const double ktwosd = 0.01;
const double konemean = 3;
const double konesd = 0.6;
const double rzero = 2;
const int noofodour = 16;
const int noofconc = 1;
const int noofodourant = 160;
const int noofoutputfileodour=16;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double kone;
	double knegone;
	double ktwo;
	double knegtwo;
	double n[noofodourant];
	double chalf;
	double maxfr[noofodourant][noofconc];
	double maxcond[noofodourant][noofconc];
	double fmax;
	double maxfmax;
	int counter;
	double temp;
	int badnesscounter;
	int totalbad=0;
	int totalbadtwo=0;
	int totalbadthree=0;
	int counterfilename;
	int on;
	srand(time(NULL));

	//generate k-values with parameters from Hill curves as contraints
	ifstream inputfileone;
	ifstream inputfiletwo;
	ifstream inputfilethree;
	ifstream inputfilefour;
	inputfileone.open("chalf.txt");
	inputfilethree.open("fmax.txt");
	inputfilefour.open("maxfmax.txt");
	inputfilefour>>maxfmax;
	inputfilefour.close();

	ofstream outputone;
	outputone.open("kone.txt");
	ofstream outputtwo;
	outputtwo.open("ktwo.txt");
	ofstream outputthree;
	outputthree.open("knegone.txt");
	ofstream outputfour;
	outputfour.open("knegtwo.txt");	
	ofstream outputfive;
	outputfive.open("totalbad.txt");
//	ofstream outputsix;
//	outputsix.open("thalfrise.txt");
//	ofstream outputseven;
//	outputseven.open("thalfdecay.txt");	
//	ofstream outputeight;
//	outputeight.open("tonetenthrise.txt");
//	ofstream outputnine;
//	outputnine.open("tonetenthdecay.txt");

	
	totalbad=0;
	totalbadtwo=0;
	
	inputfiletwo.open("nprime.txt");	
	for (i=0;i<noofodourant;i++)
	{
		inputfiletwo>>n[i];
		n[i]=n[i]/log(10);
	}	
	inputfiletwo.close();
	
	for (j=0;j<noofodour;j++)
	{
	
		for (i=0;i<noofodourant;i++)
		{
			inputfileone>>chalf;	
			inputfilethree>>fmax;
		

			ktwo=ktwomean+Norm*ktwosd;
			if (fmax>0.02)
				knegtwo=ktwo*(rzero/fmax-1);
			else knegtwo=50*(1+0.1*Norm);
			badnesscounter=0;
			do
			{
				kone = konemean+Norm*konesd;
				temp= n[i]*chalf*log(10);
				kone=kone*(1+0.25*Norm)/sqrt(exp(temp));
				knegone=kone*exp(temp)*(1+ktwo/knegtwo);
				if (knegone<0.01)
					{
						knegone=0.01*(1+0.1*RAN);
						if ((kone<5000) and (kone>0.1))
							totalbadtwo=totalbadtwo+1;
						else
							totalbadthree=totalbadthree+1;
					}				
				badnesscounter=badnesscounter+1;	
			}
			while (((kone<0.1) or (kone>5000)) and (badnesscounter<1000));
			if (badnesscounter>=2)
			{
				totalbad=totalbad+1;	
				if (badnesscounter>=1000)
					totalbad=totalbad+1000000;							
			}
			outputone<<kone<<"	";
			outputtwo<<ktwo<<"	";
			outputthree<<knegone<<"	";
			outputfour<<knegtwo<<"	";
		}
		outputone<<endl;
		outputtwo<<endl;
		outputthree<<endl;
		outputfour<<endl;
	}
	outputfive<<totalbad<<"	"<<totalbadtwo<<"	"<<totalbadthree;
	outputone.close();
	outputtwo.close();
	outputthree.close();
	outputfour.close();
	outputfive.close();
}
