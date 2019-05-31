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
const double varconst= 0.005;
const double varsd = 0.028;
const double meanone = 0.27;
const double meanvarone = 0.27;
//const double meansegre = 0.308;
const double meantwo = 0.53;
const double meanvartwo = 0.35;
const double smallres = 0.01;
const double softminthres = -0.15;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double min[noofodourant];
	double softmin[noofodourant];
	double max[noofodourant];
	double response[noofodourant][noofodour];
	double s[noofodourant][noofodourant-1];
//	double obtainings(int indexofnewodourant, int indexofpreviousodourant);	
	double getnewresponse(double responseofpreviousodourant, double s, double minprevious, double maxprevious);
//	double getnoise (double averages);
	double temprestot;
	double temprestotsq;
	double tempone;
	double temptwo;
	double var;
	double varmean;
	double normcounter;
	double noise[noofodour];
	double propofnoise[noofodourant];
			
	ifstream inputfile;
	inputfile.open("ORNdata.txt");
	for (i=0;i<knownodourant;i++)
		for (j=0;j<noofodour;j++)
		inputfile>>response[i][j];
	inputfile.close();	
	
	ofstream outzero;	
//	outzero.open("knownsets.txt");
//	for (i=0;i<knownodourant;i++)
//	{
//		for (j=0;j<noofodour;j++)
//			outzero<<response[i][j]<<" ";
//		outzero<<endl;
//	}
//	outzero.close();	
	
	for (i=0;i<knownodourant;i++)
	{	
		min[i]=1000;
		max[i]=-1000;
		for (j=0;j<noofodour;j++)
		{
			if (response[i][j]<min[i])							  //finding min and max so that it can be scaled to [0,1] and be used in eq. 4.9
				min[i] = response[i][j];
			if (response[i][j]>max[i])
				max[i] = response[i][j];	 		
		}
	}

		for (k=1;k<noofodourant;k++)						//to obtain s
		{
			for (m=0;m<k;m++)
			{
				tempone = RAN;
				if (RAN<0.16)
				{
					do 
					{
						s[k][m]= 0.88+0.075*Norm;	
					}
					while ((s[k][m]>1) or (s[k][m]<0));
				}
				else if (RAN>0.76)
				{
					do 
					{
						s[k][m]= 0.16+0.075*Norm;	
					}
					while ((s[k][m]>1) or (s[k][m]<0));
				}
				else				
				{
					do 
					{
						s[k][m]= 0.64+0.075*Norm;	
					}
					while ((s[k][m]>1) or (s[k][m]<0));
				}			
			}
			do 
			{
				propofnoise[k]= 0.9+Norm;	                        //modulate the strength of output correlation
			}
			while ((propofnoise[k]>2) or (propofnoise[k]<0));
		}
		
//		ofstream outone;	
//		outone.open("sandnoise.txt");
//		for (k=1;k<noofodourant;k++)	
//		{
//			for (m=0;m<k;m++)
//				outone<<s[k][m]<<"	";
//			outone<<endl;
//		}
//		outone<<endl;
//		for (k=1;k<noofodourant;k++)	
//			outone<<propofnoise[k]<<"	";
//		outone<<endl;
//		outone.close();		

//		s[1][0]=0.25;
//		s[2][0]=1;
//		s[2][1]=1;
//		s[3][0]=0.5;
//		s[3][1]=0.5;
//		s[3][2]=0.5;
	
		temprestot = 0;
		temprestotsq = 0;


		for (k=knownodourant;k<noofodourant;k++)
		{
			temprestot = 0;
			temprestotsq = 0;
			normcounter =0;
			for (j=0;j<noofodour;j++) 
			{
				noise[j]=0; 
				response[k][j]=0;     				 
			}
			for (m=0;m<k;m++)
			{
				normcounter = normcounter + 2*(abs(s[k][m]-0.5));
				for (j=0;j<noofodour;j++) 
				{
					response[k][j]=response[k][j]+getnewresponse(response[m][j], s[k][m],min[m],max[m]);   //1st or 3rd terms of equation (4.9)
					noise[j]=noise[j]+2*RAN00*(max[m]-min[m])*(0.5-(abs(s[k][m]-0.5)));  //2nd term of equation (4.9)
				} 
			}

			for (j=0;j<noofodour;j++) 
			{
				noise[j]=noise[j]*propofnoise[k];
				response[k][j]=response[k][j]+noise[j];
				response[k][j]=response[k][j]/(normcounter+(k-normcounter)*propofnoise[k]);  //normalization: missing in 1st or 3rd terms equation (4.9)
				temprestot = temprestot + response[k][j];			//finding mean and varaince
				temprestotsq = temprestotsq + response[k][j]*response[k][j];				
			}
			temprestot = temprestot/(noofodour);			//next 10 lines or so: same operation as in the 1st odourant
			temprestotsq = temprestotsq/(noofodour);
			var = temprestotsq - temprestot*temprestot;
			min[k] = 1000;
			softmin[k] = 1000;
			max[k] = -1000;			
			varmean=abs(varsd*Norm);
			tempone = RAN-0.5;
			temptwo = RAN;
			for (j=0;j<noofodour;j++)
			{
				if (temptwo<0.25)
					response[k][j]= (response[k][j]-temprestot)*sqrt(varmean/var)+meantwo+meanvartwo*tempone;  //scaling the variance and mean of the response
				else 
					response[k][j]= (response[k][j]-temprestot)*sqrt(varmean/var)+meanone+meanvarone*tempone;				
				if (response[k][j]<min[k])	   //finding min and max so that it can be scaled to [0,1] and be used in eq. 4.9
				{
					if (response[k][j]>softminthres)	
						softmin[k]=response[k][j];					  
					min[k] = response[k][j];
				}
				if (response[k][j]>max[k])
					max[k] = response[k][j];			
			}
		}

		for (k=knownodourant;k<noofodourant;k++)
		{
			if (softmin[k]<0)
			{
				tempone=smallres*RAN;
				for (j=0;j<noofodour;j++)
					response[k][j]=response[k][j]-softmin[k]+tempone;      //increase overall response strength a little bit
			}
		}
	
		for (k=knownodourant;k<noofodourant;k++)
			for (j=0;j<noofodour;j++)		
			{
				if (response[k][j]<0)
					response[k][j]= smallres*RAN;                    //make sure that all responses are positive
			}								

		ofstream outtwo;
		outtwo.open("responseinter.txt");
		for (j=0;j<noofodour;j++)
		{
			for (k=0;k<noofodourant;k++)
				outtwo<<response[k][j]<<"	";
			outtwo<<endl;
		}	
		outtwo.close(); 
return 1;	
}

//double obtainings(int k,int m)
//{
//	double sspec[4][3];
//	sspec[1][0]=0.25;
//	sspec[2][0]=1;
//	sspec[2][1]=1;
//	sspec[3][0]=0.5;
//	sspec[3][1]=0.5;
//	sspec[3][2]=0.5;
//	return (sspec[k][m]);
//}

double getnewresponse (double pres, double s, double min, double max)    //1st or 3rd terms of equation (4.9)
{
	double scaledpres;
	//scaledpres = (pres-min)/(max-min);		//scaled to interval between 0 and 1 
	if (s>=0.5)
		return (2*pres*(s-0.5));		//strong response if previous OR responding strongly and the two ORs are similar (factor 2 missing in thesis' equation)
	else
		return (2*(1-pres)*(0.5-s));	//strong response if previous OR responding weakly and the two ORs are dissimilar
}

//double getnoise (double s)     //2nd term of equation (4.9)
//{
//	return (2*RAN*(0.5-(abs(s-0.5))));        //a lot of noise if neither similar nor dissimilar with previous (i.e. uncorrelationed)
//}




