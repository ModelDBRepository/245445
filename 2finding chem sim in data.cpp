#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <iomanip>
using namespace std;

const int noofodour = 16;
const int knownodourant = 26;

int main()
{
	int i;
	int j;
	int k;
	int m;
	double chemsim[noofodour][noofodour-1];	
	double response[knownodourant][noofodour];
				
	ifstream inputfile;
	inputfile.open("ORNdata.txt");
	for (i=0;i<knownodourant;i++)
		for (j=0;j<noofodour;j++)
		inputfile>>response[i][j];
	inputfile.close();	

	for (j=0;j<noofodour;j++)
		for (k=0;k<noofodour-1;k++)	
			chemsim[j][k]=0;
			
	for (j=1;j<noofodour;j++)
		for (m=0;m<j;m++)
		{
			for (i=0;i<knownodourant;i++)
				chemsim[j][m]=chemsim[j][m]+(response[i][j]-response[i][m])*(response[i][j]-response[i][m]);
			chemsim[j][m]=sqrt(chemsim[j][m]/knownodourant);
		}
		
	ofstream output;
	output.open("ORNdatachemsim.txt");
	for (j=1;j<noofodour;j++)
		{
		for (m=0;m<j;m++)
			output<<chemsim[j][m]<<"	";
		output<<endl;
		}
		output.close(); 
return 1;	
}
