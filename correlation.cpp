#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

using namespace std;

const int N =10000;
double Etat[100][100];
int longueur =100, largeur=100;
const int kmax=floor(sqrt(pow(longueur-1,2)+pow(largeur-1,2)))+1;
double correl[N];
double Nk[N];

int main(){
	int i,j,i1,j1,i2,j2,k;
	double r;
	ifstream fromInd("data_Ind_tresesN_10000.dat");
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fromInd >> Etat[i][j];			
		}
	}
	fromInd.close();
	for(i1=0;i1<longueur;i1++){
		for(j1=0;j1<largeur;j1++){
			for(i2=0;i2<longueur;i2++){
				for(j2=0;j2<largeur;j2++){
					r=sqrt(pow(i1-i2,2)+pow(j1-j2,2));
					k=(int)floor(r);//rang 
					if (k<kmax+1)
					{	
						Nk[k]+=1; // nb de pts dans la couronne
						correl[k]+=Etat[i1][j1]*Etat[i2][j2];
						/*if (Etat[i1][j1]==Etat[i2][j2]){
							correl[k]+=1; //corrélation
						}
						else{
							correl[k]+=0; //decorrélation
						}*/
					}
				}
			}			
		}
	}
	ofstream toCor("Correlations.dat",ios::out | ios::trunc);
	for(i=0;i<kmax;i++){
		correl[i]/=Nk[i];
		toCor << correl[i] << endl;
	}
	toCor.close();
	
	
	return 0;
}
