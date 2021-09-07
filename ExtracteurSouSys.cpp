//Ce programme permet une extraction d'un sous systeme à partir d'une simulation lancée avec treses 2 (correspondance des fichiers de sortie). Il faut préciser plus bas le nombre d'événement avec lequel on a lancé treses 2 pour l'ecriture nM/Eve[loc] en particulier.

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
extern void seedMT2();
extern void seedMT(unsigned long int);
extern double ranMT(void);
using namespace std;

int ComptageEtat(int etat, vector<int> liste, int quest=0);

vector<double> Eve_l;
vector<int> Ind_l;
vector<vector<double> > Eve;
vector<vector<int> > Ind;
int nA=0,nB=0,nC=0;

int main (){
	int i,j,k,l, temp_i, temp_j,x,y,largeur=100;
	int Rg=3,id_centre=5050 ;// id_centre est le repère pour le centre du sous-système 
	string ligne;
	double test;
	vector<int> Sousgroupe;
	vector<int>::iterator it;


	ifstream fromEve("data_Eve_tresesN_10000.dat",ios::in);
	fromEve.precision(5);

	ifstream fromInd("data_Ind_tresesN_10000.dat",ios::in);
	fromInd.precision(5);
	
	y=floor(id_centre/largeur);
	x=id_centre-y*largeur;

	while (getline(fromEve,ligne)){			//On remplit un vecteur avec le fichier Eve généré par treses
		Eve_l.clear();
		for(i=0;i<3;i++){
			fromEve >> test;
			Eve_l.push_back(test);
			//cout << test << " ";
		}
		//cout << endl;
		Eve.push_back(Eve_l);		
	}


	while (getline(fromInd,ligne)){			//On remplit un vecteur avec le fichier Ind généré par treses
		Ind_l.clear();
		for(i=0;i<largeur;i++){
			fromInd >> test;
			Ind_l.push_back(test);
			//cout << test << " ";
		}
		//cout << endl;
		Ind.push_back(Ind_l);		
	}

	Sousgroupe.clear();
	ofstream toIndloc("data_Indloc_tresesN_10000.dat",ios::out | ios::trunc);
	ofstream toPosloc("data_Posloc_tresesN_10000.dat",ios::out | ios::trunc);
	for (i=y-Rg;i<y+Rg+1;i++){			//Calcul du pos_loc pour pouvoir faire l'affichage matlab
		for (j=x-Rg;j<x+Rg+1;j++){
			if(i!=y || j!=x){
				//périodicité horizontale
				if(j>=largeur){
					temp_j=j-largeur;
				}
				else if (j<0){
					temp_j=j+largeur;
				}
				else {
					temp_j=j;
				}
				//periodicité verticale
				if (i>=largeur){
					temp_i=i-largeur;
				}
				else if (i<0){
					temp_i=i+largeur;
				}
				else {
					temp_i=i;
				}
				Sousgroupe.push_back(temp_i*largeur+temp_j);
				toPosloc << temp_i*largeur+temp_j << " ";
			}
			else { toPosloc << id_centre << " ";}
			
		}
		toPosloc << endl;
	}
		

	l=0;
	for(it=Sousgroupe.begin();it!=Sousgroupe.end();++it){	//Creation du fichier Indloc et initialisation du nMloc
		y=floor(*it/largeur);
		x=*it-y*largeur;
		l++;
		if(l==(pow(Rg*2+1,2)-1)/2){
			toIndloc << Ind[x][y] << " " << Ind[x+1][y] << " " ;
			l++;
		}
		else{toIndloc << Ind[x][y] << " ";
		}			
		if (l%(Rg*2+1)==0){
			toIndloc << endl;
		}	
		if(Ind[x][y]==0){
			nA++;
		}
		else if(Ind[x][y]==1){
			nB++;
		}
		else if(Ind[x][y]==2){
			nC++;
		}
	}
	
	/*for (i=0;i<largeur;i++){
		nA+=ComptageEtat(0,Ind[i]);
		nB+=ComptageEtat(1,Ind[i]);
		nC+=ComptageEtat(2,Ind[i]);
	}*/


	ofstream toEveloc("data_Eveloc_tresesN_10000.dat",ios::out | ios::trunc);
	ofstream tonMloc("data_nMloc_tresesN_10000.dat",ios::out | ios::trunc);
	for (i=0;i<1000000;i++){	//Attention ici à borner par rapport au nombre d'evenement demandé dans treses; Boucle de creation de Eveloc et nMloc
		for(it=Sousgroupe.begin();it!=Sousgroupe.end();++it){
			if(*it==int(Eve[i][1])){
				toEveloc << Eve[i][0] << "\t \t" << Eve[i][1] << "\t \t" << Eve[i][2] << endl;
				if(Eve[i][2]==0){
					nA++;
					nC--;}
				else if(Eve[i][2]==1){
					nB++;
					nA--;}
				else if(Eve[i][2]==2){
					nC++;
					nB--;}
				tonMloc << Eve[i][0] << "\t \t" << nA << "\t \t" << nB << "\t \t" << nC << endl;
			}
		}
	}

	cout << "Fin du programme" << endl;

}

int ComptageEtat(int etat, vector<int> liste, int quest)
{
    int i;
    int nombreEtat=0;
    if (quest==0){
        for (i=0; i<liste.size(); i++)
        {
            //printf("%d %d\n",liste[i],etat);
            if (liste[i]==etat)
            {
                //printf("OK\n");
                nombreEtat+=1;
            }
        }
    }
    else{printf("Erreur ComptageEtat");}
    return nombreEtat;
}
