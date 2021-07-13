#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>
#include "Geometria.h"
#include "Funzioni.h"

using namespace std;
double error(double av, double av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}

 
int main (int argc, char *argv[]){

	Posizione p(1.5, 0., 0.); 
	int M = 1000000; //passi totali 
	int nblocchi = 100; //numero di blocchi
	int n = M/nblocchi; //passi di un blocco
	p.setD(0.75);
	double sum = 0.;
	H1S * f = new H1S();
	double ave=0.;
	double ave2=0.;
	int k = 0;
	ofstream fout1("rmedia1s_gauss.txt");
	ofstream fout2("pos1s_gauss.txt");
	//ofstream fout3("r_eq1s_gauss.txt"); //equilibrazione
	
	for(int i=0; i<nblocchi; i++){ //ciclo sui blocchi
		sum = 0.;
		for(int j=0; j<n; j++){ //blocco singolo
			p.Metropolis_gauss(f);
			k++;
			if(k%750 == 0){
				fout2 << p.getX() << ", " << p.getY() << ", " << p.getZ() << endl;
			}
			//if(k<250)
				//fout3 << p.getR()<<endl; //equilibrazione
			sum += p.getR();
			//cout << "r " << p.getR() << endl;
		}
		ave += sum/n;
		ave2 += pow(sum/n, 2);
		fout1 <<  (i+1) << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
  		cout <<  (i+1) << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
	}
	
	fout1.close();
	//fout2.close();
	
	cout << "n " << p.getN()/(double)M << endl;
	
return 0;
}
