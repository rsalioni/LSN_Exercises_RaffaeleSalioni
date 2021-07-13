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

	int N = 100;// num blocchi
	int M = 100000; //lanci totali
	int L = M/N; //lanci per blocco
	double sigma = 0.61;
	double mu = 0.81;
	
	PsiT * psi = new PsiT(mu, sigma);
	V * v = new V();
	H * h = new H(psi, v);
	
	Metropolis metro(0., 5.);  //Metropolis metro(0., 5.) 
  	
  //vector <double> ave;
  //vector <double> ave2;
	double ave = 0.;
	double ave2 = 0.;
	double sum = 0.;
  //double sum2 = 0.;
	ofstream fout1("8.1.txt");
	ofstream fout2("histo.txt");
  
	for(int i=0; i<N; i++){ //ciclo sui blocchi
		sum = 0.;
  		for(int j=0; j<L; j++){//singolo blocco
  			metro.Move_unif(psi);
  			sum += h->Eval(metro.getX());
  			fout2 << metro.getX() << endl;
  		}
  	ave += sum/L;
  	ave2 += pow(sum/L, 2);
  	fout1 <<  (i+1)*L << ", " << ave/(i+1)  << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  	cout <<  (i+1)*L << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  }
  cout << "n " << metro.getN()/(double)M << endl;
  
  
  fout1.close();	
  fout2.close();	
  
return 0;  
}


