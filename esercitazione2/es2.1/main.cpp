#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Integral.h"
#include "Funzioni.h"
using namespace std;

double error(double av, double av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}
 
int main (int argc, char *argv[]){
/*
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

*/

	//prima parte
	Coseno * f = new Coseno(M_PI/2., M_PI/2., 0. );
	Integral I(0.,1. , f);
	int Nb=100; //numero blocchi
	int Ntot = 1000000; //lanci totali
	int N = Ntot/Nb; //lanci per blocco
	double ave=0.;
	double ave2=0.;
	double integral;
	ofstream fout("output2.1.txt");
	
	for(int i=0; i<Nb; i++){ //ciclo sui blocchi
		integral = I.Media(N);
		ave += integral;
		ave2 += integral*integral;
		cout << i+1 << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
		fout << i+1 << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
	}
	
	
	fout.close();
	
	//seconda parte
	ave=0.;
	ave2=0.;
	
	Myfun* f1 = new Myfun();
	Integral I1(0.,1. , f1);
	ofstream fout1("output2.12.txt");
	for(int i=0; i<Nb; i++){ //ciclo sui blocchi
		integral = I1.Mediaretta(N);
		ave += integral;
		ave2 += integral*integral;
		cout << i+1 << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
		fout1 << i+1 << ", " << ave/(i+1) << ", " << error(ave/(i+1), ave2/(i+1), i+1) << endl;
	}
	
		
	
	
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   //rnd.SaveSeed();
   return 0;
}


