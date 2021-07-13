#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include <cstdlib>


using namespace std;
 
int main (int argc, char *argv[]){

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
	/*
   ofstream fout("output1.10.txt");
   ofstream fout1("output1.11.txt");
   for(int i=0; i<2000; i++){
      //cout << rnd.Exp(1.) << endl;
      fout << rnd.Exp(1.) << endl;
      fout1 << rnd.Lorentz(1., 0.) << endl;
   }
   fout.close();
   fout1.close();
	*/
	vector<double> sum = {0., 0., 0.};
	vector<int> N = {1, 2, 10, 100};
	int n = 10E4;
	vector<string> filename (3);
	vector<string> name = {"std", "exp", "lor"};
	ofstream foutstd;
	ofstream foutexp;
	ofstream foutlor;
	
	for(int j=0; j<4; j++){   //Nsize 
		for(int i=0; i<3; i++){
			filename[i] = name[i] + to_string( N[j]) + ".txt";
		}
		foutstd.open(filename[0]);
		foutexp.open(filename[1]);
		foutlor.open(filename[2]);
		for(int k=0; k<n; k++){
			sum = {0., 0., 0.} ;
			for(int i=0; i<N[j]; i++){
				sum[0] += rnd.Rannyu();
				sum[1] += rnd.Exp(1.);
				sum[2] += rnd.Lorentz(1., 0.);
			}
			foutstd << sum[0]/N[j] << endl;
			foutexp << sum[1]/N[j] << endl;
			foutlor << sum[2]/N[j] << endl;
		}
		foutstd.close(); 
		foutexp.close(); 
		foutlor.close();
		for(int i=0; i<3; i++)
			cout << "caricato " << filename[i] << endl; 
	}
	

	
   rnd.SaveSeed();
   return 0;
}


