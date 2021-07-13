#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>
using namespace std;

double error(double av, double av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}
 
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
   
   
	int Nhit = 0;
	int Ntot = 50000; //num lanci per blocco
	double d = 1.;
	double L = 0.4;
	double y1, y2, theta;
	 	
	
	 //qua inizia l'esercizio
	//vector<double> ave;
	//vector<double> ave2;
	double ave =0.;
	double ave2 = 0.;
	
	//double sum=0.;
	//double sum2=0.;
	int N = 100; //num blocchi
	ofstream fout("pi.txt");
	
	for(int j=0; j<N; j++){ //ciclo sui blocchi
		Nhit=0;
		for(int i=0; i<Ntot; i++){ //lanci di un blocco
			y1=rnd.Rannyu(0., d);
			theta=rnd.Angle();
			//theta=rnd.Rannyu(-M_PI, M_PI);
			y2=L*sin(theta)+y1;
			if(!(y2>0. and y2<d)){
				Nhit++;
			}
		}
		//ave.push_back((2*L*Ntot)/(Nhit*d)); //pi
		//ave2.push_back(ave[j]*ave[j]);
		ave += (2*L*Ntot)/(Nhit*d);
		ave2 += pow((2*L*Ntot)/(Nhit*d), 2);
		cout << "num blocchi: " << j+1 << " pi: " << ave/(j+1) << " errore: " << error(ave/(j+1), ave2/(j+1), j) << endl;
		fout << j+1 << ", " << ave/(j+1) - M_PI << ", " << error(ave/(j+1), ave2/(j+1), j) << endl; 
		
	}
	//cout << "ave e ave2 filled" << endl;	
	/*
	for(int j=10; j<N+1; j=j+10){ //ne esplicito uno ogni 10
		for(int i=j-10; i<j; i++){
			sum+=ave[i];
			sum2+=ave2[i];	
		}
		cout << "num blocchi: " << j << " pi: " << sum/j << " errore: " << error(sum/j, sum2/j, j) << endl;
		fout << j << ", " << sum/j - M_PI << ", " << error(sum/j, sum2/j, j) << endl; 
	}
	
	*/
	
	
	rnd.SaveSeed();
	return 0;
}

