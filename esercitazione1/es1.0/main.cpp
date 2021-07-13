#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
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

  //prima parte, media r
  cout << "prima parte" << endl;
  
  int N = 100;// num blocchi
  int M = 100000; //lanci totali
  int L = M/N; //lanci per blocco
  

  double ave = 0.;
  double ave2 = 0.;
  double sum = 0.;

  ofstream fout1("output1.01test.txt");
  
  for(int i=0; i<N; i++){ //ciclo sui blocchi
  	sum = 0.;
  	for(int j=0; j<L; j++){//singolo blocco
  		sum += rnd.Rannyu();
  	}
  	ave += sum/L;
  	ave2 += pow(sum/L, 2);
  	fout1 <<  (i+1)*L << ", " << ave/(i+1) - 0.5 << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  	cout <<  (i+1)*L << ", " << ave/(i+1) - 0.5 << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  }
  
  
  
  fout1.close();
 
 
  //seconda parte
  cout << "seconda parte" << endl;
  ave = 0.;
  ave2 = 0.;
  ofstream fout2("output1.02test.txt");
  for(int i=0; i<N; i++){
  	sum = 0.;
  	for(int j=0; j<L; j++){
  		sum += pow((rnd.Rannyu()-0.5), 2);
  	}
  	ave += sum/L;
  	ave2 += pow(sum/L, 2);
  	fout2 <<  (i+1)*L << ", " << ave/(i+1) - 1./12. << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  	cout <<  (i+1)*L << ", " << ave/(i+1) - 1./12. << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  }
  
  
  
  fout2.close();
  
  
  //terza parte 
  cout << "terza parte" << endl;
  
   ofstream fout3("output1.03.txt");
  M = 100;//numero intervalli
  N = 10000;//numero lanci
  double delta = 1./M;
  vector <int> n;
  int k;
  double chi = 0.;
 
  for(int i=0; i<M; i++){    //da sistemare
   		n.push_back(0);
  }
  for(int j=0; j<100; j++){
  	
  	chi = 0.;
  	for(int i=0; i<M; i++){    //da sistemare
   		n[i]=0;
  	}
   	for(int i=0; i < 10000; i++){
		k = rnd.Rannyu()/delta;
		n[k] = n[k] + 1;
		}
	for(int i=0; i<M; i++){
		chi += pow(n[i]-100,2)/(100.);
	}
	cout << chi << endl;
	fout3 << chi << endl;
  }
  
  fout3.close();
  

   rnd.SaveSeed();
   return 0;
}

