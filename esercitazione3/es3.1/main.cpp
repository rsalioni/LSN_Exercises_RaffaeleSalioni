

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

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
 /*
  for(int i=0; i<1000; i++){
  	cout << rnd.Rannyu() << endl;
  }
*/
//inizio esercizio
	//double t=0.;
	double S0 = 100.;
	double T=1.;
	double k=100.;
	double r=0.1;
	double sigma=0.25;
	double S;
	double avec=0.;
	double avec2=0.;
	double avep=0.;
	double avep2=0.;
	
	double Ntot = 100000; //lanci totali
	double N = 1000; //lanci per blocco
	double L = Ntot/N; //blocchi
	double sumc = 0.;
	double sump = 0.;
	
	//prima parte
	ofstream fout1("c1.txt");
	ofstream fout2("p1.txt");
	
	for(int i=0; i<L; i++){ //ciclo sui blocchi
		sumc = 0.;
		sump = 0.;
		for(int j=0; j<N; j++){ //singolo blocco
			S = S0*exp((r-0.5*sigma*sigma)*T + sigma*rnd.Gauss(0.,T));
			//cout << exp(-r*T)*fmax(0., S-k) << endl;
			sumc += exp(-r*T)*fmax(0., S-k); //incremento call
			sump += exp(-r*T)*fmax(0., k-S); //incremento put
		}
		avec += sumc/N;
		avec2 += pow(sumc/N, 2);
		avep += sump/N;
		avep2 += pow(sump/N, 2);
		fout1 <<  (i+1) << ", " << avec/(i+1) << ", " << error(avec/(i+1), avec2/(i+1), i) << endl;
		fout2 <<  (i+1) << ", " << avep/(i+1) << ", " << error(avep/(i+1), avep2/(i+1), i) << endl;
  		cout <<  (i+1) << ", C " << avec/(i+1)  << " err " << error(avec/(i+1), avec2/(i+1), i) << ", P " << avep/(i+1)  << " err " << error(avep/(i+1), avep2/(i+1), i) <<endl;
	}

	fout1.close();
	fout2.close();
    
    
    //seconda parte
    
    avec=0.;
	avec2=0.;
	avep=0.;
	avep2=0.;
	int npassi = 100;
	double passo = T/npassi;
    ofstream fout3("c2.txt");
	ofstream fout4("p2.txt");
	
	for(int i=0; i<L; i++){ //ciclo sui blocchi
		sumc = 0.;
		sump = 0.;
		for(int j=0; j<N; j++){ //singolo blocco
			S=S0;
			for(int k=0; k<npassi; k++){
				S = S*exp((r-0.5*sigma*sigma)*passo + sigma*rnd.Gauss(0.,1.)*sqrt(passo));
			}
			sumc += exp(-r*T)*fmax(0., S-k);
			sump += exp(-r*T)*fmax(0., k-S);
		}
		avec += sumc/N;
		avec2 += pow(sumc/N, 2);
		avep += sump/N;
		avep2 += pow(sump/N, 2);
		fout3 <<  (i+1) << ", " << avec/(i+1) << ", " << error(avec/(i+1), avec2/(i+1), i) << endl;
		fout4 <<  (i+1) << ", " << avep/(i+1) << ", " << error(avep/(i+1), avep2/(i+1), i) << endl;
  		cout <<  (i+1) << ", C " << avec/(i+1)  << " err " << error(avec/(i+1), avec2/(i+1), i) << ", P " << avep/(i+1)  << " err " << error(avep/(i+1), avep2/(i+1), i) <<endl;
	}
    
    fout3.close();
	fout4.close();
    

   rnd.SaveSeed();
   return 0;
}


