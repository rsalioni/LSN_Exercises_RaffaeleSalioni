
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>
#include "Geometria.h"

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

   
	//da qui inizio a scrivere io
	
	Posizione p(0., 0., 0.); 
	Posizione origine(0., 0., 0.);
	int Ntot = 10000; 
	int Nb = 100;	//rw per blocco
	int L = Ntot/Nb; //numero blocchi
	int N = 100; //step per rw
	cout << "n blocchi " << L << endl;
	vector <double> sum (N);
	vector <double> ave (N);
	vector <double> ave2 (N);
	ofstream fout1("output2.2.txt");
	
	for(int i=0; i<N; i++){
			ave[i]=0.;
			ave2[i]=0.;
			
		}
	
	double r;
	for(int j=0; j<L; j++){   //ciclo sui blocchi
		for(int i=0; i<N; i++){
			sum[i]=0.;
		}
		
		
		for(int i=0; i<Nb; i++){    //ciclo sui randomwalk di un blocco
			p.setX(0.);
			p.setY(0.);
			p.setZ(0.);
			for(int k=1; k<N; k++){   //singolo randomwalk
				r=rnd.Rannyu();
				if(r < 1./6.)
					p.setX(p.getX() + 1.);
				else if(1./6. < r and r < 2./6.)
					p.setX(p.getX() - 1.);
				else if(2./6. < r and r < 3./6.)
					p.setY(p.getY() + 1.);
				else if(3./6. < r and r < 4./6.)
					p.setY(p.getY() - 1.);
				else if(4./6. < r and r < 5./6.)
					p.setZ(p.getZ() + 1.);
				else if( r > 5./6.)
					p.setZ(p.getZ() - 1.);
				//cout <<k << " " << p.getX() << ", " <<  p.getY() << ", " <<  p.getZ() << endl;
				sum[k] += pow(p.getR(), 2);
			}
			
		}
		//cout << sum[30]/N << endl;  sum media di un blocco
		for(int i=0; i< N; i++){
			ave[i] += sum[i]/Nb;
			ave2[i] += pow(sum[i]/Nb, 2);
		}
	}
	/*
	for(int i=0; i<N; i++){
		cout << ave[i] << " " << ave2[i] << endl;
	}
	*/
	cout << ave[0] << " " << ave2[0] << " " << error(0., 0., N) << endl;	
	fout1 << 0 << ", " << 0 << ", " << 0 << endl;
	for(int i=1; i<N; i++){
		fout1 << i << ", " << sqrt(ave[i]/L) << ", " << error(ave[i]/L, ave2[i]/L, L)/(2*sqrt(ave[i]/L)) << endl;
		//fout1 << i << ", " << sqrt(ave[i]/N) << endl;
		cout << i << ", " << sqrt(ave[i]/N) << ", " << error(ave[i]/N, ave2[i]/N, N)/(2*sqrt(ave[i]/N)) <<endl;;	
	}
		
	fout1.close();	
	
	//parte2			
	cout << "inizio parte2" << endl;		
	
	double theta;
	double phi;			
	ofstream fout2("output2.21.txt");
	
	for(int i=0; i<N; i++){
			ave[i]=0.;
			ave2[i]=0.;
			
		}
	
	
	for(int j=0; j<L; j++){   //ciclo sui blocchi
		for(int i=0; i<N; i++){
			sum[i]=0.;
		}
		
		
		for(int i=0; i<Nb; i++){    //ciclo sui randomwalk di un blocco
			p.setX(0.);
			p.setY(0.);
			p.setZ(0.);
			for(int k=1; k<N; k++){   //singolo randomwalk
				theta=rnd.Rannyu(0., M_PI);
				phi=rnd.Rannyu(0.,2*M_PI);
				p.setX(p.getX() +sin(theta)*cos(phi));
				p.setY(p.getY() + sin(theta)*sin(phi));
				p.setZ(p.getZ() + cos(theta));
				//cout <<k << " " << p.getX() << ", " <<  p.getY() << ", " <<  p.getZ() << endl;
				sum[k] += pow(p.getR(), 2);
			}
			
		}
		//cout << sum[30]/N << endl;  sum media di un blocco
		for(int i=0; i< N; i++){
			ave[i] += sum[i]/Nb;
			ave2[i] += pow(sum[i]/Nb, 2);
		}
	}
	/*
	for(int i=0; i<N; i++){
		cout << ave[i] << " " << ave2[i] << endl;
	}
	*/
	cout << ave[0] << " " << ave2[0] << " " << error(0., 0., N) << endl;	
	fout1 << 0 << ", " << 0 << ", " << 0 << endl;
	for(int i=1; i<N; i++){
		fout2 << i << ", " << sqrt(ave[i]/L) << ", " << error(ave[i]/L, ave2[i]/L, L)/(2*sqrt(ave[i]/L)) << endl;
		//fout1 << i << ", " << sqrt(ave[i]/N) << endl;
		cout << i << ", " << sqrt(ave[i]/N) << ", " << error(ave[i]/N, ave2[i]/N, N)/(2*sqrt(ave[i]/N)) <<endl;;	
	}
		
	fout1.close();			
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	rnd.SaveSeed();
	return 0;
}

