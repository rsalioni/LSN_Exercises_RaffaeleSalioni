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
	
	//parametri
	double sigma_start= 0.15;
	double sigma;
	double mu = 0.33;
	int num_s = 100;
	int num_m = 100;
	double delta_s = 0.01;
	double delta_m = 0.01;
	
	//histo
	//int nbin = 100;
	//vector <int> histo;
	
	//for(int i=0; i<nbin+1; i++){
		//histo.push_back(0);;
	//}
	
	//double xmin = -2.;
	//double xmax = 2.;
	//double ampiezza = xmax - xmin;
	double appo;
	//int spot;
	
	cout << "test1" << endl;
	
	PsiT * psi = new PsiT(mu, sigma_start);
	V * v = new V();
	H * h = new H(psi, v);
	
	Metropolis metro(0., 4.);
  	
  //vector <double> ave;
  //vector <double> ave2;
	double ave = 0.;
	double ave2 = 0.;
	double sum = 0.;
  //double sum2 = 0.;
	ofstream fout1("8.2.txt");
	ofstream fout2("histo.txt");
	
	int count = 1;
  //int count_histo = 0;
  	
	for(int m = 0; m < num_m; m++){ //ciclo su mu
		sigma = sigma_start;
		psi->SetMu(mu);
		h->SetMu(mu);
		for(int s = 0; s < num_s; s++){ //ciclo su sigma
			//metro.SetX(0.); //inizializzo parametri
			metro.SetN(0);
			psi->SetSig(sigma);
			h->SetSig(sigma);
			
			ave = 0.;
			ave2 = 0.;
			
				for(int i=0; i<N; i++){ //ciclo sui blocchi
					sum = 0.;
  					for(int j=0; j<L; j++){//singolo blocco
  						metro.Move_unif(psi);
  						appo = metro.GetX();
  						sum += h->Eval(appo);
  						
  						/*
  						if( appo < xmax && appo > xmin){
  							count_histo++;
  							spot = nbin*(appo - xmin)/(xmax-xmin); 
  							histo[spot]++;
  						}
  						*/
  					}
  					ave += sum/L;
  					ave2 += pow(sum/L, 2);
  					//fout1 <<  (i+1)*L << ", " << ave/(i+1)  << ", " << error(ave/(i+1), ave2/(i+1), i) << endl;
  					
  				}
  			fout1 <<  mu << ", " << sigma << ", " <<  ave/N << ", " << error(ave/N, ave2/N, N) << endl;
  			cout << count << " n " << metro.GetN()/(double)M << endl;
  			
			count++;
			sigma += delta_s; //incremento sigma
		}
		mu += delta_m; //incremento mu	
	}
	 
	//histo
	/*
	psi->SetMu(0.805);
	h->SetMu(0.805);
	metro.SetN(0);
	psi->SetSig(0.613);
	h->SetSig(0.613);
	for(int i=0; i<M; i++){
		metro.Move_unif(psi);
  	appo = metro.GetX();
  	if( appo < xmax && appo > xmin){
  		count_histo++;
  		spot = nbin*(appo - xmin)/(xmax-xmin); 
  		histo[spot]++;
  	}
  }
  	*/
  	
	
	
	
	/*
			
	for(int i=0; i<nbin; i++){
		cout << ((xmax-xmin)/nbin)*(i+0.5) +xmin << ", " << i << ", " << histo[i]/(double)count_histo << endl;
		fout2 << ((xmax-xmin)/nbin)*(i+0.5) +xmin << ", " << histo[i]/((double)count_histo) << endl;
	}	
	
	cout << "count histo " << count_histo << endl;
	*/

  fout1.close();	
  //fout2.close();
  
return 0;  
}


