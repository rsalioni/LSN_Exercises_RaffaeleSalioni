#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>
//#include "Funzioni.h"

using namespace std;
Random my_rand;
int N=32;
int NN=500;
vector<vector<double>> Dist;

double error(double av, double av2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((av2- av*av)/n);
	}

int pbc(int a){
	if(a<N){
	return a;
	}
	else
	 return a%N + 1;
}

void swap (vector<int> &v, int i, int j){
	int dim = v.size();
	if(i < dim and j < dim){
		int appo;
		appo = v[i];
		v[i]=v[j];
		v[j]=appo;
	}
}



vector < vector<int> > Generate_Pop(int NN){
	vector <vector<int>> result;
	vector<int> v;
	int a, b;
	for(int i=0; i<N; i++){
		v.push_back(i+1);
	}
	result.push_back(v);
	for(int i=1; i<NN; i++){
		a=my_rand.Rannyu(1, N);
		b=my_rand.Rannyu(1, N);
		swap(v, a, b);
		result.push_back(v);
	}
return result;
}
vector<double> Generate_Pos_circ(int N){
	vector<double> theta;
	for(int i=0; i<N; i++){
		theta.push_back(my_rand.Rannyu(0., 2*M_PI));
	}
return theta;
}

vector<vector<double>> Generate_Pos_square(int N){
	vector <double> x;
	vector <double> y;
	vector<vector<double>> result;
	for(int i=0; i<N; i++){
		x.push_back(my_rand.Rannyu());
		y.push_back(my_rand.Rannyu());
	}
	result.push_back(x);
	result.push_back(y);
return result;
}
	
	

vector<vector<double>> Generate_Dist_circ(vector<double> theta, int N){
	
	vector<double> v;
	vector<vector<double>> result;
	double Mij, xi, xj, yi, yj;
		
	for(int i=0; i<N; i++){
		v.clear();
		xi = cos(theta[i]);
		yi = sin(theta[i]);
		for(int j=0; j<N; j++){
			xj = cos(theta[j]);
			yj = sin(theta[j]);
			Mij = pow((xi - xj), 2) + pow((yi - yj), 2);
			v.push_back(Mij);
		}
		
		result.push_back(v);
	}
	
return result;
}	

vector<vector<double>> Generate_Dist_square(vector<vector<double>> pos, int N){
	
	vector<double> v;
	vector<vector<double>> result;
	double Mij;
		
	for(int i=0; i<N; i++){
		v.clear();
		for(int j=0; j<N; j++){
			Mij = pow((pos[0][i] - pos[0][j]), 2) + pow((pos[1][i] - pos[1][j]), 2);
			v.push_back(Mij);
		}
		
		result.push_back(v);
	}
	
return result;
}	

void check(vector<vector<int>> v){
	
	bool check = true;
	int appo;
	for(int i=0; i<NN; i++){
		if(v[i][0] != 1){
			cout << "il primo elemento è diverso da 1!" << endl;
		}
		if(check == true){
			for(int j=0; j<N; j++){
				appo = v[i][j];
				for(int k=j+1; k<N; k++){
					if(appo == v[i][k]){
						check = false;
					}
				}
			}
		}
	}
	if(check == true)
		cout << "check ok" << endl;
	else
		cout << "check not ok" << endl;
} 
	

double L2 (vector<int> v){
	int dim = v.size();
	double sum=0.;
	int a, b; //indici posizione città
	for(int i=0; i<dim-1; i++){
		a = v[i] - 1;
		b = v[i+1] -1;
		sum += Dist[a][b];
	}
	a = v[dim-1] - 1;
	
	sum += Dist[a][0];
	
return sum;
}

void riordina(vector<vector<int>> &v){
	
	for(int i=0; i<NN; i++){
		for(int j=i; j<NN; j++){
			if(L2(v[j]) < L2(v[i])){
				v[i].swap(v[j]);
			}
		}
	}
}	


int selection(){
	int j;
	double r = my_rand.Rannyu(0., 1.);
	j = NN*pow(r, 2);
return j;
}

//mutazioni

void pair_permut(vector<int> &v){
	int N = v.size();
	int a=my_rand.Rannyu(1, N);
	int b;
	do{
	b=my_rand.Rannyu(1, N);
	}
	while(b == a);
	swap(v, a, b);
}


void m_permut(vector<int> &v){
	int N = v.size();
	int m = my_rand.Rannyu(2, N/2);
	
	int a=my_rand.Rannyu(1, N);
	int b=my_rand.Rannyu(1, N);
	for(int i=0; i<m; i++){	
		swap(v, pbc(a+i), pbc(b+i));	
	}
}
	
void inversion(vector<int> &v){
	int N = v.size();
	int m = my_rand.Rannyu(1, N+1);
	//cout << "m " << m << endl;
	int a=my_rand.Rannyu(1, N);
	
	for(int i=0; i<m/2; i++){
		swap(v, pbc(a+i), pbc(a+m-1-i));
	}
}



vector<vector<int>> New_Pop(vector<vector<int>> v){
	
	int gen1;
	int gen2;
	int pos;
	int count;
	int k;
	bool check;
	vector <int> genitore1(N);
	vector <int> genitore2(N);
	vector <int> figlio1(N);
	vector <int> figlio2(N);
	vector<vector<int>> result;
	for(int l=0; l<NN/2; l++){
		gen1 = selection();
		do{
			gen2=selection();
		}while(gen1 == gen2);
		genitore1 = v[gen1];
		genitore2 = v[gen2];
		
		if(my_rand.Rannyu() < 0.55){//probabilità crossover
			pos = my_rand.Rannyu(0, N);
			for(int i=0; i<pos; i++){
				figlio1[i]=genitore1[i];
				figlio2[i]=genitore2[i];
			}
			count = pos;
			k=0;
			check = true;
			while(count < N){
				check = true;
				for(int i=0; i<pos; i++){
					if(genitore2[k] == figlio1[i]){
						check = false;
					}
				}
				if(check == true){
					figlio1[count] = genitore2[k];
					count++;	
				}
				k++;
			}
			count = pos;
			k=0;
			while(count < N){
				check = true;
				for(int i=0; i<pos; i++){
					if(genitore1[k] == figlio2[i]){
						check = false;
					}
				}
				if(check == true){
					figlio2[count] = genitore1[k];
					count++;	
				}
				k++;
			}
		}
		else{
			figlio1 = genitore1;
			figlio2 = genitore2;
		}
		
		
		if(my_rand.Rannyu() < 0.1)
			pair_permut(figlio1);
		
		
		if(my_rand.Rannyu() < 0.1)
			m_permut(figlio1);
		
		if(my_rand.Rannyu() < 0.1)
			inversion(figlio1);
		
				
		if(my_rand.Rannyu() < 0.1)
			pair_permut(figlio2);
		
		if(my_rand.Rannyu() < 0.1)
			m_permut(figlio2);
			
		if(my_rand.Rannyu() < 0.1)
			inversion(figlio2);
			
		
		result.push_back(figlio1);
		result.push_back(figlio2);
		
		}
return result;
}
		



//main

int main (int argc, char *argv[]){

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
	         my_rand.SetRandom(seed,p1,p2);
	      }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	//cout << my_rand.Rannyu() << endl;
	//cout << my_rand.Rannyu() << endl;
	
	
	ofstream fout("L2_square.txt");
	ofstream fout2("path_square.txt");
	ofstream fout3("ave_square.txt");
	vector<vector<int>> Pop = Generate_Pop(NN);
	//vector<double> theta = Generate_Pos_circ(N);
	check(Pop);
	vector <vector<double>> Pos = Generate_Pos_square(N);
	
	//Dist = Generate_Dist_circ(theta, N);
	Dist = Generate_Dist_square(Pos, N);
	double sum, sum2;
	
	riordina(Pop);
	
	
	for(int i=0; i<500; i++){ //300
		Pop = New_Pop(Pop);
		cout << i << " " << endl;
		riordina(Pop);
		//check(Pop);
		sum = 0.;
		sum2 = 0.;
		for(int j=0; j<10; j++){
			sum += L2(Pop[j]);
			sum2 += pow(L2(Pop[j]),2);
		}
		sum = sum/10.;
		sum2 = sum2/10.;
			
		fout << i << ", " << L2(Pop[0]) << endl;
		fout3 << i << ", " << sum << ", " << error(sum, sum2, 10) << endl;
		
	}

	riordina(Pop);
	
	int city;
	for(int i=0; i<32; i++){
		city = Pop[0][i];
		cout << city << " ";
		fout2 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;
	}
	
	
	 fout.close();
	 fout2.close();
	 fout3.close();
	 
	 my_rand.SaveSeed();

return 0;
}

	
	
	

