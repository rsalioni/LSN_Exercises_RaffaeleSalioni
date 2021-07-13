#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;
Random my_rand;
int N=32;
//int NN=500;
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

//mutazioni

vector<int> pair_permut(vector<int> v){
	int N = v.size();
	int a=my_rand.Rannyu(1, N);
	int b;
	do{
	b=my_rand.Rannyu(1, N);
	}
	while(b == a);
	swap(v, a, b);
return v;
}


vector<int> m_permut(vector<int> v){
	int N = v.size();
	int m = my_rand.Rannyu(2, N/2);
	
	int a=my_rand.Rannyu(1, N);
	int b=my_rand.Rannyu(1, N);
	for(int i=0; i<m; i++){	
		swap(v, pbc(a+i), pbc(b+i));	
	}
return v;
}
	
vector<int> inversion(vector<int> v){
	int N = v.size();
	int m = my_rand.Rannyu(1, N+1);
	//cout << "m " << m << endl;
	int a=my_rand.Rannyu(1, N);
	
	for(int i=0; i<m/2; i++){
		swap(v, pbc(a+i), pbc(a+m-1-i));
	}
return v;
}

vector<int> shift(vector<int> v){
	int N = v.size();
	int m = my_rand.Rannyu(1, N-1);
	int start=my_rand.Rannyu(1, N);
	int n = my_rand.Rannyu(1, N-1);
	vector<int> result(N);
	
	if(start + m < N){
		for(int i=0; i<start; i++){
			result[i] = v[i];
		}
		for(int i=start; i<N; i++){
			if(i+n<N){
				result[i+n] = v[i];
			}
			else
				result[(i+n)%N + start] = v[i];
		}
	}
	else{
		result[0] = v[0];
		for(int i=1; i<N; i++){
			result[pbc(i+n)] = v[i];
		}
	}
return result;
}


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

	
	vector<int> path;
	vector<int> new_path;
	double temp = 10.0;
	double delta_temp = 0.001; //0.001
	double delta_L;
	double r;
	int city;
	int kk = 0;
	
	ofstream fout("L2_square.txt");
	ofstream fout2("path_square.txt");
	
	//inizializzo path
	for(int i=0; i<N; i++){
		path.push_back(i+1);
	}
	for(int i=1; i<N; i++){
		swap(path, my_rand.Rannyu(1, N), my_rand.Rannyu(1, N));
	}
	
	
	vector <vector<double>> Pos = Generate_Pos_square(N);
	Dist = Generate_Dist_square(Pos, N);
	
	
	while( temp > 0.){
		for(int i=0; i<100; i++){ //temperatura costante
			//pair permut
			new_path = pair_permut(path);
			delta_L = L2(new_path) - L2(path);
			if(delta_L > 0.){
				r = my_rand.Rannyu();
				if(r < exp(-delta_L/temp)){
					path = new_path;
				}
			}
			else{
				path = new_path;
			}
			
			//m permut
			new_path = m_permut(path);
			delta_L = L2(new_path) - L2(path);
			if(delta_L > 0.){
				r = my_rand.Rannyu();
				if(r < exp(-delta_L/temp)){
					path = new_path;
				}
			}
			else{
				path = new_path;
			}
			
			//inversion
			new_path = inversion(path);
			delta_L = L2(new_path) - L2(path);
			if(delta_L > 0.){
				r = my_rand.Rannyu();
				if(r < exp(-delta_L/temp)){
					path = new_path;
				}
			}
			else{
				path = new_path;
			}
		}
		
		fout << kk << ", " << L2(path) << endl;
		cout << temp << endl;
		temp = temp - delta_temp;
		kk++;
	}
	
	
	for(int i=0; i<N; i++){
		city = path[i];
		cout << city << " ";
		fout2 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;
	}
	
	
	fout.close();
	fout2.close();
	
return 0;
}


