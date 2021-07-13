#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>
#include <cstdlib>
#include "mpi.h"
//#include "Funzioni.h"

using namespace std;
Random my_rand;
Random my_rand0;
int N=32;
int NN=500;
int Nmigr = 25;
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
		a=my_rand0.Rannyu(1, N);
		b=my_rand0.Rannyu(1, N);
		swap(v, a, b);
		result.push_back(v);
	}
return result;
}
vector<double> Generate_Pos_circ(int N){
	vector<double> theta;
	for(int i=0; i<N; i++){
		theta.push_back(my_rand0.Rannyu(0., 2*M_PI));
	}
return theta;
}

vector<vector<double>> Generate_Pos_square(int N){
	vector <double> x;
	vector <double> y;
	vector<vector<double>> result;
	for(int i=0; i<N; i++){
		x.push_back(my_rand0.Rannyu());
		y.push_back(my_rand0.Rannyu());
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

	int size, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat1, stat2;
	
	
	int seed0[4];
	int p01, p02;
	ifstream Primes0("Primes");
	if (Primes0.is_open()){
		Primes0 >> p01 >> p02;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes0.close();

	ifstream input0("seed.in");
	string property0;
	if (input0.is_open()){
	   while ( !input0.eof() ){
	      input0 >> property0;
	      if( property0 == "RANDOMSEED" ){
	         input0 >> seed0[0] >> seed0[1] >> seed0[2] >> seed0[3];
	         my_rand0.SetRandom(seed0,p01,p02);
	      }
	   }
	   input0.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	for(int i=0; i<4; i++){
		if(rank == i){
			ifstream input("seed" + to_string(i) + ".in");
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
			
		}
	}
	
	
	//cout << my_rand.Rannyu() << endl;
	//cout << my_rand.Rannyu() << endl;
	
	
	ofstream fout_rk0("L2_squarerk0.txt");
	ofstream fout2_rk0("path_squarerk0.txt");
	ofstream fout3_rk0("ave_squarerk0.txt");
	
	ofstream fout_rk1("L2_squarerk1.txt");
	ofstream fout2_rk1("path_squarerk1.txt");
	ofstream fout3_rk1("ave_squarerk1.txt");
	
	ofstream fout_rk2("L2_squarerk2.txt");
	ofstream fout2_rk2("path_squarerk2.txt");
	ofstream fout3_rk2("ave_squarerk2.txt");
	
	ofstream fout_rk3("L2_squarerk3.txt");
	ofstream fout2_rk3("path_squarerk3.txt");
	ofstream fout3_rk3("ave_squarerk3.txt");	
	
	
	
	vector<vector<int>> Pop = Generate_Pop(NN);
	//vector<double> theta = Generate_Pos_circ(N);
	check(Pop);
	vector <vector<double>> Pos = Generate_Pos_square(N);
	
	//Dist = Generate_Dist_circ(theta, N);
	Dist = Generate_Dist_square(Pos, N);
	double sum, sum2;
	
	int* imesg = new int[N]; 
	int* imesg2 = new int[N];
	int itag=1; int itag2=2;
	vector <int> swap (4);
	vector <bool> preso(4);
	//int jj;
	
	riordina(Pop);
	for(int i=0; i<300; i++){
		
		
		Pop = New_Pop(Pop);
		for(int j=0; j<4; j++){
			if(rank == j){
				cout <<"rk " << j << " " <<  i << " " << endl;
			}
		}
		for(int j=0; j<4; j++){
			preso[j] = false;
		}
		if(i%Nmigr == 0 && i != 0){
			swap[0] = my_rand0.Rannyu(0., 4.);
			preso[swap[0]] = true;
			do{
				swap[1] = my_rand0.Rannyu(0., 4.);
			}while(swap[0] == swap[1]);
			preso[swap[1]] = true;
			
			for(int j=0; j<4; j++){
				if(preso[j] == false){
					swap[2] = j;
				}
			}
			swap[3] = 6 - swap[0] - swap[1] - swap[2];
			cout << "swap ";
			for(int j=0; j<4; j++){
			 	cout << swap[j] << " ";
			 }
			 cout << endl; 
			 
			
			//scambio i migliori delle prime due popolazioni
			for(int j=0; j<N; j++){	
			 	imesg[j] = Pop[0][j];
			 	imesg2[j] = Pop[0][j];
			}
			 
			if(rank==swap[1]){
				MPI_Send(&imesg[0],N,MPI_INTEGER,swap[0],itag,MPI_COMM_WORLD);
				MPI_Recv(&imesg2[0],N,MPI_INTEGER,swap[0],itag2, MPI_COMM_WORLD,&stat2);
				//cout<<"messaggio1 = "<<imesg2[0]<<endl;
			}
			else if(rank==swap[0]){
				MPI_Send(&imesg2[0],N,MPI_INTEGER,swap[1],itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesg[0],N,MPI_INTEGER,swap[1],itag, MPI_COMM_WORLD,&stat1);
				//cout<<"messaggio = "<<imesg[0]<<endl;
			}
			
			for(int j=0; j<N; j++){
				if(rank==swap[1]){
					Pop[0][j] = imesg2[j];
				}
				else if(rank==swap[0]){
					Pop[0][j] = imesg[j];
				}
			}
			
			
			//scambio i migliori delle seconde due popolazioni
			for(int j=0; j<N; j++){	
			 	imesg[j] = Pop[0][j];
			 	imesg2[j] = Pop[0][j];
			}
			 
			if(rank==swap[3]){
				MPI_Send(&imesg[0],N,MPI_INTEGER,swap[2],itag,MPI_COMM_WORLD);
				MPI_Recv(&imesg2[0],N,MPI_INTEGER,swap[2],itag2, MPI_COMM_WORLD,&stat2);
				//cout<<"messaggio1 = "<<imesg2[0]<<endl;
			}
			else if(rank==swap[2]){
				MPI_Send(&imesg2[0],N,MPI_INTEGER,swap[3],itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesg[0],N,MPI_INTEGER,swap[3],itag, MPI_COMM_WORLD,&stat1);
				//cout<<"messaggio = "<<imesg[0]<<endl;
			}
			
			for(int j=0; j<N; j++){
				if(rank==swap[3]){
					Pop[0][j] = imesg2[j];
				}
				else if(rank==swap[2]){
					Pop[0][j] = imesg[j];
				}
			}
				 
			 		
		}
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
		
		if(rank == 0){
			fout_rk0 << i << ", " << L2(Pop[0]) << endl;
			fout3_rk0 << i << ", " << sum << ", " << error(sum, sum2, 10) << endl;
		}
		else if(rank == 1){
			fout_rk1 << i << ", " << L2(Pop[0]) << endl;
			fout3_rk1 << i << ", " << sum << ", " << error(sum, sum2, 10) << endl;
		}
		else if(rank == 2){
			fout_rk2 << i << ", " << L2(Pop[0]) << endl;
			fout3_rk2 << i << ", " << sum << ", " << error(sum, sum2, 10) << endl;
		}
		else if(rank == 3){
			fout_rk3 << i << ", " << L2(Pop[0]) << endl;
			fout3_rk3 << i << ", " << sum << ", " << error(sum, sum2, 10) << endl;
		}
		
		
	}

	riordina(Pop);
	int city;
	
	
	for(int i=0; i<32; i++){
		city = Pop[0][i];
		cout << city << " ";
		if(rank == 0)
			fout2_rk0 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;
		else if(rank == 1)
			fout2_rk1 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;
		else if(rank == 2)
			fout2_rk2 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;
		else if(rank == 3)
			fout2_rk3 << Pos[0][city-1] << ", " << Pos[1][city-1] << endl;						
			
	}
	
	 fout_rk0.close();
	 fout2_rk0.close();
	 fout3_rk0.close();
	 
	 fout_rk1.close();
	 fout2_rk1.close();
	 fout3_rk1.close();
	 
	 //my_rand.SaveSeed();
	 
	 MPI_Finalize();

return 0;
}

