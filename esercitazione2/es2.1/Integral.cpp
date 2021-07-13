#include<iostream>
#include"Integral.h"
#include <fstream>
#include"random.h"
//#include"Funzioni.h"

using namespace std;

/*
Integral::Integral( double a, double b, const FunzioneBase * f ){
	if ( b > a ){
		m_sign = 1;
		m_a = a;
		m_b = b;
	} else {
		m_sign = -1;	
		m_b = a;
		m_a = b;
	}
	m_f = f;
	m_sum = 0;
	m_integral = 0;
	m_h = 0;
}
*/

Integral::Integral( double a, double b, const FunzioneBase * f){
	if ( b > a ){
		m_sign = 1;
		m_a = a;
		m_b = b;
	} else {
		m_sign = -1;	
		m_b = a;
		m_a = b;
	}
	m_f = f;
	m_sum = 0;
	m_integral = 0;
	m_h = 0;
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
	         m_rnd.SetRandom(seed,p1,p2);
	      }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

/*
double Integral::Midpoint( unsigned int nstep){
	m_h = ( m_b - m_a ) / nstep;
	int n = 0;
	while( n < nstep ){
		m_sum += (m_f->Eval(m_a + ( n + 0.5)*m_h ) )* m_h;
		n++;
	}
	return m_integral = m_sign*m_sum;
}

double Integral::Simpson( unsigned int nstep){
	if ( nstep % 2 == 1 ){
		cout << "Numero dispari!!" << endl;
		exit(-1);
	}
	m_h = ( m_b - m_a ) / nstep;
	int n = 1;
	while ( n  < nstep ){
		if ( n % 2 == 1 ){
			m_sum += 4.*m_f->Eval(m_a + n*m_h);
		} else {
			m_sum += 2.*m_f->Eval(m_a + n*m_h); 
		}
		n++;
	}
	m_sum += ( m_f->Eval(m_a) + m_f->Eval(m_b) );
	m_integral = m_sign*(m_h/3.)*m_sum;
	return m_integral;
}

double Integral::Trapezi( unsigned int nstep ){
	m_h = (m_b - m_a) / nstep;
	m_sum = 0.5*m_h*(m_f->Eval(m_a) + m_f->Eval(m_b));
	for( int count = 1; count < nstep; count++)
		m_sum+=m_f->Eval(m_a + count*m_h);
	m_sum *= m_h;
	return m_integral = m_sign*m_sum;
}

double Integral::Trapezi_prec( double prec ){
	double sum = m_f -> Eval(m_a) + m_f -> Eval(m_b);
	double integral = sum*(m_b-m_a);
	m_sum = sum; m_integral = integral;
	int i = 1;
	double eps = prec + 1;
	while(eps > prec){
		integral = m_integral;
		sum = m_sum;
		m_h = (m_b-m_a)/(double)pow(2,i);
		for(int j = 0; j < pow(2,i-1); j++)
			m_sum += m_f->Eval(m_h + 2*j*m_h);
		m_integral = m_sum * (m_b-m_a) / (double) pow(2,i);	
		eps = 4./3.*fabs(m_integral - integral);
		i++;
	}
	m_integral = integral * m_sign;
	return m_integral;
}
*/
//VERSIONE PESANTE
/*double Integral::Trapezi_prec( double prec ){
	int n = 1;
	double int1 = 0;
	double int2 = 0;
	double err = 2*prec;
	Integral * Integ = new Integral( m_a, m_b, m_f);
	while ( err > prec ){
		int1 = Integ->Trapezi(n);
		int2 = Integ->Trapezi(2*n);
		err = (4./3.)*fabs( int2 - int1);	
		n++;
	}
	m_integral = m_sign*int1;
	return m_integral;
}*/

double Integral::Media(int nstep ){
	double sum = 0;
	for( int i = 0; i < nstep; i++ )
		sum += m_f -> Eval( m_rnd.Rannyu(m_a, m_b) );
	m_integral = m_sign * sum / (double) nstep * ( m_b - m_a );
	return m_integral;
}

double Integral::Mediaretta(int nstep ){
	double sum = 0;
	for( int i = 0; i < nstep; i++ )
		sum += m_f -> Eval( m_rnd.Retta() );
	m_integral = m_sign * sum / (double) nstep * ( m_b - m_a );
	return m_integral;
}
/*
double Integral::HitMiss( unsigned int nstep, double fmax ){
	double x = 0, y = 0;
	unsigned int N = 0; 
	for( int i = 0; i < nstep; i++ ){
		x = m_myrand->Unif(m_a, m_b);
		y = m_myrand->Unif(0, fmax);
		if ( y < m_f -> Eval(x) )
			N++;
	}
	m_integral = m_sign * fmax*(m_b - m_a)*N/nstep;
	return m_integral;
}

double Integral::Media_prec( double prec, double fmed ){
	unsigned int N = 2;
	double sum = 0;
	double eps  = prec + 2;
	while ( eps > prec ){
		for(int j = 0; j < N; j++){
			sum += pow( m_f -> Eval( m_myrand -> Unif( m_a, m_b) ) - fmed, 2);
		}
		eps = sum / (N - 1);
		eps = sqrt( fabs( eps/N*(m_b - m_a) ) );
		cout << N << "\t" << eps << endl;
		N = N+10000;
		sum = 0;
	}
	cout << eps << endl;
	m_integral = Media(N);
	return m_integral;
}
*/
