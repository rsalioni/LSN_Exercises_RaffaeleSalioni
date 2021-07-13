#include<iostream>
#include<cmath>
#include<cstdlib>
#include <fstream>
#include <string>
#include"Geometria.h"
using namespace std;
// Posizione

Posizione::Posizione() {
	m_x = 0;
	m_y = 0;
	m_y = 0;
	m_n = 0;
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
};

Posizione::Posizione(double x, double y, double z) {
	m_x = x;
	m_y = y;
	m_z = z;
	m_n = 0;
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
};

Posizione::~Posizione() {

};

double Posizione::getR() const {
	return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
};

double Posizione::getTheta() const {
	return acos(m_z/(sqrt(m_x*m_x + m_y*m_y + m_z*m_z)));
};

double Posizione::getPhi() const {
	return asin(m_y/sqrt(m_x*m_x + m_y*m_y));
};

double Posizione::Distanza(const Posizione& p) const {
	return sqrt( (m_x - p.m_x)*(m_x - p.m_x) + (m_y - p.m_y)*(m_y - p.m_y) + (m_z - p.m_z)*(m_z - p.m_z));
};

double Posizione::getRho() const {
	return sqrt( m_x*m_x + m_y*m_y );
};



void Posizione::RW_continuo(double a){
	m_theta=m_rnd.Rannyu(0., M_PI);
	m_phi=m_rnd.Rannyu(0.,2*M_PI);
	m_x += a*sin(m_theta)*cos(m_phi);
	m_y += a*sin(m_theta)*sin(m_phi);
	m_z += a*cos(m_theta);
};

void Posizione::RW_discreto(double a){
				m_temp=m_rnd.Rannyu();
				if(m_temp < 1./6.)
					m_x += a;
				else if(1./6. < m_temp and m_temp < 2./6.)
					m_x = m_x - a;
				else if(2./6. < m_temp and m_temp < 3./6.)
					m_y = m_y + a;
				else if(3./6. < m_temp and m_temp < 4./6.)
					m_y = m_y - a;
				else if(4./6. < m_temp  and m_temp < 5./6.)
					m_z = m_z + a;
				else if( m_temp > 5./6.)
					m_z = m_z -a;
};


void Posizione::Metropolis_unif(FunzioneScalareBase* f){
	m_x1 = m_rnd.Rannyu(-m_d/2, m_d/2) + m_x;
	m_y1 = m_rnd.Rannyu(-m_d/2, m_d/2) + m_y;
	m_z1 = m_rnd.Rannyu(-m_d/2, m_d/2) + m_z;
	m_a = fmin(1, f->Eval(m_x1, m_y1, m_z1 )/(f->Eval(m_x, m_y, m_z)));
	
	m_temp = m_rnd.Rannyu();
	//cout << "a " <<  m_a << "temp "<< m_temp << endl;
	if( m_temp < m_a){
		m_x = m_x1;
		m_y = m_y1;
		m_z = m_z1;
		m_n++;
	}
};


/*

// Campo Vettoriale

CampoVettoriale::CampoVettoriale(const Posizione& p) : Posizione(p) {
	m_Fx = 0;
	m_Fy = 0;
	m_Fz = 0;
};

CampoVettoriale::CampoVettoriale(const Posizione& p, double Fx, double Fy, double Fz ) : Posizione(p) {
	m_Fx = Fx;
	m_Fy = Fy;
	m_Fz = Fz;
};

CampoVettoriale::CampoVettoriale(double x, double y, double z, double Fx, double Fy, double Fz ) : Posizione(x, y, z) {
	m_Fx = Fx;
	m_Fy = Fy;
	m_Fz = Fz;
};

CampoVettoriale CampoVettoriale::operator+( const CampoVettoriale & v ) const{
	CampoVettoriale sum(v);
	sum.m_Fx += getFx();
	sum.m_Fy += getFy();
	sum.m_Fz += getFz();
	return sum;
};

CampoVettoriale & CampoVettoriale::operator+=( const CampoVettoriale & v) {
	return (*this) = (*this)+v;
};

double CampoVettoriale::Modulo(){
	return sqrt(m_Fx*m_Fx + m_Fy*m_Fy + m_Fz*m_Fz);
};

void CampoVettoriale::Print(){
	std::cout << m_Fx << "\t" << m_Fy << "\t" << m_Fz << std::endl;
}

*/
