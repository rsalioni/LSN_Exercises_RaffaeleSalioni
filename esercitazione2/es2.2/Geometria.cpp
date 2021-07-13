#include<iostream>
#include<cmath>
#include<cstdlib>

#include"Geometria.h"

// Posizione

Posizione::Posizione() {
	m_x = 0;
	m_y = 0;
	m_y = 0;
};

Posizione::Posizione(double x, double y, double z) {
	m_x = x;
	m_y = y;
	m_z = z;
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
