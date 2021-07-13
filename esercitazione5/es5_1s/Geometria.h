#ifndef __Geometria_h__
#define __Geometria_h__

#include "random.h"
#include "Funzioni.h"
// Posizione

class Posizione {

public:

 	// costruttori
	Posizione();
	Posizione(double x, double y, double z); 
	// distruttore
	~Posizione();
	// metodi
	double getX() const { return m_x;};       // Coordinate cartesiane
	double getY() const { return m_y;};
	double getZ() const { return m_z;};
  	int getN() const { return m_n;};
  	
	void setX( double x) { m_x = x;};
	void setY( double y) { m_y = y;};  
	void setZ( double z) { m_z = z;};
	void setD( double d) { m_d = d;};
  	void setN(int n) { m_n = n;};
  	
  	
	double getR() const;       // Coordinate sferiche
	double getPhi() const;
	double getTheta() const;
	double getRho() const;  // raggio delle coordinate cilindriche
	   
  
	double Distanza(const Posizione&) const; // distanza da un altro punto
	
	void RW_continuo(double a);
	void RW_discreto(double a);
	
	void Metropolis_unif(FunzioneScalareBase* f);
	
private:

	double m_x, m_y, m_z, m_theta, m_phi, m_temp;
	double m_x1, m_y1, m_z1;
	double m_d; //larghezza distribuzione unif
	Random m_rnd;  
	double m_a;
	int m_n;
};


/*
// Campo Vettoriale

class CampoVettoriale : public Posizione {

public:

  CampoVettoriale(const Posizione&);
  CampoVettoriale(const Posizione&, double Fx, double Fy, double Fz);
  CampoVettoriale(double x, double y, double z, double Fx, double Fy, double Fz);
  ~CampoVettoriale() {};

  CampoVettoriale & operator+=( const CampoVettoriale & ) ;
  CampoVettoriale operator+( const CampoVettoriale & ) const;
  
  double getFx() const {return m_Fx;}
  double getFy() const {return m_Fy;}
  double getFz() const {return m_Fz;}

  double Modulo();
  
  void Print();

private:

  double m_Fx, m_Fy, m_Fz;

};

*/

#endif // __Geometria_h__





