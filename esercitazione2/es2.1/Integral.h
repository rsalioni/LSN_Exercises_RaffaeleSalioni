#ifndef __Integral_h__
#define __Integral_h__

#include<iostream>

#include"Funzioni.h"
#include"random.h"

class Integral {
	
	public:
		Integral() {m_a = 0; m_b = 0; m_f = nullptr;};
		//Integral( double a, double b, const FunzioneBase * f );
		Integral( double a, double b, const FunzioneBase * f);
		~Integral(){;};
		
		//double Midpoint( unsigned int nstep );
		//double Simpson( unsigned int nstep );
		//double Trapezi( unsigned int nstep );
		//double Trapezi_prec( double prec );
		double Media(int nstep );
		double Mediaretta(int nstep );
		//double Media_prec( double prec, double fmed );
		//double HitMiss( unsigned int nstep, double fmax );


	private:
		double m_a, m_b;
		const FunzioneBase * m_f;
		double m_integral;
		double m_sum;
		double m_h;
		int m_sign;
		Random m_rnd;
		
};

#endif
