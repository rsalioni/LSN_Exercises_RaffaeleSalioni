#include<iostream>
#include<cstdlib>
#include<cmath>
#include"Funzioni.h"

using namespace std;

/*
double Bisezione::CercaZeri( double xmin, double xmax, const FunzioneBase * f, double eps, double nmax ){
	if(xmin < xmax){
		m_a=xmin;
		m_b=xmax;
	} else {
		m_a=xmax;
		m_b=xmin;
	}
	if (f->Eval(m_a)*f->Eval(m_b) > 0){
		cerr << "no teo zeri" << endl;
		found = false;
		exit(-1);
	}
	int n=0;
	double prec = m_b - m_a;
	double t = (m_b+m_a)/2;
	while(n < nmax && prec > eps){
	if(t==0){
		found = true;
		return t;
	}
	if (f->Eval(m_a)*f->Eval(t) < 0){
		m_b = t;
	}
	else { 
	m_a = t;
	}
	t= (m_b+m_a) / 2;
	n++;
	prec = m_b - m_a;
	}
	prec_eff = 0.5 * prec;
	found = true;
	return t;
}

double Secante::CercaZeri( double xmin, double xmax, const FunzioneBase * f, double eps, double nmax ){
	if ( xmax < xmin ) {
		double temp = xmax;
		xmax = xmin;
		xmin = temp;
	}
	double x3 = xmax - ( xmax - xmin )/( f->Eval(xmax) - f->Eval(xmin) )*( f -> Eval(xmax) );
	double delta = (xmax - xmin ) / 2;
	int i = 0;
	double vero = ( xmax + xmin ) / 2;;
	while ( ( i < nmax ) && ( delta > eps ) ) {
		if ( x3*f->Eval(xmax) < 0) {
			double temp = x3 - ( x3 - xmax )/( f->Eval(x3) - f->Eval(xmax) )*( f -> Eval(x3) );
			xmin = xmax;
			xmax = x3;
			x3 = temp;
		}	else {
			double temp = x3 - ( x3 - xmin )/( f->Eval(x3) - f->Eval(xmin) )*( f -> Eval(x3) );
			xmax =x3;
			x3 = temp;									
		}
		delta = abs( xmax - xmin ) / 2;
		i++;
		vero = ( xmax + xmin ) / 2;
	}
	return vero;
};
*/
