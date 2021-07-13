#ifndef __Funzione_Base__
#define __Funzione_Base__

#include<cmath>

class FunzioneBase {
	public:
		virtual double Eval( double x ) const = 0;
		
};

class Parabola : public FunzioneBase{
	public:
		
		//costruttori
		Parabola() {m_a = 0; m_b = 0; m_c = 0;};
		~Parabola() {};
		Parabola( double a, double b, double c ) { m_a = a; m_b = b; m_c = c;};
		
		//setparametri
		void SetA( double a ) { m_a = a; };
		double GetA() const { return m_a; };
		void SetB( double a ) { m_b = a; };
		double GetB() const { return m_b; };
		void SetC( double a ) { m_c = a; };
		double GetC() const { return m_c; };
	
		virtual double Eval (double x) const { return m_a*x*x+m_b*x+m_c; };
	
	
	private:
		double m_a, m_b, m_c;

};

class Coseno : public FunzioneBase{
	public:
		Coseno(double a, double b, double c){m_a = a; m_b=b; m_c=c;};
		~Coseno(){;};
		virtual double Eval( double x ) const { return m_a*cos(m_b*x+m_c); };
	private:
	double m_a, m_b, m_c;		
};

class Myfun : public FunzioneBase{
	public:
		Myfun() {};
		~Myfun(){;};
		virtual double Eval( double x ) const { return ((M_PI/2.)*cos((M_PI/2.)*x))/(2*(1-x)); };
};

class Solutore {
	public:
		virtual double CercaZeri( double xmin, double xmax, const FunzioneBase * f, double eps, double nmax ) = 0;
		virtual bool Trovato() = 0;
		virtual double Incertezza() = 0;
	
	protected:
		double m_a, m_b;
		bool found;
		double prec_eff;

};

class Bisezione : public Solutore {
	public:
		Bisezione(){};
		~Bisezione(){};
		
		virtual double CercaZeri( double xmin, double xmax, const FunzioneBase * f, double eps, double nmax );
		virtual bool Trovato() { return found; };
		virtual double Incertezza() {return prec_eff; };
};

#endif
