#ifndef __Funzione_Base__
#define __Funzione_Base__

#include<cmath>
#include<cstdlib>
#include<vector>


using namespace std;

class FunzioneBase {
	public:
		virtual double Eval( double x ) const = 0;
		
};
/*
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

class Seno : public FunzioneBase{
	public:
		Seno(){};
		~Seno(){;};
		virtual double Eval( double x ) const { return sin(x); };		
};
*/

class FunzioneScalareBase {
	public:
		virtual double Eval( double x, double y, double z ) const = 0;
};

class H1S : public FunzioneScalareBase{ //(1./sqrt(M_PI))*exp(-2*sqrt(x*x + y*y + z*z))
	public:
		H1S(){};
		~H1S() {};
		virtual double Eval( double x, double y, double z ) const { return pow((1./sqrt(M_PI))*exp(-sqrt(x*x + y*y + z*z)), 2); };
};

class H2P : public FunzioneScalareBase{ //(1./sqrt(M_PI))*exp(-2*sqrt(x*x + y*y + z*z))
	public:
		H2P(){};
		~H2P() {};
		virtual double Eval( double x, double y, double z ) const { return pow((1./sqrt(32*M_PI))*z*exp(-0.5*sqrt(x*x + y*y + z*z)), 2); };
};


/*
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

class Secante : public Solutore {
	public:
		Secante(){};
		~Secante(){};
		
		virtual double CercaZeri( double xmin, double xmax, const FunzioneBase * f, double eps, double nmax );
		virtual bool Trovato() { return found; };
		virtual double Incertezza() {return prec_eff; };
};
*/
#endif
