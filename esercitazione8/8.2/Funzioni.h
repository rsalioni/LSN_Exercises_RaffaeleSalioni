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


class PsiT : public FunzioneBase{
	public:
		PsiT(double mu, double sigma){
			m_mu = mu;
			m_sig = sigma;
		};
		~PsiT() {};
		virtual double Eval( double x ) const {
			return exp(- (pow(x-m_mu,2)) / (2*m_sig*m_sig)) + exp(- (pow(x+m_mu,2)) / (2*m_sig*m_sig));
		};
		
		//double GetMu() const { return m_mu; };
		//double GetSig() const { return m_sig; };
		
		void SetMu( double mu ) { m_mu = mu; };
		double GetMu() const { return m_mu; };
		void SetSig( double sig ) { m_sig = sig; };
		double GetSig() const { return m_sig; };		
		
	private:
		double m_mu, m_sig;
};


class V : public FunzioneBase{
	public:
		V(){};
		~V() {};
		virtual double Eval( double x ) const { return pow(x, 4) - (5./2.)*x*x; };
		
	};

class H :  public FunzioneBase{
	public:
		H(PsiT * Psi, V * V){
			m_Psi = Psi;
			m_V = V;
			m_mu = m_Psi -> GetMu();
			m_sig = m_Psi -> GetSig();
		};
		~H() {};
		virtual double Eval( double x ) const {
			double psi = m_Psi->Eval(x);
			double v = m_V->Eval(x);
			double psisec = -( exp(- (pow(x-m_mu,2)) / (2*m_sig*m_sig))*(-1. + (pow(x-m_mu,2)) / (m_sig*m_sig) ) + exp(- (pow(x+m_mu,2)) / (2*m_sig*m_sig))*(-1. + (pow(x+m_mu,2)) / (m_sig*m_sig) ) ) / (2*m_sig*m_sig);
		return psisec/psi + v ;
	};
	
		void SetMu( double mu ) { m_mu = mu; };
		double GetMu() const { return m_mu; };
		void SetSig( double sig ) { m_sig = sig; };
		double GetSig() const { return m_sig; };
	
	private:
		PsiT * m_Psi;
		V * m_V;		
		double m_mu, m_sig;
		
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
