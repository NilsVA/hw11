#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx, const double omega, const double xmin);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.;
  const double xmax = 40.;
	const double Tend = 10. * M_PI;
	const double dx = (xmax-xmin)/Nx;
	const double dt =  dx/20. ;
  double t = 0.;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10.;
  const double omega = 0.2;
        const double alpha = sqrt(omega);
        
  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
        cmplx* psi1 = new cmplx[Nx];
        cmplx* h;
	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) 
                {
                 step(psi0,psi1,dt,dx,Nx,omega,xmin);
                 h=psi0;
                 psi0 = psi1;
                 psi1= h;
                 t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}

void step(cmplx* const psi0, cmplx* const psi1, const double dt,
	  const double dx, const int Nx, const double omega, const double xmin)
{
  
const cmplx a = cmplx(0.0,1.0)*dt/(dx*dx*4.) * (-1.0);
const cmplx a2 = a*a;
cmplx* d=new cmplx[Nx];
double* V = new double[Nx];

for(int i=0;i<Nx;i++)
{
    V[i] = (omega*omega)/2 * (xmin + i*dx)*(xmin + i*dx);
    d[i] = 1.0 + cmplx(0.0,1.0)*((-1.)*2.*a + (dt/2.) * V[i]);
}
for(int i=1;i<Nx;i++)
  {
      
  d[i] -= a2/d[i-1];
  psi0[i] -= -a/d[i-1]*psi0[i-1];

  }
  
  psi1[Nx-1] = psi0[Nx-1]/d[Nx-1];

  for(int i=Nx-2;i>0; i--)

  psi1[i] = (psi0[i] + a*psi1[i+1])/d[i];

  delete[] V;
  delete[] d;
}