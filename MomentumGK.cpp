//**********************************************************************************************************************
//**********************************************************************************************************************
//
//   Program MomentumGK.cpp:
//
//   This program computes the spectral linear momentum dp/dw transferred to a 
//   NP due a fast swfit electron. We calculate the closed surface integral of
//   the Maxwell stress tensor via Gauss Konrod quadrature, the real part of
//   the the dpdw is the linear momentum, and its imaginary part 
//   the estimated error.   
//
//   The index in the dpdw array correpond to:
//
//    0 -> Electric int       1 -> Magnetic int
//    2 -> Electric scat      3 -> Magnetic scat
//    4 -> Electric ext       5 -> Magnetic ext
//
//
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (27/06/18)
//
//**********************************************************************************************************************
//**********************************************************************************************************************



#include <iostream>                                                // Standart i/o C++ library
#include <complex>                                                 // Compĺex numbers
#include <boost/math/special_functions/bessel.hpp>                 // BOOST LIBRARIES:  1. BesselK in external fields
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>               // Lengendre Plm
#include <boost/math/quadrature/gauss_kronrod.hpp>                 // Gauss Konrod Quadrature for surface integral of T.da
#include <fstream>

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.


double wp     = 0.555;
double Gamma  = 0.00555;
double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;



//************************************************
// Maximum number of multipole to take into account

const int Lmax = 3;


// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|

const int Nt = 25;       // Must be odd integer!!!
const int Nf = 31;       

//************************************************



//**********************************************************************
// *********************** Dielectric function ***************************
//**********************************************************************
dcomplex eps(double w){
return 1. - pow(wp,2.)/(w*(w + 1i*Gamma));
} 

//**********************************************************************
// Common functions
//********************************************************************** 

double BesselK(int n, double x){
	return boost::math::cyl_bessel_k(n,x);
}                                                  

double LP(int l, int m, double x){
	return boost::math::legendre_p(l,m,x);
}

double factorial2(int n)
{
    if (n == 0 || n==1){
      return 1.;
    }else if(n > 1){
    return n*factorial2(n-2);
    }
    else{
        return pow(-1.,(n-1.)/2.)*n/factorial2(abs(n));
    }
}


double factorial(int n)
{
    if (n == 0 || n==1){
      return 1.;
    }else{
    return n*factorial(n-1);
    }
}


//********************************************************************************
//*********************** External fields (cartesian) **********************
//********************************************************************************

dcomplex  Ez(double w, double beta, double b, double x, double y, double z){
double v = beta*Cspeed;
double gam = 1./sqrt(1.0 - pow(beta,2.0));
return 2.0*1i*w*exp(1i*w*z/v)*(pow(v,-2.0) - pow(Cspeed,-2.0))*
BesselK(0,w*sqrt(pow(y,2.) + pow(x-b,2.0))/(v*gam));
}


dcomplex   Ey(double w, double beta, double b, double x, double y, double z){
double v = beta*Cspeed;
double gam = 1./sqrt(1.0 - pow(beta,2.0));
return -2.0*y*exp(1i*w*z/v)*w/(gam*pow(v,2.0)*sqrt(pow(y,2.0) + pow(x-b,2.0)))*
BesselK(1,w*sqrt(pow(y,2.) + pow(x-b,2.0))/(v*gam));
}


dcomplex   Ex(double w, double beta, double b, double x, double y, double z){
double v = beta*Cspeed;
double gam = 1./sqrt(1.0 - pow(beta,2.0));
return -2.0*(x-b)*exp(1i*w*z/v)*w/(gam*pow(v,2.0)*sqrt(pow(y,2.) + pow(x-b,2.0)))*
BesselK(1,w*sqrt(pow(y,2.) + pow(x-b,2.0))/(v*gam));
}


dcomplex   Hy(double w, double beta, double b, double x, double y, double z){
double l = 2.*Pi*Cspeed/w;
double v = beta*Cspeed;
double gam = 1./sqrt(1.0 - pow(beta,2.0));
return -4.0*Pi*(x-b)*exp(1i*w*z/v)/(gam*v*l*sqrt(pow(y,2.) + pow(x-b,2.0)))*
BesselK(1,w*sqrt(pow(y,2.) + pow(x-b,2.0))/(v*gam));
}

dcomplex   Hx(double w, double beta, double b, double x, double y, double z){
double l = 2.*Pi*Cspeed/w;
double v = beta*Cspeed;
double gam = 1./sqrt(1.0 - pow(beta,2.0));
return 4.0*Pi*y*exp(1i*w*z/v)/(gam*v*l*sqrt(pow(y,2.) + pow(x-b,2.0)))*
BesselK(1,w*sqrt(pow(y,2.) + pow(x-b,2.0))/(v*gam));
}
//********************************************************************************



//********************************************************************************
// ****************** Gauss Konrod Quadrature weights and abscisas ********************
//********************************************************************************
// GKQt for theta, and GKQf for phi

void gauss_konrod(double GKQt[2*Nt + 1][3], double GKQf[2*Nf + 1][3]){
   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*Nt+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*Nt+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, Nt>::weights();

   auto XGKf = boost::math::quadrature::gauss_kronrod<double, 2*Nf+1>::abscissa();
   auto WGKf = boost::math::quadrature::gauss_kronrod<double, 2*Nf+1>::weights();
   auto WGf =  boost::math::quadrature::gauss<double, Nf>::weights();

   double WGLt[2*Nt+1], WGLf[2*Nf+1];


// Changes the n order array of Gaussian quadrature to a 2n+1 array

for(int i = 0; i < XGKt.size(); ++i){
      if (i % 2 == 0){
           WGLt[i] = WGt[i/2];}
      else WGLt[i] = 0.;
}

for(int i = 0; i < XGKf.size(); ++i){
      if (i % 2 == 0){
           WGLf[i] = WGf[i/2];}
      else WGLf[i] = 0.;
}

// writtes [0] abscisas, [1] Konrod weights, [2] Gauss weigths

for(int i = 0; i < 2*XGKt.size() - 1; ++i){
      if (i <= Nt){ 
         GKQt[i][0] = -1.*XGKt[Nt - i];
         GKQt[i][1] = WGKt[Nt - i];                
         GKQt[i][2] = WGLt[Nt - i];
      }
      else{
         GKQt[i][0] = XGKt[-Nt + i];
         GKQt[i][1] = WGKt[-Nt + i];
         GKQt[i][2] = WGLt[-Nt + i];
      }
}

for(int i = 0; i < 2*XGKf.size() - 1; ++i){
      if (i <= Nf){ 
         GKQf[i][0] = -1.*XGKf[Nf - i];
         GKQf[i][1] = WGKf[Nf - i];
         GKQf[i][2] = WGLf[Nf - i];
      }
      else{
         GKQf[i][0] = XGKf[-Nf + i];
         GKQf[i][1] = WGKf[-Nf + i];
         GKQf[i][2] = WGLf[-Nf + i];
      }
}

}
//***********************************************************************



//**************************************************************************************
// Recursive functions and parameter needed for calculate the scatterd fields
//*************************************************************************************

double II(int l, int m, int i1, int i2){

if (l == m - 2 || l == m - 1){return 0.;} 
     else if(l == m){ 	if (i2 % 2 == 0){
     		               return pow(-1.,m)*factorial2(2*m - 1)*boost::math::beta((i1 + m + 2)/2., (i2 + 1)/2.,1.);}
     		            else{return 0.;}}
     	else{return (1./(l - m))*((2.*l - 1.)*II(l - 1, m,i1, i2 + 1) - (l + m - 1.)*II(l - 2, m,i1, i2));}
              }


double alm(int l, int m){
	if (abs(m)<=l)
	{
		return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
	}
	else{return 0.;}
}


dcomplex  A(int l, int m, double betav){

double gam = 1./sqrt(1.0 - pow(betav,2.0));
dcomplex  res = 0.;

if(m >= 0){
for (int j = m; j <= l; ++j){
	res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*alm(l,m)*II(l,m,j,l-j)
	/(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
	}
return res;
 } 
else {
 	return pow(-1.,abs(m))*A(l,abs(m),betav);
 }
}



dcomplex  B(int l, int m, double betav){
return A(l,m+1,betav)*sqrt((l+m+1)*(l-m)) - A(l,m-1,betav)*sqrt((l-m+1)*(l+m));
}

dcomplex  PsiE(int l, int m, double w, double b,double betav){
double gam = 1./sqrt(1.0 - pow(betav,2.0));
return -2.*Pi*pow(1i,1-l)*w*B(l,m,betav)*BesselK(m,w*b/(betav*Cspeed*gam))/(pow(Cspeed,2.)*gam*l*(l+1));
}

dcomplex PsiM(int l, int m, double w, double b,double betav){
double gam = 1./sqrt(1.0 - pow(betav,2.0));
return -(4.*m)*Pi*pow(1i,1-l)*betav*w*A(l,m,betav)*BesselK(m,w*b/(betav*Cspeed*gam))/(pow(Cspeed,2)*l*(l+1));
}


dcomplex j(int n, dcomplex z){    
return sph_besselJ(n,z);
}

dcomplex h(int n, dcomplex z){                              
return sph_hankelH1(n,z);
}


dcomplex djj(int n, dcomplex x){
	return j(n,x) + 0.5*x*(j(n-1,x)-j(n+1,x) - j(n,x)/x);
}

dcomplex dhh(int n, dcomplex x){
	return h(n,x) + 0.5*x*(h(n-1,x)-h(n+1,x) - h(n,x)/x);
}


// *********** Polarizabilities *********************
dcomplex tE(int l, double w, double a){
double x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -1i*( eps(w)*j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - eps(w)*j(l,xi)*dhh(l,x0) ); 
}

dcomplex tM(int l, double w, double a){
double x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -1i*( j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - j(l,xi)*dhh(l,x0) );
}



dcomplex fs(int l, dcomplex z){
double dl = 1.*l;
	return (dl + 1.)*1i*h(l,z)/z  - 1i*h(l+1,z);
}
//******************************************************************************************************
//******************************************************************************************************



//*********************************************************************************************************
// Constructs the maxwell stress tensor and calculates its surface integral, writes to dpdwx (transversal)
// and dpdwz (longitudinal). Real part is its value and imaginary its error (konrod integration). 
// [0] Electric contribution. [1] Magnetic contribution.
//*********************************************************************************************************
//*********************************************************************************************************
void dpdw(double r, double w, double vv, double b, double a, dcomplex dpdwx[6], dcomplex dpdwz[6]){

double k0 = w/Cspeed;
double dm;
double dl;

dpdwx[0] = 0.;
dpdwz[0] = 0.;
dpdwx[1] = 0.;
dpdwz[1] = 0.;

dpdwx[3] = 0.;
dpdwz[3] = 0.;
dpdwx[2] = 0.;
dpdwz[2] = 0.;

dpdwx[4] = 0.;
dpdwz[4] = 0.;
dpdwx[5] = 0.;
dpdwz[5] = 0.;

dcomplex CM[Lmax][2*Lmax +1];
dcomplex DE[Lmax][2*Lmax +1];

dcomplex dle, clm;
dcomplex dhl[Lmax], dfl[Lmax];

// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    dle = pow(1i,l)*tE(l,w,a);
    clm = pow(1i,l)*tM(l,w,a);
    dhl[l-1] = h(l,k0*r);
    dfl[l-1] = fs(l,k0*r);


    for (int m = -l; m <= l; ++m){
    
        DE[l-1][m + Lmax] = alm(l,m)*dle*PsiE(l,m,w,b,vv);
        CM[l-1][m + Lmax] = alm(l,m)*clm*PsiM(l,m,w,b,vv);
    
    }
}   


dcomplex Ers, Ets, Efs;
dcomplex Hrs, Hts, Hfs;

dcomplex Eex, Eey, Eez;
dcomplex Hey, Hex, Hez;

dcomplex Esx, Esy, Esz, Hsx, Hsy, Hsz;

double TEesx, TEesz, THesx, THesz;
double TEssx, TEssz, THssx, THssz;
double TEeex, TEeez, THeex, THeez;
double theta, phi;

// Calls the Gauss - Konrod function

double GKQf[2*Nf + 1][3], GKQt[2*Nt + 1][3];
gauss_konrod(GKQt,GKQf);


//************************** GENERATION OF GRID *******************************

for(int j=0;  j < 2*Nt + 1;  ++j) {
      theta = Pi*(GKQt[j][0] + 1.0)/2.0;

    for(int k=0; k < 2*Nf + 1;  ++k)  {
       phi= Pi*(GKQf[k][0] + 1.0);

Ers = 0.;
Ets = 0.;
Efs = 0.;

Hrs = 0.;
Hts = 0.;
Hfs = 0.;




// ************** Scatterrd fields ****************
for (int l1 = 1; l1 <= Lmax; ++l1){
        dl = 1.*l1;
        int l = l1; 
        for (int m = -l; m <= l; ++m){
          dm = 1.*m;
          
             Ers += exp(1i*dm*phi)*DE[l-1][m+Lmax]*dl*(dl+1.)*LP(l,m,cos(theta))*(1i*dhl[l-1])/(k0*r);

             Ets += -exp(1i*dm*phi)*(CM[l-1][m+Lmax]*LP(l,m,cos(theta))*(1i*dm*dhl[l-1])/sin(theta) + 
                      DE[l-1][m+Lmax]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*dfl[l-1]);

             Efs += 1i*exp(1i*dm*phi)*(DE[l-1][m+Lmax]*LP(l,m,cos(theta))*(dm*dfl[l-1])/sin(theta) + 
                      CM[l-1][m+Lmax]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*1i*dhl[l-1]);



             Hrs += exp(1i*dm*phi)*CM[l-1][m+Lmax]*dl*(dl+1.)*LP(l,m,cos(theta))*(1i*dhl[l-1])/(k0*r);

             Hts += exp(1i*dm*phi)*(DE[l-1][m+Lmax]*LP(l,m,cos(theta))*(1i*dm*dhl[l-1])/sin(theta) - 
                      CM[l-1][m+Lmax]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*dfl[l-1]);

             Hfs += 1i*exp(1i*dm*phi)*(CM[l-1][m+Lmax]*LP(l,m,cos(theta))*(dm*dfl[l-1])/sin(theta) - 
                      DE[l-1][m+Lmax]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*1i*dhl[l-1]);
        }

}




// ************ External fields ******************
Eex = Ex(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Eey = Ey(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Eez = Ez(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));

Hex = Hx(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Hey = Hy(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Hez = 0.;


Esx = Ers*sin(theta)*cos(phi) + Ets*cos(theta)*cos(phi) - Efs*sin(phi);
Esy = Ers*sin(theta)*sin(phi) + Ets*cos(theta)*sin(phi) + Efs*cos(phi);
Esz = Ers*cos(theta)          - Ets*sin(theta);

Hsx = Hrs*sin(theta)*cos(phi) + Hts*cos(theta)*cos(phi) - Hfs*sin(phi);
Hsy = Hrs*sin(theta)*sin(phi) + Hts*cos(theta)*sin(phi) + Hfs*cos(phi);
Hsz = Hrs*cos(theta)          - Hts*sin(theta);

// ***************+ Maxwell stress tensor *********************

TEesx = ((Eex*conj(Esx) - Eey*conj(Esy) - Eez*conj(Esz))*sin(theta)*cos(phi) + (Eex*conj(Esy) + Esx*conj(Eey))*sin(theta)*sin(phi)  +    (Eex*conj(Esz) + Esx*conj(Eez))*cos(theta)).real();
THesx = ((Hex*conj(Hsx) - Hey*conj(Hsy) - Hez*conj(Hsz))*sin(theta)*cos(phi) + (Hex*conj(Hsy) + Hsx*conj(Hey))*sin(theta)*sin(phi)  +    (Hex*conj(Hsz) + Hsx*conj(Hez))*cos(theta)).real();

TEesz = ((Eez*conj(Esz) - Eey*conj(Esy) - Eex*conj(Esx))*cos(theta) + (Eez*conj(Esy) + Esz*conj(Eey))*sin(theta)*sin(phi)  +    (Eez*conj(Esx) + Esz*conj(Eex))*sin(theta)*cos(phi)).real();
THesz = ((Hez*conj(Hsz) - Hey*conj(Hsy) - Hex*conj(Hsx))*cos(theta) + (Hez*conj(Hsy) + Hsz*conj(Hey))*sin(theta)*sin(phi)  +    (Hez*conj(Hsx) + Hsz*conj(Hex))*sin(theta)*cos(phi)).real();

//************************************


TEssx = (0.5*(Esx*conj(Esx) - Esy*conj(Esy) - Esz*conj(Esz))*sin(theta)*cos(phi) + Esx*conj(Esy)*sin(theta)*sin(phi)  +    Esx*conj(Esz)*cos(theta)).real();
THssx = (0.5*(Hsx*conj(Hsx) - Hsy*conj(Hsy) - Hsz*conj(Hsz))*sin(theta)*cos(phi) + Hsx*conj(Hsy)*sin(theta)*sin(phi)  +    Hsx*conj(Hsz)*cos(theta)).real();

TEssz = (0.5*(Esz*conj(Esz) - Esy*conj(Esy) - Esx*conj(Esx))*cos(theta) + Esz*conj(Esy)*sin(theta)*sin(phi)  +    Esz*conj(Esx)*sin(theta)*cos(phi)).real();
THssz = (0.5*(Hsz*conj(Hsz) - Hsy*conj(Hsy) - Hsx*conj(Hsx))*cos(theta) + Hsz*conj(Hsy)*sin(theta)*sin(phi)  +    Hsz*conj(Hsx)*sin(theta)*cos(phi)).real();


//************************************

TEeex = real(0.5*(Eex*conj(Eex) - Eey*conj(Eey) - Eez*conj(Eez))*sin(theta)*cos(phi) + 
                  Eey*conj(Eex)*sin(theta)*sin(phi) +
                  Eez*conj(Eex)*cos(theta));

TEeez = real(0.5*(Eez*conj(Eez) - Eey*conj(Eey) - Eex*conj(Eex))*cos(theta) + 
                  Eez*conj(Eex)*sin(theta)*cos(phi) +
                  Eez*conj(Eey)*sin(theta)*sin(phi));

THeex = real(0.5*(Hex*conj(Hex) - Hey*conj(Hey))*sin(theta)*cos(phi) + Hey*conj(Hex)*sin(theta)*sin(phi));

THeez = real(0.5*( - Hey*conj(Hey) - Hex*conj(Hex))*cos(theta)  );

//***********************************

    // ************************* SURFACE CLOSED INTEGRAL ***************************   

        dpdwx[0] += TEesx*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwz[0] += TEesz*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwx[1] += THesx*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;  
        dpdwz[1] += THesz*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;

        dpdwx[2] += TEssx*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwz[2] += TEssz*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwx[3] += THssx*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;  
        dpdwz[3] += THssz*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;

        dpdwx[4] += TEeex*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwz[4] += TEeez*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;
        dpdwx[5] += THeex*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;  
        dpdwz[5] += THeez*(GKQf[k][1]*GKQt[j][1]*(1. + 1i) - 1i*GKQf[k][2]*GKQt[j][2])*sin(theta)/8.0;


     }
 }


}
//*********************************************************************************************************
//*********************************************************************************************************


//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

double vv = 0.5;
double w = 0.32;
double b = 1.5*nm;                        // Defines parameter electron speed vv, frequency a.u,  
double a = 1.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 1.05*nm;

dcomplex dpdwx[6], dpdwz[6];
dpdw(r, w, vv, b, a, dpdwx, dpdwz); 


cout.precision(17);
cout << endl;
cout << "Spectral contribution to the linear momentum (Gauss - Kronrod) :" << endl;
cout << endl;
cout << "Lmax = " << Lmax << endl;
cout << endl;
cout << "a = " << a/nm << "nm" << endl;
cout << "b = " << b/nm << "nm" << endl;
cout << "v = " << vv << "c" << endl;
cout << "w = " << w << " har" << endl;
cout << endl;
cout << "dpE/dpx : " << dpdwx[0] << endl;
cout << "dpH/dpx : " << dpdwx[1] << endl;                // prints result
cout << endl;
cout << "dpE/dpz : " << dpdwz[0] << endl;
cout << "dpH/dpz : " << dpdwz[1] << endl;
cout << endl;
cout << "dpEs/dpx : " << dpdwx[2] << endl;
cout << "dpHs/dpx : " << dpdwx[3] << endl;                // prints result
cout << endl;
cout << "dpEs/dpz : " << dpdwz[2] << endl;
cout << "dpHs/dpz : " << dpdwz[3] << endl;
cout << endl;
cout << "dpEe/dpx : " << dpdwx[4] << endl;
cout << "dpHe/dpx : " << dpdwx[5] << endl;
cout << endl;
cout << "dpEe/dpz : " << dpdwz[4] << endl;
cout << "dpHe/dpz : " << dpdwz[5] << endl;

}

//**************************************************** END ****************************************************