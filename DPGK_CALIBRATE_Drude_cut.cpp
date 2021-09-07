//*****************************************************************************************************************
//
//   DPGK_CALIBRATE.cpp: NOT TESTED
//
//    This program calculates the spectral linear momentum dp/dw transferred to a 
//   NP due a fast swfit electron. It calculate the closed surface integral of
//   the Maxwell stress tensor via Gauss Konrod quadrature, the real part of               
//   the the dpdw void function represents the momentum, and the imaginary part 
//   the estimated error. The integral in frequencies is calculated via Gauss - Kronrod 
//   adaptative algorithm, for low frecuencies, and for high frecuencies an Exp Sinh
//   automatic algorithm is used. Both included in the BOOST libraries.
//   
//   Warning: The frecuency partition set in this program may not work very well in realistic eps(w) materials
//            or for very big NP > 10nm. This can be fixed by locating the resonance region.  
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (25/10/2018)
//
//     THIS PROGRAM IS USED FOR CALIBRATING MORE EFICIENT (ALTOUGH LESS ACCURATE) QUADRATURES AND OMEGA INTERVALS!!
//
//*****************************************************************************************************************
//*****************************************************************************************************************



#include <iostream>                                                // Standart i/o C++ library
#include <complex>                                                 // Compĺex numbers
#include <boost/math/special_functions/bessel.hpp>                 // BOOST LIBRARIES:  1. BesselK in external fields
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>               // Lengendre Plm
#include <boost/math/quadrature/gauss_kronrod.hpp>                 // Gauss Konrod Quadrature for surface integral of T.da
#include <fstream>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include "IN11.h"

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.


// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|

const int Lmax = 1;



//************************************************

double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;
double wp     = 0.555;
double Gamma  = 0.00555;

//************************************************

static double IM[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IN[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IU[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IV[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IW[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IX[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IY[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double IZ[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double ID[LSmax][LSmax][2*LSmax + 1][2*LSmax + 1];
static double III[LSmax][LSmax+1][LSmax+1];

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

dcomplex  PsiE(int l, int m, double w, double b,double betav, dcomplex BB){
double gam = 1./sqrt(1.0 - betav*betav);
return -2.*Pi*pow(1i,1-l)*w*BB*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*gam*l*(l+1));
}

dcomplex PsiM(int l, int m, double w, double b,double betav, dcomplex AA){
double gam = 1./sqrt(1.0 - betav*betav);
return -(4.*m)*Pi*pow(1i,1-l)*betav*w*AA*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*l*(l+1));
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

dcomplex fe(int l, dcomplex z){
double dl = 1.*l;
  return (dl + 1.)*j(l,z)/z  - j(l+1,z);
}
//******************************************************************************************************
//******************************************************************************************************


//*********************************************************************************************************
// Linear momentum spectrum:
//
// [0] Electric contribution. [1] Magnetic contribution.  Cross
// [2] Electric contribution. [3] Magnetic contribution.  Scat
// [4] Electric contribution. [5] Magnetic contribution.  Ext
//*********************************************************************************************************
//*********************************************************************************************************
void dpdw(double r, double w, double vv, double b, double a, double dpdwx[6], double dpdwz[6]){

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErrx[4], IErrz[4], IHrrx[4], IHrrz[4];
dcomplex IErtx[4], IErtz[4], IHrtx[4], IHrtz[4];
dcomplex IErf[4], IHrf[4];
dcomplex IEttx[4], IEttz[4], IHttx[4], IHttz[4];
dcomplex IEffx[4], IEffz[4], IHffx[4], IHffz[4];

double dl1, dl2, dm1, dm2;

double IN1, IV1, IW1, IW2, IW3, IU1, IU2, IU3, IU4;
double IM1, IZ1, IX1, IX2, IX3, IY1, IY2, IY3, IY4, ID1, ID2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);



dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}



double k0 = w/Cspeed;


// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = 1i*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(1i,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dpdwx[rr] = 0.;
dpdwz[rr] = 0.;
}



for (int l1 = 1; l1 <= Lmax; ++l1){
     dl1 = 1.*l1; 
for (int l2 = 1; l2 <= Lmax; ++l2){
         dl2 = 1.*l2; 
         for (int m1 = -l1; m1 <= l1; ++m1){
              dm1 = 1.*m1;
              for (int m2 = -l2; m2 <= l2; ++m2){
                   dm2 = 1.*m2;



                    if(m2 == m1+1 || m2 == m1-1){

                                    // Radial - Radial 

                        IN1 = IN[l1-1][l2-1][m1+l1][m2+l2];
                        IV1 = IV[l1-1][l2-1][m1+l1][m2+l2];

                        IW1 = IW[l1-1][l2-1][m1+l1][m2+l2];
                        IW2 = IW[l1][l2-1][m1+l1+1][m2+l2];   
                        IW3 = IW[l1-1][l2][m1+l1][m2+l2+1];

                        IU1 = IU[l1-1][l2-1][m1+l1][m2+l2];
                        IU2 = IU[l1][l2-1][m1+l1+1][m2+l2];
                        IU3 = IU[l1-1][l2][m1+l1][m2+l2+1];
                        IU4 = IU[l1][l2][m1+l1+1][m2+l2+1];

                        /* Integrants, 0 -> Scat
                                       1 -> Ext
                                       2 -> Ext-Scat 
                                       3 -> Scat-Ext (second is conjugated)*/


                    IErrx[0] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 

                    IHrrx[0] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 


                    IErrx[1] = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 

                    IHrrx[1] = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; 



                    IErrx[2] = DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G

                    IHrrx[2] = CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G



                    IErrx[3] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G

                    IHrrx[3] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G


                                   // Zenital - Radial

                    IErtx[0] = -(tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[0] =  (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);



                    IErtx[1] = -(DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[1] =  (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);



                    IErtx[2] = -(DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[2] =  (- CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);




                    IErtx[3] = -(tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                               + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx[3] =  (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);




                                    // Azimutal - Radial


                    IErf[0] =  (dm2-dm1)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[0] =  (dm2-dm1)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);  

                    IErf[1] =  (dm2-dm1)*( CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[1] =  (dm2-dm1)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);  

                    IErf[2] =  (dm2-dm1)*( CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[2] =  (dm2-dm1)*( -DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);

                    IErf[3] =  (dm2-dm1)*( tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf[3] =  (dm2-dm1)*( -tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                          + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);


                                        // Azimutal - Azimutal

                    IEffx[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //   

                    IEffx[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  



                    IEffx[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  

                    IEffx[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //  




                                        //  Zenital - Zenital

                    IEttx[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 



                    IEttx[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);  


                    IEttx[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 


                    IEttx[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IU1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); }


                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IErrx[rr] = 0.;
                     IHrrx[rr] = 0.;
                     IErtx[rr] = 0.;
                     IHrtx[rr] = 0.;
                     IEffx[rr] = 0.;
                     IHffx[rr] = 0.;
                     IErf[rr] = 0.;
                     IHrf[rr] = 0.;
                     IEttx[rr] = 0.;
                     IHttx[rr] = 0.;
                    }
                   }




                if(m2 == m1){
                                              // Radial - Radial 

                        IM1 = IM[l1-1][l2-1][m1+l1][m2+l2];
                        IZ1 = IZ[l1-1][l2-1][m1+l1][m2+l2];

                        ID1 = ID[l1-1][l2-1][m1+l1][m2+l2];
                        ID2 = ID[l1][l2-1][m1+l1+1][m2+l2];

                        IX1 = IX[l1-1][l2-1][m1+l1][m2+l2];
                        IX2 = IX[l1][l2-1][m1+l1+1][m2+l2];
                        IX3 = IX[l1-1][l2][m1+l1][m2+l2+1];

                        IY1 = IY[l1-1][l2-1][m1+l1][m2+l2];
                        IY2 = IY[l1][l2-1][m1+l1+1][m2+l2];
                        IY3 = IY[l1-1][l2][m1+l1][m2+l2+1];
                        IY4 = IY[l1][l2][m1+l1+1][m2+l2+1];



                    IErrz[0] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[0] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[1] = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[1] = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[2] = DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[2] = CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                    IErrz[3] = tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz[3] = tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                                              // Zenital - Radial

                    IErtz[0] = -( tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[0] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[1] = -( DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[1] = (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[2] = -( DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[2] = (- CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(df[1][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*(dz[1][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                    IErtz[3] = -( tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz[3] = (- tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(df[0][l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             +    tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*(dz[0][l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                                           // Azimutal - Azimutal

                    IEffz[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G




                    IEffz[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G




                    IEffz[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G





                    IEffz[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G



                                                                      //  Zenital - Zenital

                    IEttz[0] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[0] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 




                    IEttz[1] =  CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[1] =  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  





                    IEttz[2] =  CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[2] =  DE[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*dz[1][l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*dz[1][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - CM[l1-1][m1+l1]*conj(tEl[l2-1]*DE[l2-1][m2+l2]*dz[0][l2-1])*df[1][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + CM[l1-1][m1+l1]*conj(tMl[l2-1]*CM[l2-1][m2+l2]*df[0][l2-1])*df[1][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  


                    IEttz[3] =  tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz[3] =  tEl[l1-1]*DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*dz[0][l1-1]*dm1*dm2*IY1 
                              - tEl[l1-1]*DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*dz[0][l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - tMl[l1-1]*CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dz[1][l2-1])*df[0][l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + tMl[l1-1]*CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*df[1][l2-1])*df[0][l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); }  
                                              
                else{
                  for (int rr = 0; rr < 4; ++rr){
                     IErrz[rr] = 0.;
                     IHrrz[rr] = 0.;
                     IErtz[rr] = 0.;
                     IHrtz[rr] = 0.;
                     IEttz[rr] = 0.;
                     IHttz[rr] = 0.;
                     IEffz[rr] = 0.;
                     IHffz[rr] = 0.;
                    }
                   }




                     
                    dpdwx[0] += (1./(4.*Pi))*((0.5*(IErrx[2] - IEttx[2] - IEffx[2]) + IErtx[2] - IErf[2]).real()
                                            + (0.5*(IErrx[3] - IEttx[3] - IEffx[3]) + IErtx[3] - IErf[3])).real();
                    dpdwz[0] += (1./(2.*Pi))*((0.5*(IErrz[2] - IEttz[2] - IEffz[2]) - IErtz[2]).real()
                                            + (0.5*(IErrz[3] - IEttz[3] - IEffz[3]) - IErtz[3])).real();

                    dpdwx[1] +=  (1./(4.*Pi))*((0.5*(IHrrx[2] - IHttx[2] - IHffx[2]) + IHrtx[2] - IHrf[2]).real()
                                             + (0.5*(IHrrx[3] - IHttx[3] - IHffx[3]) + IHrtx[3] - IHrf[3])).real();
                    dpdwz[1] +=  (1./(2.*Pi))*((0.5*(IHrrz[2] - IHttz[2] - IHffz[2]) - IHrtz[2]).real()
                                             + (0.5*(IHrrz[3] - IHttz[3] - IHffz[3]) - IHrtz[3])).real();

                    dpdwx[2] +=  (1./(4.*Pi))*((0.5*(IErrx[0] - IEttx[0] - IEffx[0]) + IErtx[0] - IErf[0])).real();
                    dpdwz[2] +=  (1./(2.*Pi))*((0.5*(IErrz[0] - IEttz[0] - IEffz[0]) - IErtz[0])).real();

                    dpdwx[3] +=  (1./(4.*Pi))*((0.5*(IHrrx[0] - IHttx[0] - IHffx[0]) + IHrtx[0] - IHrf[0])).real();
                    dpdwz[3] +=  (1./(2.*Pi))*((0.5*(IHrrz[0] - IHttz[0] - IHffz[0]) - IHrtz[0])).real();

                    dpdwx[4] +=  (1./(4.*Pi))*((0.5*(IErrx[1] - IEttx[1] - IEffx[1]) + IErtx[1] - IErf[1])).real();
                    dpdwz[4] +=  (1./(2.*Pi))*((0.5*(IErrz[1] - IEttz[1] - IEffz[1]) - IErtz[1])).real();

                    dpdwx[5] +=  (1./(4.*Pi))*((0.5*(IHrrx[1] - IHttx[1] - IHffx[1]) + IHrtx[1] - IHrf[1])).real();
                    dpdwz[5] +=  (1./(2.*Pi))*((0.5*(IHrrz[1] - IHttz[1] - IHffz[1]) - IHrtz[1])).real();





            } // for m2
        }  // for m1

} // for l2
} //for l1

}



double dp(double b, double vv, double x){

double a = 1.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 1.05*nm;

double dpdwx[6], dpdwz[6];

if( BesselK(0,x) == 0. ){
  return 0.;}
else{
dpdw(r, x, vv, b, a, dpdwx, dpdwz);      // Calls function momentum dpdw
return (dpdwx[0] + dpdwx[1] + dpdwx[2] + dpdwx[3]);}
}



//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

double vv = 0.95;
double b = 1.5*nm;

double lf = 260.;  // 40, v0.3 b1.5;  65, v0.5 b1.5;  100, v0.3 b1.5;  160, v0.9 b1.5; 



auto ff = [b,vv](double x) { return dp(b,vv,x) ; };


cout << "Calibrating Al NP" << endl;        // NEEDED TO SET BY HAND !!!!
cout << "v = " << vv << endl;
cout << "b = " << b/nm << endl;
cout << "Lmax = "<< Lmax << endl;
cout << "wc "<< lf << " Har" << endl;
cout << endl;
cout.precision(17);

double l1 = 0.;
double l2 = .3;
double l3 = .4;
double l4 = 2.;
double l5 = std::numeric_limits<double>::infinity();
double termination = sqrt(std::numeric_limits<double>::epsilon());
double error;
double L1;


double Q  = boost::math::quadrature::gauss_kronrod<double, 101>::integrate(ff,  l1, l2, 5, 5.e-14, &error);
double Q1 = boost::math::quadrature::gauss_kronrod<double, 201>::integrate(ff, l2, l3, 5, 5.e-14, &error);
double Q2 = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(ff,  l3, l4, 5, 5.e-14, &error);



boost::math::quadrature::exp_sinh<double> integrator;
double Q3 = integrator.integrate(ff, l4, l5, termination, &error, &L1);


double II = Q + Q1 + Q2 + Q3;


cout << "I = " << II <<endl<< endl;





int i, NN = 24;
double qq[NN], QQ, lcut[NN];


for ( i = 1; i <= NN ; ++i){

lcut[i] = l4 + (lf-l4)*i/NN;

QQ  = boost::math::quadrature::gauss_kronrod<double, 151>::integrate(ff, l4, lcut[i], 0, 0, &error);

QQ = Q + Q1 + Q2 + QQ;
qq[i] = abs(1. - QQ/II) ;

cout << "{ " << lcut[i] <<","<<  qq[i] << " },"<<endl;

}



}

//**************************************************** END ****************************************************