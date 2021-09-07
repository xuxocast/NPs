//*****************************************************************************************************************
//
//   DPGK_SingleL_Drude.cpp:
//
//    This program calculates the spectral linear momentum dp/dw transferred to a 
//   NP due a fast swfit electron. It calculate the closed surface integral of
//   the Maxwell stress tensor via Gauss Konrod quadrature, the real part of               
//   the the dpdw void function represents the momentum, and the imaginary part 
//   the estimated error. The spectrum is written to a file named "dpdw*.dat".
//   Also, the integral in frequency space is calculated via Gauss - Kronrod using
//   a partition set to fast convergence. The results are written to a file named "DP*.dat". 
//
//
//          By Jesús Castrejon, jcastrejon@ciencias.unam.mx (27/10/2018)
//
//*****************************************************************************************************************
//*****************************************************************************************************************

#include "IN11.h"

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel                                                                          

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.


//double wp     = 0.555;               // Drude Al
//double Gamma  = 0.00555;

double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;
double wp     = 0.555;
double Gamma  = 0.00555;


//************************************************
// Maximum number of multipole to take into account

const int Lmax = 8;

// Order N of the gaussian (Legendre) quadrature (must be even)
//const int Nt = 120;       
//const int Nf = 130;  

const int Nt = 25;      
const int Nf = 30;        
//*******************************************





const int nw1 = 51;   
const int nw2 = 201;   
const int nw3 = 51;                   // Drude Al
const int nw4 = 51;                   //for .99 81 is 7 order error, .9 61 is good, .5 41 is good, .1 21 is good, .3 31 is very good
const int NN = nw1 + nw2 + nw3 + nw4;

double w1 = 0.;
double w2 = .3;
double w3 = .4;
double w4 = 2.;
double w5 = 60.;

int iw1 = 2*nw1 + 1;
int iw2 = 2*(nw1 + nw2) + 2;
int iw3 = 2*(nw1 + nw2 +nw3) + 3;
int iw4 = 2*NN + 4;


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



//********************************************************************************
// ****************** Gauss Konrod Quadrature weights and abscisas ********************
//********************************************************************************
void gauss_legendre(double GLQt[Nt][2], double GLQf[Nf][2]){
   auto XGt = boost::math::quadrature::gauss<double, Nt>::abscissa();       
   auto WGt =  boost::math::quadrature::gauss<double, Nt>::weights();

   auto XGf = boost::math::quadrature::gauss<double, Nf>::abscissa();
   auto WGf =  boost::math::quadrature::gauss<double, Nf>::weights();


for(int i = 0; i < Nt; ++i){
         if (i < Nt/2){ 
         GLQt[i][0] = -XGt[-i + Nt/2 - 1];                
         GLQt[i][1] = WGt[-i + Nt/2 - 1];
          }
      else{
         GLQt[i][0] = XGt[i - Nt/2];                
         GLQt[i][1] = WGt[i - Nt/2];   
      }
 }

if (Nf % 2 == 0){
for(int i = 0; i < Nf; ++i){
         if (i < Nf/2){ 
         GLQf[i][0] = -XGf[-i + Nf/2 - 1];                
         GLQf[i][1] = WGf[-i + Nf/2 - 1];
          }
      else{
         GLQf[i][0] = XGf[i - Nf/2];                
         GLQf[i][1] = WGf[i - Nf/2];   
      }
 }
}

}
//***********************************************************************
void gauss_konrod_w(double GKQt[2*nw1 + 1][3], double GKQf[2*nw2 + 1][3]){

  int const N1 = nw1;
  int const N2 = nw2;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   auto XGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::abscissa();
   auto WGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::weights();
   auto WGf =  boost::math::quadrature::gauss<double, N2>::weights();

   double WGLt[2*N1+1], WGLf[2*N2+1];


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
      if (i <= N1){ 
         GKQt[i][0] = -1.*XGKt[N1 - i];
         GKQt[i][1] = WGKt[N1 - i];                
         GKQt[i][2] = WGLt[N1 - i];
      }
      else{
         GKQt[i][0] = XGKt[-N1 + i];
         GKQt[i][1] = WGKt[-N1 + i];
         GKQt[i][2] = WGLt[-N1 + i];
      }
}

for(int i = 0; i < 2*XGKf.size() - 1; ++i){
      if (i <= N2){ 
         GKQf[i][0] = -1.*XGKf[N2 - i];
         GKQf[i][1] = WGKf[N2 - i];
         GKQf[i][2] = WGLf[N2 - i];
      }
      else{
         GKQf[i][0] = XGKf[-N2 + i];
         GKQf[i][1] = WGKf[-N2 + i];
         GKQf[i][2] = WGLf[-N2 + i];
      }
}

}





void gauss_konrod_w2(double GKQt[2*nw3 + 1][3], double GKQf[2*nw4 + 1][3]){

  int const N1 = nw3;
  int const N2 = nw4;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   auto XGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::abscissa();
   auto WGKf = boost::math::quadrature::gauss_kronrod<double, 2*N2+1>::weights();
   auto WGf =  boost::math::quadrature::gauss<double, N2>::weights();

   double WGLt[2*N1+1], WGLf[2*N2+1];


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
      if (i <= N1){ 
         GKQt[i][0] = -1.*XGKt[N1 - i];
         GKQt[i][1] = WGKt[N1 - i];                
         GKQt[i][2] = WGLt[N1 - i];
      }
      else{
         GKQt[i][0] = XGKt[-N1 + i];
         GKQt[i][1] = WGKt[-N1 + i];
         GKQt[i][2] = WGLt[-N1 + i];
      }
}

for(int i = 0; i < 2*XGKf.size() - 1; ++i){
      if (i <= N2){ 
         GKQf[i][0] = -1.*XGKf[N2 - i];
         GKQf[i][1] = WGKf[N2 - i];
         GKQf[i][2] = WGLf[N2 - i];
      }
      else{
         GKQf[i][0] = XGKf[-N2 + i];
         GKQf[i][1] = WGKf[-N2 + i];
         GKQf[i][2] = WGLf[-N2 + i];
      }
}

}




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
//******************************************************************************************************
//******************************************************************************************************







void Omegas(double xi[2*NN + 4], double xk[2*NN + 4], double xg[2*NN + 4]){


double GKQw1[2*nw1 + 1][3], GKQw2[2*nw2 + 1][3];
gauss_konrod_w(GKQw1,GKQw2);
double GKQw3[2*nw3 + 1][3], GKQw4[2*nw4 + 1][3];
gauss_konrod_w2(GKQw3,GKQw4);

double w;

int NN1 = 2*nw1 + 1;
int NN2 = 2*nw2 + 1;
int NN3 = 2*nw3 + 1;
int NN4 = 2*nw4 + 1;

for(int i=0; i < NN1;  ++i){
xi[i] = ((w2-w1)/2.)*GKQw1[i][0] + ((w2+w1)/2.);;
xk[i] = ((w2-w1)/2.)*GKQw1[i][1];
xg[i] = ((w2-w1)/2.)*GKQw1[i][2];
}

for(int i=0; i < NN2;  ++i){
xi[i + NN1] = ((w3-w2)/2.)*GKQw2[i][0] + ((w3+w2)/2.);;
xk[i + NN1] = ((w3-w2)/2.)*GKQw2[i][1];
xg[i + NN1] = ((w3-w2)/2.)*GKQw2[i][2];
}

for(int i=0; i < NN3;  ++i){
xi[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][0] + ((w4+w3)/2.);;
xk[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][1];
xg[i + NN2+ NN1] = ((w4-w3)/2.)*GKQw3[i][2];
}

for(int i=0; i < NN4;  ++i){
xi[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][0] + ((w5+w4)/2.);;
xk[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][1];
xg[i + NN3+ NN1+ NN2] = ((w5-w4)/2.)*GKQw4[i][2];
}

}




void DP(double r, double vv, double b, double a){

double dm;
double dl;

double Errx[4], Errz[4];


dcomplex dpdwx[4], dpdwz[4];
dcomplex DPx[4], DPz[4];

dcomplex clm, dle;
dcomplex CM[Lmax][2*Lmax +1];
dcomplex DE[Lmax][2*Lmax +1];
dcomplex dhl[Lmax], dfs[Lmax];


dcomplex Ers, Ets, Efs;
dcomplex Hrs, Hts, Hfs;

dcomplex Eex, Eey, Eez;
dcomplex Hey, Hex;

double TEesx, TEesz, THesx, THesz;
double theta, phi;

// Calls the Gauss - Konrod function
double GKQf[Nf][2], GKQt[Nt][2];
gauss_legendre(GKQt,GKQf);

double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];
Omegas(xi, xk, xg);

double w, k0;


dcomplex IErrx, IErrz, IHrrx, IHrrz;
dcomplex IErtx, IErtz, IHrtx, IHrtz;
dcomplex IErf, IHrf;
dcomplex IEttx, IEttz, IHttx, IHttz;
dcomplex IEffx, IEffz, IHffx, IHffz;

double Tex, Thx, Tez, Thz;

double dl1, dl2, dm1, dm2;

double IN1, IV1, IW1, IW2, IW3, IU1, IU2, IU3, IU4;
double IM1, IZ1, IX1, IX2, IX3, IY1, IY2, IY3, IY4, ID1, ID2;


FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


double dPHsz, dPHsx, dPEsx, dPEsz;



char filename[sizeof "DP_a1nm_v0.99c_b1.5nm_Lmax_10.dat"];
sprintf(filename, "DP_a%.2gnm_v%.2g_b%.2gnm_Lmax_%d.dat", a/(1.*nm), vv, b/(1.*nm),Lmax);
FILE *fpp = fopen(filename,"w+");
fprintf(fpp,"Total momentum transfered, a: %.2gnm    v: %.2gc   b: %.2gnm   \n", a/(1.*nm), vv, b/(1.*nm));
fprintf(fpp,"\n");
fprintf(fpp,"           DPEx                    DPHx                  DPEz                  DPHz            DPEsx                   DPHsx                 DPEsz                 DPHsz\n");


char filenamer[sizeof "DP_a1nm_v0.99c_b1.5nm_error_Lmax_10.dat"];
sprintf(filenamer, "DP_a%.2gnm_v%.2g_b%.2gnm_error_Lmax_%d.dat", a/(1.*nm), vv, b/(1.*nm),Lmax);
FILE *fppe = fopen(filenamer,"w+");
fprintf(fppe,"Total momentum transfered, a: %.2gnm    v: %.2gc   b: %.2gnm   \n", a/(1.*nm), vv, b/(1.*nm));
fprintf(fppe,"\n");
fprintf(fppe,"      errDPEx                   errDPHx                 errDPEz                errDPHz        errDPEsx               errDPHsx               errDPEsz              errDPHsz\n");


char filenamex[sizeof "dpdw_a1nm_v0.99c_b1.5nm_Lmax_%d.dat"];
sprintf(filenamex, "dpdw_a%.2gnm_v%.2g_b%.2gnm_Lmax_%d.dat", a/(1.*nm), vv, b/(1.*nm),Lmax);
FILE *fpx = fopen(filenamex,"w+");
fprintf(fpx,"          w(au)                   dpEdwx                 dpHdwx                 dpEdwz                dpHdwz                  dpEsdwx                dpHsdwx                dpEsdwz                dpHsdwz\n");





dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}

dcomplex aa, bb;





DPx[0] = 0.;
DPz[0] = 0.;
DPx[1] = 0.;
DPz[1] = 0.;

DPx[2] = 0.;
DPz[2] = 0.;
DPx[3] = 0.;
DPz[3] = 0.;



for (int i = 0; i < iw4; ++i){     // Here we can select the Omega interval to integrate "2*NN + 4" for total.
    w = xi[i];
    k0 = w/Cspeed;


// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    dle = pow(1i,l)*tE(l,w,a);
    clm = pow(1i,l)*tM(l,w,a);
    dhl[l-1] = h(l,k0*r);
    dfs[l-1] = fs(l,k0*r);


    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = alm(l,m)*dle*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = alm(l,m)*clm*PsiM(l,m,w,b,vv, aa);
    
    }
} 


dpdwx[0] = 0.;
dpdwz[0] = 0.;
dpdwx[1] = 0.;
dpdwz[1] = 0.;

dpdwx[2] = 0.;
dpdwz[2] = 0.;
dpdwx[3] = 0.;
dpdwz[3] = 0.;



for (int l1 = 1; l1 <= Lmax; ++l1){
     dl1 = 1.*l1; 

for(int j=0;  j < Nt ;  ++j) {
      theta = Pi*(GKQt[j][0] + 1.0)/2.0;

    for(int k=0; k < Nf ;  ++k){
       phi= Pi*(GKQf[k][0] + 1.0);

Ers = 0.;
Ets = 0.;
Efs = 0.;

Hrs = 0.;
Hts = 0.;
Hfs = 0.;

// ************** Scatterrd fields ****************

        dl = 1.*l1;
        int l = l1;
        for (int m = -l; m <= l; ++m){
             dm = 1.*m;
          
             Ers += exp(1i*dm*phi)*DE[l-1][m+l]*dl*(dl+1.)*LP(l,m,cos(theta))*(1i*dhl[l-1])/(k0*r);

             Ets += -exp(1i*dm*phi)*(CM[l-1][m+l]*LP(l,m,cos(theta))*(1i*dm*dhl[l-1])/sin(theta) + 
                      DE[l-1][m+l]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*dfs[l-1]);

             Efs += 1i*exp(1i*dm*phi)*(DE[l-1][m+l]*LP(l,m,cos(theta))*(dm*dfs[l-1])/sin(theta) + 
                      CM[l-1][m+l]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*1i*dhl[l-1]);



             Hrs += exp(1i*dm*phi)*CM[l-1][m+l]*dl*(dl+1.)*LP(l,m,cos(theta))*(1i*dhl[l-1])/(k0*r);

             Hts += exp(1i*dm*phi)*(DE[l-1][m+l]*LP(l,m,cos(theta))*(1i*dm*dhl[l-1])/sin(theta) - 
                      CM[l-1][m+l]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*dfs[l-1]);

             Hfs += 1i*exp(1i*dm*phi)*(CM[l-1][m+l]*LP(l,m,cos(theta))*(dm*dfs[l-1])/sin(theta) - 
                      DE[l-1][m+l]*( (dl+1.)*LP(l,m,cos(theta))/tan(theta) - (dl - dm +1.)*LP(l+1,m,cos(theta))/sin(theta))*1i*dhl[l-1]);
        }



// ************ External fields ******************
Eex = Ex(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Eey = Ey(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Eez = Ez(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));

Hex = Hx(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
Hey = Hy(w,vv,b,r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));


// ***************+ Maxwell stress tensor *********************
TEesx = (Eex*conj(Ers) - conj(Efs)*(Eey*sin(theta) +  Eez*cos(theta)*sin(phi)) + Eez*conj(Ets)*cos(phi)).real();
TEesz = (Eez*conj(Ers) - conj(Ets)*(Eex*cos(phi) + Eey*sin(phi)) + conj(Efs)*(Eex*cos(theta)*sin(phi) - Eey*cos(theta)*cos(phi))).real();
THesx = (Hex*conj(Hrs) - conj(Hfs)*Hey*sin(theta)).real();
THesz = (- conj(Hts)*(Hex*cos(phi) + Hey*sin(phi)) + conj(Hfs)*(Hex*cos(theta)*sin(phi) - Hey*cos(theta)*cos(phi))).real();


    // ************************* SURFACE CLOSED INTEGRAL ***************************   

        dpdwx[0] += TEesx*(GKQf[k][1]*GKQt[j][1])*sin(theta)/8.0;
        dpdwz[0] += TEesz*(GKQf[k][1]*GKQt[j][1])*sin(theta)/8.0;
        dpdwx[1] += THesx*(GKQf[k][1]*GKQt[j][1])*sin(theta)/8.0;  
        dpdwz[1] += THesz*(GKQf[k][1]*GKQt[j][1])*sin(theta)/8.0;

     }
 }






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


                    IErrx = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G

                    IHrrx = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IN1; //G


                                   // Zenital - Radial

                    IErtx = -(  DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1 - (dl1-dm1+1.)*IW2)
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);
     
                    IHrtx =  (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IV1- (dl1-dm1+1.)*IW2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IW1);



                                    // Azimutal - Radial


                    IErf =  (dm2-dm1)*( CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                      + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1); 

                    IHrf =  (dm2-dm1)*( -DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                                       +CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*IU1);    


                                        // Azimutal - Azimutal

                    IEffx =     DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); //G

                    IHffx =     CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);   //   




                                        //  Zenital - Zenital

                    IEttx =     CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*dm1*dm2*IU1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4); 

                    IHttx =     DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*dm1*dm2*IU1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm1*((dl2+1.)*IW1 - (dl2-dm2+1.)*IU3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm2*((dl1+1.)*IW1 - (dl1-dm1+1.)*IU2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*((dl1+1.)*(dl2+1.)*IV1 - (dl2+1.)*(dl1-dm1+1.)*IW2
                              - (dl1+1.)*(dl2-dm2+1.)*IW3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IU4);                                      


               }
                else{IErrx = 0.;
                     IHrrx = 0.;
                     IErtx = 0.;
                     IHrtx = 0.;
                     IEffx = 0.;
                     IHffx = 0.;
                     IErf = 0.;
                     IHrf = 0.;
                     IEttx = 0.;
                     IHttx = 0.;}




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



                    IErrz = DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; // G of "Good" cheked with TestGKL.cpp

                    IHrrz = CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*k0*r*r))*dl1*dl2*(dl1+1.)*(dl2+1.)*IM1; //G


                                              // Zenital - Radial

                    IErtz = -( DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1);  //G
   
                    IHrtz = (- CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dfs[l1-1]/(k0*r))*dl2*(dl2+1.)*((dl1+1.)*IM1 - (dl1-dm1+1.)*ID2)
                             + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*(dhl[l1-1]/(k0*r))*dm1*dl2*(dl2+1.)*ID1); //G

                                           // Azimutal - Azimutal

                    IEffz =     DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                    IHffz =     CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  +  (dl2-dm2+1.)*(dl1-dm1+1.)*IY4);  //G

                                                                      //  Zenital - Zenital

                    IEttz =     CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*dm1*dm2*IY1 
                              + CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              + DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); 

                    IHttz =     DE[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*dhl[l2-1])*dhl[l1-1]*dm1*dm2*IY1 
                              - DE[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*1i*dhl[l1-1]*dm1*((dl2+1.)*IX1 - (dl2-dm2+1.)*IY3)
                              - CM[l1-1][m1+l1]*conj(DE[l2-1][m2+l2]*1i*dhl[l2-1])*dfs[l1-1]*dm2*((dl1+1.)*IX1 - (dl1-dm1+1.)*IY2)
                              + CM[l1-1][m1+l1]*conj(CM[l2-1][m2+l2]*dfs[l2-1])*dfs[l1-1]*((dl1+1.)*(dl2+1.)*IZ1 - (dl2+1.)*(dl1-dm1+1.)*IX2
                              - (dl1+1.)*(dl2-dm2+1.)*IX3  + (dl2-dm2+1.)*(dl1-dm1+1.)*IY4); }  

                                              
                else{IErrz = 0.;
                     IHrrz = 0.;
                     IErtz = 0.;
                     IHrtz = 0.;
                     IEttz = 0.;
                     IHttz = 0.;
                     IEffz = 0.;
                     IHffz = 0.;}

                     Tex = (0.5*(IErrx - IEttx - IEffx) + IErtx - IErf).real();
                     Thx = (0.5*(IHrrx - IHttx - IHffx) + IHrtx - IHrf).real();

                     Tez = (0.5*(IErrz - IEttz - IEffz) - IErtz).real();
                     Thz = (0.5*(IHrrz - IHttz - IHffz) - IHrtz).real();


                 dpdwx[2] += (1./(4.*Pi))*Tex;
                 dpdwx[3] += (1./(4.*Pi))*Thx;
                 dpdwz[2] += (1./(2.*Pi))*Tez;      // Sum l2,m2,m1
                 dpdwz[3] += (1./(2.*Pi))*Thz;



            } // for m2
        }  // for m1

} // for l2
} //for l1




for (int rr = 0; rr < 4; ++rr){
  DPx[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwx[rr].real();
  DPz[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwz[rr].real();
}




// Here print the dpdw's, for each, w, l
fprintf(fpx,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",w,dpdwx[0].real(),dpdwx[1].real(),dpdwz[0].real(),dpdwz[1].real(), dpdwx[2].real(), dpdwx[3].real(), dpdwz[2].real(), dpdwz[3].real());

cout << "In   " << i + 1 << "  of   " << 2*NN + 4 << endl;

}// for w




// Here print the total momentum
fprintf(fpp,"%.17g %.17g %.17g %.17g  %.17g %.17g %.17g %.17g \n",DPx[0].real(),DPx[1].real(),DPz[0].real(),DPz[1].real(), DPx[2].real(),DPx[3].real(),DPz[2].real(),DPz[3].real());
fprintf(fppe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",DPx[0].imag(),DPx[1].imag(),DPz[0].imag(),DPz[1].imag(), DPx[2].imag(),DPx[3].imag(),DPz[2].imag(),DPz[3].imag());


cout << endl;
cout << "DPEx : " << DPx[0] << endl;
cout << "DPHx : " << DPx[1] << endl;                // prints result
cout << endl;
cout << "DPEz : " << DPz[0] << endl;
cout << "DPHz : " << DPz[1] << endl;
cout << endl;
cout << "DPEsx : " << DPx[2] << endl;
cout << "DPHsx : " << DPx[3] << endl;                // prints result
cout << endl;
cout << "DPEsz : " << DPz[2] << endl;
cout << "DPHsz : " << DPz[3] << endl;
cout << endl;
cout << "DPx : " << DPx[0] + DPx[1] + DPx[2] + DPx[3] << endl;
cout << "DPz : " << DPz[0] + DPz[1] + DPz[2] + DPz[3] << endl;
cout << endl;


} //end void











//*********************************************************************************************************
//********************************************** MAIN PROGRAM *********************************************
int main(void){

double vv = 0.5;
double b = 1.5*nm;                        // Defines parameter electron speed vv, frequency a.u,  
double a = 1.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 1.05*nm;



cout.precision(17);
cout << endl;
cout << "Momentum Gauss - Kronrod :" << endl;
cout << endl;
cout << "Lmax = " << Lmax << endl;
cout << endl;
cout << "a = " << a/nm << "nm." << endl;
cout << "b = " << b/nm << "nm." << endl;
cout << "v = " << vv << "c." << endl;
cout << endl;
DP(r, vv, b, a); 

 
return(0);


}

//**************************************************** END ****************************************************