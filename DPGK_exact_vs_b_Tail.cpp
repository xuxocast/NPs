#include <iostream>      
#include <fstream>   
#include <boost/math/special_functions/bessel.hpp>                 // BOOST LIBRARIES:  1. BesselK in external fields
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>               // Lengendre Plm
#include <boost/math/quadrature/gauss_kronrod.hpp>    
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

// Spherical Bessel (Hankel) funcions with *complex argument* (very important) for the scatterred fields. 
// Only found in SLATEC project (FORTRAN 77).

#include <complex_bessel_bits/sph_besselFunctions.h>  

// Source:       https://github.com/valandil/complex_bessel     



#include "IN51.h"


using namespace std;  


    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;
    using std::sin;
    using std::cos;
    using std::conj;


typedef boost::multiprecision::cpp_complex<50> dcomplex;          // Defines complex number with double entries.
typedef boost::multiprecision::cpp_bin_float_50  dfloat;


//typedef boost::multiprecision::cpp_complex_quad dcomplex;
//typedef boost::multiprecision::cpp_bin_float_quad dfloat;


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

//const int Lmax = 3;

//************************************************
dfloat Cspeed = 137.035999084;
dfloat Pi     = boost::math::constants::pi<dfloat>();           // Parameter for Drude model of the NP and atomic units
dfloat nm     = dfloat(100)/5.2918;
dfloat wp     = 0.555;
dfloat Gamma  = 0.00555;
dcomplex z1{0, 1};



// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|

const int nw1 = 15;  

dfloat w1 = 0;
dfloat w2 = .1;

int iw1 = 2*nw1 + 1;




//**********************************************************************
// Common functions
//********************************************************************** 

//dcomplex eps(dfloat w ){
//return 1 - wp*wp/(w*(w + Gamma*z1 ));
//} 



dcomplex LL(dfloat w, dfloat w0, dfloat G, dfloat A) {
  return A/(w0*w0 - w*w - z1*w*G );
}


dcomplex eps(dfloat w){
return dfloat(1) + LL(w, 0.                , 0.007349864496171 , 0.152742986686889)
          + LL(w, 0.146997289923419 , 0.055123983721282 , 0.060232866544963)
          + LL(w, 0.268270054110239 , 0.12127276418682  , 0.074008096113541)
          + LL(w, 0.47039132775494  , 0.433642005274085 , 0.249709798748062)
          + LL(w, 0.694562194888153 , 2.60920189614068  , 0.983308298910026)
          + LL(w, 0.731311517369008 , 0.106573035194478 , 0.088728684574082)
          + LL(w, 1.0620554196967   , 0.143322357675333 , 0.067525635140093)
          + LL(w, 1.42219878000908  , 0.477741192251111 , 0.100883298899298)
          + LL(w, 2.36298143551895  , 1.90728983675636  , 0.734678910324206);
} 


dfloat BesselK(int n, dfloat x){
    return boost::math::cyl_bessel_k(n,x);
}                                                  



dfloat LP(int l, int m, dfloat x){
  return boost::math::legendre_p(l,m,x);
}

dfloat factorial2(int n){
   return boost::math::double_factorial<dfloat>(n);
}

dfloat factorial(int n){
  return  boost::math::factorial<dfloat>(n);
}



//**************************************************************************************
// Recursive functions and parameter needed for calculate the scatterd fields
//*************************************************************************************

dfloat alm(int l, int m){
  if (abs(m)<=l)
  {
    return pow((2.*l+1.)*factorial(l-m)/(4.*Pi*factorial(l+m)),0.5);
  }
  else{return 0.;}
}


dcomplex  A(int l, int m, dfloat betav){

dfloat gam = 1./sqrt(1.0 - pow(betav,2.0));
dcomplex  res = 0.;

if(m >= 0){
for (int j = m; j <= l; ++j){
  res += pow(betav,-(l+1))*pow(z1,l-j)*factorial2(2*l+1)*alm(l,m)*III[l-1][m][j]
  /(pow(2.*gam,j)*factorial(l-j)*factorial((j-m)/2)*factorial((j+m)/2));
  }
return res;
 } 
else {
  return pow(-1.,abs(m))*A(l,abs(m),betav);
 }
}



dcomplex  B(int l, int m, dfloat betav){
return A(l,m+1,betav)*sqrt((l+m+1)*(l-m)) - A(l,m-1,betav)*sqrt((l-m+1)*(l+m));
}





dcomplex  PsiE(int l, int m, dfloat w, dfloat b,dfloat betav, dcomplex BB){
dfloat gam = dfloat(1)/sqrt(dfloat(1) - betav*betav);
return -dfloat(2)*Pi*pow(z1,1-l)*w*BB*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*gam*l*(l+1));
}

dcomplex PsiM(int l, int m, dfloat w, dfloat b,dfloat betav, dcomplex AA){
dfloat gam = dfloat(1)/sqrt(dfloat(1) - betav*betav);
return -(dfloat(4)*m)*Pi*pow(z1,1-l)*betav*w*AA*BesselK(m,w*b/(betav*Cspeed*gam))/(Cspeed*Cspeed*l*(l+1));
}



complex<double> double_j(int n, complex<double>  z){   
if( z == 0. ){
    if(n==0){ return 1.;}
      else{ return 0.; }
}
 else{return sp_bessel::sph_besselJ(n,z);}
}


complex<double> double_h(int n, complex<double>  z){   
if( z == 0. ){
      if(n==0){ return 1.;}
          else{ return 0.; }}
 else{return sp_bessel::sph_hankelH1(n,z);}
}


dcomplex j(int n, dcomplex  zz){  
   complex<double> z = static_cast<complex<double>>(zz);
   complex<double> jj =  double_j(n, z);
   return static_cast<dcomplex>(jj);
}


dcomplex h(int n, dcomplex  zz){  
   complex<double> z = static_cast<complex<double>>(zz);
   complex<double> hh =  double_h(n, z);
   return static_cast<dcomplex>(hh);
}


dcomplex djj(int n, dcomplex x){
  return j(n,x) + 0.5*x*(j(n-1,x)-j(n+1,x) - j(n,x)/x);
}

dcomplex dhh(int n, dcomplex x){
  return h(n,x) + 0.5*x*(h(n-1,x)-h(n+1,x) - h(n,x)/x);
}


// *********** Polarizabilities *********************


dcomplex tE(int l, dfloat w, dfloat a){
dfloat x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -z1*( eps(w)*j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - eps(w)*j(l,xi)*dhh(l,x0) ); 
}

dcomplex tM(int l, dfloat w, dfloat a){
dfloat x0 = w*a/Cspeed;
dcomplex xi = x0*sqrt(eps(w));
return -z1*( j(l, xi)*djj(l,x0) - j(l,x0)*djj(l,xi) )/( h(l,x0)*djj(l,xi) - j(l,xi)*dhh(l,x0) );
}



dcomplex fs(int l, dcomplex z){
dfloat dl = dfloat(1)*l;
  return (dl + dfloat(1))*z1*h(l,z)/z  - z1*h(l+1,z);
}

dcomplex fe(int l, dcomplex z){
dfloat dl = dfloat(1)*l;
  return (dl + dfloat(1))*j(l,z)/z  - j(l+1,z);
}












//********************************************************************************
// ****************** Gauss Konrod Quadrature weights and abscisas ********************
//********************************************************************************

//***********************************************************************
void gauss_konrod_w(double GKQt[2*nw1 + 1][3]){

  int const N1 = nw1;

   auto XGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::abscissa();    
   auto WGKt = boost::math::quadrature::gauss_kronrod<double, 2*N1+1>::weights();      
   auto WGt =  boost::math::quadrature::gauss<double, N1>::weights();

   double WGLt[2*N1+1];

// Changes the n order array of Gaussian quadrature to a 2n+1 array
for(int i = 0; i < XGKt.size(); ++i){
      if (i % 2 == 0){
           WGLt[i] = WGt[i/2];}
      else WGLt[i] = 0.;
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
}



void Omegas(dfloat xi[2*nw1 + 1], dfloat xk[2*nw1 + 1], dfloat xg[2*nw1 + 1]){

double GKQw1[2*nw1 + 1][3];
gauss_konrod_w(GKQw1);

double w;
int NN1 = 2*nw1 + 1;

for(int i=0; i < NN1;  ++i){
xi[i] = ((w2-w1)/2.)*GKQw1[i][0] + ((w2+w1)/2.);;
xk[i] = ((w2-w1)/2.)*GKQw1[i][1];
xg[i] = ((w2-w1)/2.)*GKQw1[i][2];
 }
}








void DP(int Lmax, dfloat r, dfloat vv, dfloat b, dfloat a, dcomplex DPx[6], dcomplex DPz[6]){

dfloat dm;
dfloat dl;

dfloat dpdwx[6], dpdwz[6];

dcomplex CM[Lmax][2*Lmax +1], tMl[Lmax];
dcomplex DE[Lmax][2*Lmax +1], tEl[Lmax];
dcomplex dz[2][Lmax], df[2][Lmax];

dcomplex IErrx[4], IErrz[4], IHrrx[4], IHrrz[4];
dcomplex IErtx[4], IErtz[4], IHrtx[4], IHrtz[4];
dcomplex IErf[4], IHrf[4];
dcomplex IEttx[4], IEttz[4], IHttx[4], IHttz[4];
dcomplex IEffx[4], IEffz[4], IHffx[4], IHffz[4];

dfloat dl1, dl2, dm1, dm2;

double IN1, IV1, IW1, IW2, IW3, IU1, IU2, IU3, IU4;
double IM1, IZ1, IX1, IX2, IX3, IY1, IY2, IY3, IY4, ID1, ID2;
FUN(IM,IN,IU,IV,IW,IX,IY,IZ,ID,III);


dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;


// Calls the Gauss - Konrod function
dfloat xi[2*nw1 + 1], xk[2*nw1 + 1], xg[2*nw1 + 1];
Omegas(xi, xk, xg);

dfloat w, k0;




for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}


for (int i = 0; i < 6; ++i){
DPx[i] = dfloat(0);
DPz[i] = dfloat(0);
}




for (int i = 0; i < iw1; ++i){     // Here we can select the Omega interval to integrate "2*NN + 4" for total.
    w = dfloat(1)*xi[i];
    k0 = w/Cspeed;


// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = z1*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(z1,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(z1,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dpdwx[rr] = dfloat(0);
dpdwz[rr] = dfloat(0);
}




// Calculates  common parameters of the scatterred fields
for (int l = 1; l <= Lmax; l++){
    tEl[l-1] = tE(l,w,a);
    tMl[l-1] = tM(l,w,a);
    
    dz[0][l-1] = z1*h(l,k0*r);
    df[0][l-1] = fs(l,k0*r);
    dz[1][l-1] = j(l,k0*r);
    df[1][l-1] = fe(l,k0*r);

    for (int m = -l; m <= l; ++m){

        aa = AA[l-1][m+l];
        bb = BB[l-1][l+m];
    
        DE[l-1][m + l] = pow(z1,l)*alm(l,m)*PsiE(l,m,w,b,vv, bb);
        CM[l-1][m + l] = pow(z1,l)*alm(l,m)*PsiM(l,m,w,b,vv, aa);
    
    }
} 



for (int rr = 0; rr < 6; ++rr){
dpdwx[rr] = dfloat(0);
dpdwz[rr] = dfloat(0);
}



for (int l1 = 1; l1 <= Lmax; ++l1){
     dl1 = dfloat(1)*l1; 
for (int l2 = 1; l2 <= Lmax; ++l2){
         dl2 = dfloat(1)*l2; 
         for (int m1 = -l1; m1 <= l1; ++m1){
              dm1 = dfloat(1)*m1;
              for (int m2 = -l2; m2 <= l2; ++m2){
                   dm2 = dfloat(1)*m2;

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
                     IErrx[rr] = dfloat(0);
                     IHrrx[rr] = dfloat(0);
                     IErtx[rr] = dfloat(0);
                     IHrtx[rr] = dfloat(0);
                     IEffx[rr] = dfloat(0);
                     IHffx[rr] = dfloat(0);
                     IErf[rr] = dfloat(0);
                     IHrf[rr] = dfloat(0);
                     IEttx[rr] = dfloat(0);
                     IHttx[rr] = dfloat(0);
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
                     IErrz[rr] = dfloat(0);
                     IHrrz[rr] = dfloat(0);
                     IErtz[rr] = dfloat(0);
                     IHrtz[rr] = dfloat(0);
                     IEttz[rr] = dfloat(0);
                     IHttz[rr] = dfloat(0);
                     IEffz[rr] = dfloat(0);
                     IHffz[rr] = dfloat(0);
                    }
                   }




                     
                    dpdwx[0] += (float(1)/(4*Pi))*((0.5*(IErrx[2] - IEttx[2] - IEffx[2]) + IErtx[2] - IErf[2]).real()
                                            + (0.5*(IErrx[3] - IEttx[3] - IEffx[3]) + IErtx[3] - IErf[3])).real();
                    dpdwz[0] += (float(1)/(2*Pi))*((0.5*(IErrz[2] - IEttz[2] - IEffz[2]) - IErtz[2]).real()
                                            + (0.5*(IErrz[3] - IEttz[3] - IEffz[3]) - IErtz[3])).real();

                    dpdwx[1] +=  (float(1)/(4*Pi))*((0.5*(IHrrx[2] - IHttx[2] - IHffx[2]) + IHrtx[2] - IHrf[2]).real()
                                             + (0.5*(IHrrx[3] - IHttx[3] - IHffx[3]) + IHrtx[3] - IHrf[3])).real();
                    dpdwz[1] +=  (float(1)/(2*Pi))*((0.5*(IHrrz[2] - IHttz[2] - IHffz[2]) - IHrtz[2]).real()
                                             + (0.5*(IHrrz[3] - IHttz[3] - IHffz[3]) - IHrtz[3])).real();

                    dpdwx[2] +=  (float(1)/(4*Pi))*((0.5*(IErrx[0] - IEttx[0] - IEffx[0]) + IErtx[0] - IErf[0])).real();
                    dpdwz[2] +=  (float(1)/(2*Pi))*((0.5*(IErrz[0] - IEttz[0] - IEffz[0]) - IErtz[0])).real();

                    dpdwx[3] +=  (float(1)/(4*Pi))*((0.5*(IHrrx[0] - IHttx[0] - IHffx[0]) + IHrtx[0] - IHrf[0])).real();
                    dpdwz[3] +=  (float(1)/(2*Pi))*((0.5*(IHrrz[0] - IHttz[0] - IHffz[0]) - IHrtz[0])).real();

                    dpdwx[4] +=  (float(1)/(4*Pi))*((0.5*(IErrx[1] - IEttx[1] - IEffx[1]) + IErtx[1] - IErf[1])).real();
                    dpdwz[4] +=  (float(1)/(2*Pi))*((0.5*(IErrz[1] - IEttz[1] - IEffz[1]) - IErtz[1])).real();

                    dpdwx[5] +=  (float(1)/(4*Pi))*((0.5*(IHrrx[1] - IHttx[1] - IHffx[1]) + IHrtx[1] - IHrf[1])).real();
                    dpdwz[5] +=  (float(1)/(2*Pi))*((0.5*(IHrrz[1] - IHttz[1] - IHffz[1]) - IHrtz[1])).real();





            } // for m2
        }  // for m1

} // for l2
} //for l1




for (int rr = 0; rr < 6; ++rr){ 
  DPx[rr] += (xk[i]*(dfloat(1) + z1) - z1*xg[i])*dpdwx[rr];
  DPz[rr] += (xk[i]*(dfloat(1) + z1) - z1*xg[i])*dpdwz[rr];
}



}// for w

std::cout << std::setprecision(std::numeric_limits<typename dcomplex::value_type>::digits10);

cout << "b : " << b/nm << endl;
cout << "DPx : " << DPx[0] + DPx[1] + DPx[2] + DPx[3] << endl;
cout << "DPz : " << DPz[0] + DPz[1] + DPz[2] + DPz[3] << endl;
cout << endl;


} //end void





int main(void){


double nm1 = 100./5.2918;



double v = 0.5;                    // Defines parameter electron speed vv, frequency a.u,  
double a = 50;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 50.05;

dcomplex DPx[6], DPz[6];
dfloat Px, Pz;


int Lmax=50;




char filenamex[sizeof "DP_a1nm_v0.5c_vs_b_Tail_Au_X.dat"];
sprintf(filenamex, "DP_a%.2gnm_v%gc_vs_b_Tail_Au_X.dat", a/nm, v);
FILE *fppx = fopen(filenamex,"a+");
fprintf(fppx,"Total momentum transfered X, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppx,"\n");
fprintf(fppx,"     b        DPx     DPEesX     DPHesX     DPEssX     DPHssX     DPEeeX     DPHeeX\n");

char filenamez[sizeof "DP_a1nm_v0.5c_vs_b_Tail_Au_Z.dat"];
sprintf(filenamez, "DP_a%.2gnm_v%gc_vs_b_Tail_Au_Z.dat", a/nm, v);
FILE *fppz = fopen(filenamez,"a+");
fprintf(fppz,"Total momentum transfered Z, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppz,"\n");
fprintf(fppz,"     b        DPz     DPEesZ     DPHesZ     DPEssZ     DPHssZ     DPEeeZ     DPHeeZ\n");

char filenameex[sizeof "DP_a1nm_v0.5c_vs_b_Tail_Au_X_err.dat"];
sprintf(filenameex, "DP_a%.2gnm_v%gc_vs_b_Tail_Au_X_err.dat", a/nm, v);
FILE *fppex = fopen(filenameex,"a+");
fprintf(fppex,"Total momentum transfered err X, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppex,"\n");
fprintf(fppex,"     b        DPx     DPEesX     DPHesX     DPEssX     DPHssX     DPEeeX     DPHeeX\n");

char filenameez[sizeof "DP_a1nm_v0.5c_vs_b_Tail_Au_Z_err.dat"];
sprintf(filenameez, "DP_a%.2gnm_v%gc_vs_b_Tail_Au_Z_err.dat", a/nm, v);
FILE *fppez = fopen(filenameez,"a+");
fprintf(fppez,"Total momentum transfered err Z, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppez,"\n");
fprintf(fppez,"     b        DPz     DPEesZ     DPHesZ     DPEssZ     DPHssZ     DPEeeZ     DPHeeZ\n");




int k, rr;
double ddPx, ddPz, dpx[4],dpz[4],idpx[4],idpz[4];



double b, bb[] = {50.3, 50.5,51.,52.,53.,54.,55.,56.} ;


cout <<  "Starting ... "<< endl << endl;

#pragma omp parallel for shared(r, bb, a,v) private(b,k,rr,Px,Pz,DPx,DPz) num_threads(8)
for( k=0; k < 8; ++k){

b = bb[k];

DP(Lmax, r*nm, dfloat(1)*v, b*nm, a*nm, DPx, DPz); 

Px=0.;
Pz=0.;

for (rr = 0; rr < 6; ++rr){ 
  Px += DPx[rr].real();
  Pz += DPz[rr].real();

  dpx[rr] = static_cast<double>(DPx[rr].real());
  dpz[rr] = static_cast<double>(DPz[rr].real());

  idpx[rr] = static_cast<double>(DPx[rr].imag());
  idpz[rr] = static_cast<double>(DPz[rr].imag());

}


ddPx = static_cast<double>(Px);
ddPz = static_cast<double>(Pz);




fprintf(fppx ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b,ddPx,    dpx[0], dpx[1], dpx[2], dpx[3], dpx[4], dpx[5]);
fprintf(fppex ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b,     idpx[0],idpx[1],idpx[2],idpx[3],idpx[4],idpx[5]);

fprintf(fppz ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b,ddPz,dpz[0], dpz[1], dpz[2], dpz[3], dpz[4],dpz[5]);
fprintf(fppez ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b,idpz[0], idpz[1],idpz[2],idpz[3],idpz[4],idpz[5]);


}
            



return(0);
}
//**************************************************** END ****************************************************