
#include "IN41.h"
                                                            

using namespace sp_bessel;
using namespace std;                  

typedef complex<double> dcomplex;          // Defines complex number with double entries.

double Cspeed = 137.035999139;
double Pi     = boost::math::constants::pi<double>();           // Parameter for Drude model of the NP and atomic units
double nm     = 100./5.2918;
double wp     = 0.555;
double Gamma  = 0.00555;


//************************************************
// Maximum number of multipole to take into account

const int Lmax = 38;

// Order N of the gaussian quadrature (must be odd), the integral is
// performed using a 2*N +1 Konrod quadrature, the error is |K2n+1 - Gn|
// Omega partition ******************************
const int nw1 = 121;   
const int nw2 = 101;   
const int nw3 = 31;                   // Werner Au
const int nw4 = 15;  
const int NN = nw1 + nw2 + nw3 + nw4;


double w1 = 0.;
double w2 = .8;
double w3 = 2.;
double w4 = 5.;
double w5 = 50.;   // This last frequency strongly depends of impact parameter and speed!!
                   // In Au, b=1.5nm a=1nm, for v=0.5 appx 35,  v=0.8 appx    ,v=.99 appx 180..

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
dcomplex LL(double w, double w0, double G, double A) {

  return A/(w0*w0 - w*w - 1i*w*G );
}


dcomplex eps(double w){
return 1. + LL(w, 0.                , 0.007349864496171 , 0.152742986686889)
          + LL(w, 0.146997289923419 , 0.055123983721282 , 0.060232866544963)
          + LL(w, 0.268270054110239 , 0.12127276418682  , 0.074008096113541)
          + LL(w, 0.47039132775494  , 0.433642005274085 , 0.249709798748062)
          + LL(w, 0.694562194888153 , 2.60920189614068  , 0.983308298910026)
          + LL(w, 0.731311517369008 , 0.106573035194478 , 0.088728684574082)
          + LL(w, 1.0620554196967   , 0.143322357675333 , 0.067525635140093)
          + LL(w, 1.42219878000908  , 0.477741192251111 , 0.100883298899298)
          + LL(w, 2.36298143551895  , 1.90728983675636  , 0.734678910324206);
} 



//**********************************************************************


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
	res += pow(betav,-(l+1))*pow(1i,l-j)*factorial2(2*l+1)*alm(l,m)*III[l-1][m][j]
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




void DP(double r, double vv, double b, double a, dcomplex DPx[6], dcomplex DPz[6]){

double dm;
double dl;

double dpdwx[6], dpdwz[6];

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


// Calls the Gauss - Konrod function
double xi[2*NN + 4], xk[2*NN + 4], xg[2*NN + 4];
Omegas(xi, xk, xg);

double w, k0;




char filenamexee[sizeof "dpdw_a1nm_v0.99c_b1.5nm_Au_exact_L20.dat"];
sprintf(filenamexee, "dpdw_a%.2gnm_v%.2g_b%.2gnm_Au_exact_L%d.dat", a/(1.*nm), vv, b/(1.*nm),Lmax);
FILE *fpx = fopen(filenamexee,"w+");
fprintf(fpx,"Momentum Spectrum, a: %.2gnm    v: %.2gc   b: %.2gnm    Lmax: %d \n", a/(1.*nm), vv, b/(1.*nm),Lmax);
fprintf(fpx,"\n");
fprintf(fpx,"         w(au)                   dpEdwx                 dpHdwx                 dpEdwz                dpHdwz                  dpEsdwx                dpHsdwx                dpEsdwz                dpHsdwz              dpEedwx                dpHedwx                dpEedwz                dpHedwz\n");




dcomplex AA[Lmax][2*Lmax+1], BB[Lmax][2*Lmax+1];
dcomplex aa, bb;

for (int l = 1; l <= Lmax; ++l){
    for (int m = -l; m <= l; ++m){ 
    AA[l-1][l+m] = A(l,m,vv);
    BB[l-1][l+m] = B(l,m,vv);
    }
}


for (int i = 0; i < 6; ++i){
DPx[i] = 0.;
DPz[i] = 0.;
}


for (int i = 0; i < iw4; ++i){     // Here we can select the Omega interval to integrate "2*NN + 4" for total.
    w = xi[i];
    k0 = w/Cspeed;


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




fprintf(fpx,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",w,dpdwx[0],dpdwx[1],dpdwz[0],dpdwz[1], dpdwx[2], dpdwx[3], dpdwz[2], dpdwz[3],dpdwx[4],dpdwx[5],dpdwz[4],dpdwz[5]);


for (int rr = 0; rr < 6; ++rr){ 
  DPx[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwx[rr];
  DPz[rr] += (xk[i]*(1. + 1i) - 1i*xg[i])*dpdwz[rr];
}



}// for w


cout << "b : " << b/nm << " nm" <<endl;
cout << "DPx : " << DPx[0] + DPx[1] + DPx[2] + DPx[3] + DPx[4] + DPx[5] << endl;
cout << "DPz : " << DPz[0] + DPz[1] + DPz[2] + DPz[3] + DPx[4] + DPx[5] << endl;
cout << endl;


} //end void



int main(void){

double v = 0.5;;                    // Defines parameter electron speed vv, frequency a.u,  
double a = 50.*nm;                        // impact parameter b, NP ratio, Integration ratio r.
double r = 50.05*nm;

dcomplex DPx[6], DPz[6];
double Px, Pz;


char filenamex[sizeof "DP_a1nm_v0.5c_vs_b_comp_Au_X.dat"];
sprintf(filenamex, "DP_a%.2gnm_v%gc_vs_b_comp_Au_X.dat", a/nm, v);
FILE *fppx = fopen(filenamex,"a+");
fprintf(fppx,"Total momentum transfered X, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppx,"\n");
fprintf(fppx,"     b        DPx     DPEesX     DPHesX     DPEssX     DPHssX     DPEeeX     DPHeeX\n");



char filenamez[sizeof "DP_a1nm_v0.5c_vs_b_comp_Au_Z.dat"];
sprintf(filenamez, "DP_a%.2gnm_v%gc_vs_b_comp_Au_Z.dat", a/nm, v);
FILE *fppz = fopen(filenamez,"a+");
fprintf(fppz,"Total momentum transfered Z, a: %.2gnm     v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppz,"\n");
fprintf(fppz,"     b        DPz     DPEesZ     DPHesZ     DPEssZ     DPHssZ     DPEeeZ     DPHeeZ\n");



char filenamexe[sizeof "DP_a1nm_v0.5c_vs_b_comp_Au_X_err.dat"];
sprintf(filenamexe, "DP_a%.2gnm_v%gc_vs_b_comp_Au_X_err.dat", a/nm, v);
FILE *fppxe = fopen(filenamexe,"a+");
fprintf(fppxe,"Total momentum transfered X err, a: %.2gnm      v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppxe,"\n");
fprintf(fppxe,"     b        errPx     errPEesX     errPHesX     errPEssX     errPHssX     errPEeeX     errPHeeX\n");


char filenamezerr[sizeof "DP_a1nm_v0.5c_vs_b_comp_Au_Z_err.dat"];
sprintf(filenamezerr, "DP_a%.2gnm_v%gc_vs_b_comp_Au_Z_err.dat", a/nm, v);
FILE *fppze = fopen(filenamezerr,"a+");
fprintf(fppze,"Total momentum transfered Z err, a: %.2gnm      v: %.2gc    Lmax: %d  \n", a/nm, v, Lmax);
fprintf(fppze,"\n");
fprintf(fppze,"     b        errPz     errPEesZ     errPHesZ     errPEssZ     errPHssZ     errPEeeZ     errPHeeZ\n");






double b, bb[] = {56.} ;

int k, rr;


cout <<  "Starting ... "<< endl << endl;



#pragma omp parallel for shared(r, bb, a,v) private(b,k,Px,Pz,DPx,DPz) num_threads(8)
for( k=0; k < 1; ++k){

b = bb[k]*nm;

DP(r, v, b, a, DPx, DPz); 

Px=0.;
Pz=0.;

for ( rr = 0; rr < 6; ++rr){ 
  Px += DPx[rr].real();
  Pz += DPz[rr].real();}

fprintf(fppx ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b/nm,Px,DPx[0].real(),DPx[1].real(), DPx[2].real(),DPx[3].real(),DPx[4].real(),DPx[5].real());
fprintf(fppxe,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",      b/nm,   DPx[0].imag(),DPx[1].imag(), DPx[2].imag(),DPx[3].imag(),DPx[4].imag(),DPx[5].imag());

fprintf(fppz ,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",b/nm,Pz,DPz[0].real(),DPz[1].real(), DPz[2].real(),DPz[3].real(),DPz[4].real(),DPz[5].real());
fprintf(fppze,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g \n",      b/nm,   DPz[0].imag(),DPz[1].imag(), DPz[2].imag(),DPz[3].imag(),DPz[4].imag(),DPz[5].imag());

}

return(0);
}

//**************************************************** END ****************************************************
//**************************************************** END ****************************************************