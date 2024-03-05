#include <iostream> 
#include <stdlib.h>
#include <math.h>
#include "../lib/pdf.h"
#include "../lib/jet.h"
#include "../lib/utils.h"
#include "../lib/vegas.h"
#include "../lib/vector4.h"
#include "../lib/interface.h"
#include "../lib/histogram.h"

using namespace std;

model MODEL = kT;
unsigned int PDF = KMR;

double ml = 0.0005;
double S = sqr(1800.0);
double scale_factor = 1.0;

double Mmin = 11.0;
double Mmax = 1500.0;
double qtmin = 0.0;
double qtmax = 200.0;

const char* prefix = "-p2-uds.dat";

histogram hMCDF = histogram(40.0,600.0,17,"M-CDF",prefix);
histogram hMyCDF = histogram(11.0,150.0,9,"My-CDF",prefix); 
histogram hyCDFb = histogram(0.0,2.4,8,"y-CDFb",prefix);

void save1800(void)
{
 
 hMCDF.saveashistogram();
 hMyCDF.saveashistogram();
 hyCDFb.saveashistogram();
}

void put1800(vector4 k1, vector4 k2, vector4 p1, vector4 p2, double dsigma)
{
   vector4 q, pl, pnl;
choose2(p1,p2,pl,pnl);

  q.set(0,p1.get(0) + p2.get(0));
  q.set(1,p1.get(1) + p2.get(1));
  q.set(2,p1.get(2) + p2.get(2));
  q.set(3,p1.get(3) + p2.get(3));

  double pt = q.transverse_momentum();
  double y = q.rapidity();
  double eta = q.pseudo_rapidity();
  double etal = pl.pseudo_rapidity();
  double etanl = pnl.pseudo_rapidity();
  double Etl = pl.transverse_energy();
  double Etnl = pnl.transverse_energy();
  double ptl = pl.transverse_momentum();
  double ptnl = pnl.transverse_momentum();
  double M = inv_mass(p1,p2);
  double dphi = phi(p1,p2);

  double F = 1.0;// + pow(M,3);
  dsigma = dsigma/F;

 
    
 if ( (M >= 40.0)&&
  (fabs(etal) <= 4.2) &&(fabs(etanl) <= 4.2)&&
  (Etl >= 15.0) &&(Etnl >= 15.0)&&
  (ptl >= 10.0) &&(ptnl >= 10.0))
 {
   hMCDF.put(M,dsigma);
 }


if ( (M >= 11.0) && (M <= 150.0)&& 
  (fabs(etal) <= 1.0)&& (fabs(etanl) <= 1.0) && 
  (ptl >=4.8)&& (ptnl >= 4.8)&& 
  (Etl >= 5.0)&& (Etnl >= 5.0))
  {
  hMyCDF.put(M,dsigma/2.0);
  }

 
if ( (M >= 116.0) &&  (fabs(y) <= 2.4)
      && (Etl >= 15.0)&& (Etnl >= 15.0 ))
 { 
  
   hyCDFb.put(y,dsigma);
 }        
}

void format(void)
{

hMCDF.format(0,40.0,50.0);
hMCDF.format(1,50.0,60.0);
hMCDF.format(2,60.0,70.0);
hMCDF.format(3,70.0,78.0);
hMCDF.format(4,78.0,86.0);
hMCDF.format(5,86.0,88.0);
hMCDF.format(6,88.0,90.0);
hMCDF.format(7,90.0,92.0);
hMCDF.format(8,92.0,94.0);
hMCDF.format(9,94.0,100.0);
hMCDF.format(10,100.0,105.0);
hMCDF.format(11,105.0,125.0);
hMCDF.format(12,125.0,150.0);
hMCDF.format(13,150.0,200.0);
hMCDF.format(14,200.0,300.0);
hMCDF.format(15,300.0,400.0);
hMCDF.format(16,400.0,600.0);

hMyCDF.format(0,11.0,15.0);
hMyCDF.format(1,15.0,20.0);
hMyCDF.format(2,20.0,30.0);
hMyCDF.format(3,30.0,40.0);
hMyCDF.format(4,40.0,50.0);
hMyCDF.format(5,50.0,60.0);
hMyCDF.format(6,60.0,70.0);
hMyCDF.format(7,70.0,110.0);
hMyCDF.format(8,110.0,150.0);
  
}

double K(double pt, double M)
{
  double CF = 4.0/3.0;
  double mu2 = pow(pt,4.0/3.0)*pow(M,2.0/3.0);

  return exp(CF*aQCDSS(mu2)*sqr(PI)/(2.0*PI));
}

