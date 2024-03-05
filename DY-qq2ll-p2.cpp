#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "../lib/pdf.h"
#include "../lib/kmr.h"
#include "../lib/nrupdf.h"
#include "../lib/updf.h"
#include "../lib/jet.h"
#include "../lib/utils.h"
#include "../lib/vegas.h"
#include "../lib/vector4.h"
#include "../lib/interface.h"

using namespace std;

extern model MODEL;
extern unsigned int PDF;
extern double ml, S, scale_factor;
extern double Mmin, Mmax;
extern double qtmin, qtmax;

nrupdf uv(DOUBLE_SCALE);
nrupdf dv(DOUBLE_SCALE);

void title()
{
  display(MODEL,PDF);

  cout << "       Ecm = " << sqrt(S) << " GeV" << endl;
  cout << "       scale factor = " << scale_factor << endl;
  cout << "       number of flavours = " << NF << endl;
  cout << "       QCD parameter = " << LQCD*1000.0 << " MeV" << endl;
}

pdf grv94(GRV94);
pdf mstw(MSTW2008,"mstw2008lo.00.dat");
pdf *KMR_INPUT;
extern double K(double pt2, double M);
extern void save1800(void);
extern void put1800(vector4 k1, vector4 k2, vector4 p1, vector4 p2, double dsigma);
extern void format(void);

extern double qg2llq(vector4 k1, vector4 k2,vector4 l1, vector4 p1, vector4 p2, vector4 p3,  vector4 e2, double ml2, double I3Q, double eQ, double muR2);
vector4 beam1, beam2;
double calcqq2ll(double x[], double wgt)
{
  double pdfs;
  double M, ml2, muR2, muF2, J, qt;
  double pt12, phi1, y1, mt1;
  double pt22, phi2, y2, mt2;
  double pt32, y3, mt3;
  double kt12, theta1;
  double kt22, theta2;
  double x1, x2;
  double u1, u2, d1, d2;
  double _u , _d, _s, _c,_b,_u1,_d1, _c1,_s1,_ub1,_db1,_cb1,_sb1, _sea1, _sea2, _t,_u2,_d2,_c2,_s2,_ub2,_db2,_cb2,_sb2,_gluon2,_gluon1;
  double s,s2,s1, t, u;
  vector4 p1, p2,p3, k1, k2,kt2, l1, l2,e2, q;

  double dsigma = 0.0;
  double F = 1.0;
 
  pt12 = sqr(exp(x[1]));
  pt22 = sqr(exp(x[2]));
  y1 = x[3];
  y2 = x[4];
  y3 = x[4];

  kt12 = 1.0e-10;
  kt22 = 1.0e-10;
  J = 4.0*pt12*pt22;  

  if (  MODEL == kT ) 
  {
    kt12 = sqr(exp(x[6]));
    kt22 = sqr(exp(x[7]));

    J = 4.0*kt12*kt22*J;
  }  

  ml2 = sqr(ml);
  mt1 = sqrt(ml2 + pt12);
  mt2 = sqrt(ml2 + pt22);

  phi1 = random(0.0,2.0*PI);
  phi2 = random(0.0,2.0*PI);
  theta1 = random(0.0,2.0*PI);
  theta2 = random(0.0,2.0*PI);

  p1.set(mt1*cosh(y1),sqrt(pt12)*cos(phi1),sqrt(pt12)*sin(phi1),mt1*sinh(y1));
  p2.set(mt2*cosh(y2),sqrt(pt22)*cos(phi2),sqrt(pt22)*sin(phi2),mt2*sinh(y2));

  k1.set(1,sqrt(kt12)*cos(theta1));
  k1.set(2,sqrt(kt12)*sin(theta1));

  k2.set(1,sqrt(kt22)*cos(theta2));
  k2.set(2,sqrt(kt22)*sin(theta2));

  p3.set(1,k1.get(1) + k2.get(1) - p1.get(1) - p2.get(1));
  p3.set(2,k1.get(2) + k2.get(2) - p1.get(2) - p2.get(2));
 
  pt32 = p3.transverse_momentum2();
  mt3 = sqrt( pt32);

  p3.set(0,mt3*cosh(y3));
  p3.set(3,mt3*sinh(y3));

  q.set(0,p1.get(0) + p2.get(0));
  q.set(1,p1.get(1) + p2.get(1));
  q.set(2,p1.get(2) + p2.get(2));
  q.set(3,p1.get(3) + p2.get(3));

  qt = q.transverse_momentum();

  M = inv_mass(p1,p2);

  if ( M < Mmin ) goto Exit;
  if ( M > Mmax ) goto Exit;

  if (qt < qtmin) goto Exit;
  if (qt > qtmax) goto Exit; 

  x1 = ( mt1*exp( + y1) + mt2*exp( + y2) + mt3*exp( + y3) )/sqrt(S); if (x1 >= 1.0 - 1.0e-10) goto Exit;
  x2 = ( mt1*exp( - y1) + mt2*exp( - y2) + mt3*exp( - y3) )/sqrt(S); if (x2 >= 1.0 - 1.0e-10) goto Exit;

  k1.set(0,x1*sqrt(S)/2.0);
  k1.set(3,x1*sqrt(S)/2.0);

  k2.set(0,x2*sqrt(S)/2.0);
  k2.set(3, - x2*sqrt(S)/2.0);

  e2.set(0.0,cos(theta2),sin(theta2),0.0);

  l1.set(x1*0.5*sqrt(S),0.0,0.0,x1*0.5*sqrt(S));
  l2.set(x2*0.5*sqrt(S),0.0,0.0, - x2*0.5*sqrt(S));

  s =  dot(k1,k1) + dot(k2,k2) + 2.0*dot(k1,k2);
  s2 = dot(p1,p1) + dot(p2,p2) + 2.0*dot(p1,p2);
  s1 = dot(k2,k2) - 2.0*dot(k2,p3);//p3p3=0

  if (sqrt(s2) < 0.0) goto Exit;
  if (sqrt(s2) > sqrt(s)) goto Exit;
  if (G(s,s1,s2,dot(k1,k1),dot(k2,k2),0.0) > 0.0) goto Exit;
  if ( L(s,dot(k1,k1),dot(k2,k2)) < 0.0 ) goto Exit;
  if ( L(s,s2,0.0) < 0.0 ) goto Exit;

  if (isolated_cone(p1,p2) == FALSE) goto Exit;
  if (separated_cone(p1,p2) == FALSE) goto Exit;

  if (isolated_cone(p1,p3) == FALSE) goto Exit;
  if (isolated_cone(p2,p3) == FALSE) goto Exit;

  muR2 = scale_factor*sqr(M); 
  muF2 = scale_factor*sqr(M); 

  if ( (MODEL == kT) AND (PDF == KMR) ) 
  { 
    KMR_INPUT = &mstw; 

_u1 = kmr_u(x1,kt12,muF2);_d1 = kmr_d(x1,kt12,muF2);_c1 = kmr_c(x1,kt12,muF2);_s1 = kmr_s(x1,kt12,muF2); 
	
_ub1 = kmr_ub(x1,kt12,muF2);_db1 =kmr_db(x1,kt12,muF2);_cb1 = kmr_cb(x1,kt12,muF2);_sb1= kmr_sb(x1,kt12,muF2);
_gluon1 = kmr_gluon(x1,kt12,muF2);_gluon2 = kmr_gluon(x2,kt22,muF2);
	
_u2 = kmr_u(x2,kt22,muF2);_d2 = kmr_d(x2,kt22,muF2);_c2 = kmr_c(x2,kt22,muF2);_s2 = kmr_s(x2,kt22,muF2); 
	
_ub2 = kmr_ub(x2,kt22,muF2);_db2 =kmr_db(x2,kt22,muF2);_cb2 = kmr_cb(x2,kt22,muF2);_sb2= kmr_sb(x2,kt22,muF2);
    _u = kmr_ub(x1,kt12,muF2) + kmr_u(x1,kt12,muF2); 
    _d = kmr_db(x1,kt12,muF2) + kmr_d(x1,kt12,muF2);
    _c = kmr_cb(x1,kt12,muF2) + kmr_c(x1,kt12,muF2);
    _s = kmr_sb(x1,kt12,muF2) + kmr_s(x1,kt12,muF2);
    _b = kmr_bb(x1,kt12,muF2) + kmr_b(x1,kt12,muF2);
    _t = kmr_tb(x1,kt12,muF2) + kmr_t(x1,kt12,muF2);
    
  }

  dsigma = J*(_u * _gluon2*qg2llq(k1,k2,l1,p1,p2,p3,e2,ml2,I3uct,eU,muR2)+_d * _gluon2*qg2llq(k1,k2,l1,p1,p2,p3,e2,ml2,I3dsb,eD,muR2)+_s * _gluon2*qg2llq(k1,k2,l1,p1,p2,p3,e2,ml2,I3dsb,eS,muR2))/(256.0*pow(PI,3)*sqr(x1*x2*S)); 


  dsigma = 2.0*dsigma; // qg + gq
  dsigma = K(q.transverse_momentum(),M)*GeV2pb*dsigma;
     
  //F = 1.0 + pow(M,3);
  if (init > 0) 
  {
    if ( S == sqr(1800.0) ) put1800(k1,k2,p1,p2,F*wgt*dsigma/itmx);

  }

Exit:
  return F*dsigma;
}

main(void)
{
  title();
  randomize();
  format();

  if ( MODEL == LO )
  {
    ndim = 5;

    regn[1] = log(1.0e-03); regn[6] = log(1.0e+03);
    regn[2] = log(1.0e-03); regn[7] = log(1.0e+03);
    regn[3] = - 6.5;        regn[8] = 6.5;
    regn[4] = - 6.5;        regn[9] = 6.5;
    regn[5] = - 6.5;        regn[10] = 6.5;
  }
    else
  {
    ndim = 7;

    regn[1] = log(1.0e-03); regn[8] = log(1.0e+03);
    regn[2] = log(1.0e-03); regn[9] = log(1.0e+03);
    regn[3] = - 6.5;        regn[10] = 6.5;
    regn[4] = - 6.5;        regn[11] = 6.5;
    regn[5] = - 6.5;        regn[12] = 6.5;
    regn[6] = log(0.001);   regn[13] = log(1.0e+03);
    regn[7] = log(0.001);   regn[14] = log(1.0e+03);
  }

  init = 0;
  itmx = 30;
  calls = 200000;

  vegas(regn,ndim,calcqq2ll,init,calls,itmx,0,&tgral,&sd,&chi2a);
  
  init = 1;
  itmx = 60;
  calls = 1000000;

  vegas(regn,ndim,calcqq2ll,init,calls,itmx,0,&tgral,&sd,&chi2a);
  cout << "\nqq2ll contribution: " << tgral << " (sd = " << sd << ", chi2a = " << chi2a << ")\n";


  if ( S == sqr(1800.0) ) save1800();
}


