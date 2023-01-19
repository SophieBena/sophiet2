
#include "duanzhang.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define R 8.3143

/* pure EOSterms for 0 to 0.2 GPa */
static const double H2OLa1  =  4.38269941E-02;
static const double H2OLa2  = -1.68244362E-01;
static const double H2OLa3  = -2.36923373E-01;
static const double H2OLa4  =  1.13027462E-02;
static const double H2OLa5  = -7.67764181E-02;
static const double H2OLa6  =  9.71820593E-02;
static const double H2OLa7  =  6.62674916E-05;
static const double H2OLa8  =  1.06637349E-03;
static const double H2OLa9  = -1.23265258E-03;
static const double H2OLa10 = -8.93953948E-06;
static const double H2OLa11 = -3.88124606E-05;
static const double H2OLa12 =  5.61510206E-05;
static const double H2OLa   =  7.51274488E-03;
static const double H2OLb   =  2.51598931E+00;
static const double H2OLc   =  3.94000000E-02;

static const double CO2La1  =  1.14400435E-01;
static const double CO2La2  = -9.38526684E-01;
static const double CO2La3  =  7.21857006E-01;
static const double CO2La4  =  8.81072902E-03;
static const double CO2La5  =  6.36473911E-02;
static const double CO2La6  = -7.70822213E-02;
static const double CO2La7  =  9.01506064E-04;
static const double CO2La8  = -6.81834166E-03;
static const double CO2La9  =  7.32364258E-03;
static const double CO2La10 = -1.10288237E-04;
static const double CO2La11 =  1.26524193E-03;
static const double CO2La12 = -1.49730823E-03;
static const double CO2La   =  7.81940730E-03;
static const double CO2Lb   = -4.22918013E+00;
static const double CO2Lc   =  1.58500000E-01;

/* pure EOS terms for 0.2 to 10 GPa */
static const double H2OHa1  =  4.68071541E-02;
static const double H2OHa2  = -2.81275941E-01;
static const double H2OHa3  = -2.43926365E-01;
static const double H2OHa4  =  1.10016958E-02;
static const double H2OHa5  = -3.86603525E-02;
static const double H2OHa6  =  9.30095461E-02;
static const double H2OHa7  = -1.15747171E-05;
static const double H2OHa8  =  4.19873848E-04;
static const double H2OHa9  = -5.82739501E-04;
static const double H2OHa10 =  1.00936000E-06;
static const double H2OHa11 = -1.01713593E-05;
static const double H2OHa12 =  1.63934213E-05;
static const double H2OHa   = -4.49505919E-02;
static const double H2OHb   = -3.15028174E-01;
static const double H2OHc   =  1.25000000E-02;

static const double CO2Ha1  =  5.72573440E-03;
static const double CO2Ha2  =  7.94836769E+00;
static const double CO2Ha3  = -3.84236281E+01;
static const double CO2Ha4  =  3.71600369E-02;
static const double CO2Ha5  = -1.92888994E+00;
static const double CO2Ha6  =  6.64254770E+00;
static const double CO2Ha7  = -7.02203950E-06;
static const double CO2Ha8  =  1.77093234E-02;
static const double CO2Ha9  = -4.81892026E-02;
static const double CO2Ha10 =  3.88344869E-06;
static const double CO2Ha11 = -5.54833167E-04;
static const double CO2Ha12 =  1.70489748E-03;
static const double CO2Ha   = -4.13039220E-01;
static const double CO2Hb   = -8.47988634E+00;
static const double CO2Hc   =  2.80000000E-02;

/* H2O, critical constants, K, bars, J/bar */

static const double H2OTc = 647.25;
static const double H2OPc = 221.19;
#define H2OVc (8.314467*H2OTc/H2OPc)

/* CO2 */

static const double CO2Tc = 304.1282;
static const double CO2Pc = 73.773;
#define CO2Vc (8.314467*CO2Tc/CO2Pc)

#define H2O 0
#define CO2 1

static const double idealCoeff[13][2] = {
    {  3.10409601236035e+01, -0.18188731e+01 },
    { -3.91422080460869e+01,  0.12903022e+02 },
    {  3.79695277233575e+01, -0.96634864e+01 },
    { -2.18374910952284e+01,  0.42251879e+01 },
    {  7.42251494566339e+00, -0.10421640e+01 },
    { -1.38178929609470e+00,  0.12683515e+00 },
    {  1.08807067571454e-01, -0.49939675e-02 },
    { -1.20771176848589e+01,  0.24950242e+01 },
    {  3.39105078851732e+00, -0.82723750e+00 },
    { -5.84520979955060e-01,  0.15372481e+00 },
    {  5.89930846488082e-02, -0.15861243e-01 },
    { -3.12970001415882e-03,  0.86017150e-03 },
    {  6.57460740981757e-05, -0.19222165e-04 }
};

static double powSum(double a, double fa, double b, double fb) {
    double sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    return (sum >= 0.0) ? pow(sum/(fa+fb), 3.0) : -pow(-sum/(fa+fb), 3.0);
}

static double DpowSum(double a, double fa, double b, double fb, double da, double db) {
    double sum = 0.0, dsum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum = (sum >= 0.0) ? pow(sum/(fa+fb), 2.0) : pow(-sum/(fa+fb), 2.0);
    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    return 3.0*sum*dsum/(fa+fb);
}

static double D2powSum(double a, double fa, double b, double fb, double da, double db, double d2a, double d2b) {
    double sum = 0.0, dsum = 0.0, d2sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum /= fa + fb;

    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    dsum /= fa + fb;

    d2sum += (a >= 0.0) ? d2a*fa/pow(a, 2.0/3.0)/3.0 : d2a*fa/pow(-a, 2.0/3.0)/3.0;
    d2sum += (b >= 0.0) ? d2b*fb/pow(b, 2.0/3.0)/3.0 : d2b*fb/pow(-b, 2.0/3.0)/3.0;
    d2sum += (a >= 0.0) ? -2.0*da*da*fa/pow(a, 5.0/3.0)/9.0 : 2.0*da*da*fa/pow(-a, 5.0/3.0)/9.0;
    d2sum += (b >= 0.0) ? -2.0*db*db*fb/pow(b, 5.0/3.0)/9.0 : 2.0*db*db*fb/pow(-b, 5.0/3.0)/9.0;
    d2sum /= fa + fb;

    return 6.0*sum*dsum*dsum + 3.0*sum*sum*d2sum;
}

static double D3powSum(double a, double fa, double b, double fb, double da, double db, double d2a, double d2b, double d3a, double d3b) {
    double sum = 0.0, dsum = 0.0, d2sum = 0.0, d3sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    sum /= fa + fb;

    dsum += (a >= 0.0) ? da*fa/pow(a, 2.0/3.0)/3.0 : da*fa/pow(-a, 2.0/3.0)/3.0;
    dsum += (b >= 0.0) ? db*fb/pow(b, 2.0/3.0)/3.0 : db*fb/pow(-b, 2.0/3.0)/3.0;
    dsum /= fa + fb;

    d2sum += (a >= 0.0) ? d2a*fa/pow(a, 2.0/3.0)/3.0 : d2a*fa/pow(-a, 2.0/3.0)/3.0;
    d2sum += (b >= 0.0) ? d2b*fb/pow(b, 2.0/3.0)/3.0 : d2b*fb/pow(-b, 2.0/3.0)/3.0;
    d2sum += (a >= 0.0) ? -2.0*da*da*fa/pow(a, 5.0/3.0)/9.0 : 2.0*da*da*fa/pow(-a, 5.0/3.0)/9.0;
    d2sum += (b >= 0.0) ? -2.0*db*db*fb/pow(b, 5.0/3.0)/9.0 : 2.0*db*db*fb/pow(-b, 5.0/3.0)/9.0;
    d2sum /= fa + fb;

    d3sum += (a >= 0.0) ? d3a*fa/pow(a, 2.0/3.0)/3.0 : d3a*fa/pow(-a, 2.0/3.0)/3.0;
    d3sum += (b >= 0.0) ? d3b*fb/pow(b, 2.0/3.0)/3.0 : d3b*fb/pow(-b, 2.0/3.0)/3.0;
    d3sum += (a >= 0.0) ? -6.0*da*d2a*fa/pow(a, 5.0/3.0)/9.0 : 6.0*da*d2a*fa/pow(-a, 5.0/3.0)/9.0;
    d3sum += (b >= 0.0) ? -6.0*db*d2b*fb/pow(b, 5.0/3.0)/9.0 : 6.0*db*d2b*fb/pow(-b, 5.0/3.0)/9.0;
    d3sum += (a >= 0.0) ? 10.0*da*da*da*fa/pow(a, 8.0/3.0)/27.0 : 10.0*da*da*da*fa/pow(-a, 8.0/3.0)/27.0;
    d3sum += (b >= 0.0) ? 10.0*db*db*db*fb/pow(b, 8.0/3.0)/27.0 : 10.0*db*db*db*fb/pow(-b, 8.0/3.0)/27.0;
    d3sum /= fa + fb;

    return 6.0*dsum*dsum*dsum + 18.0*sum*dsum*d2sum + 3.0*sum*sum*d3sum;
}

static void BVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *bv, double bvPrime[2], double *dbvdt, double *d2bvdt2, double *d3bvdt3, double dbvPrimedt[2], double d2bvPrimedt2[2], double d3bvPrimedt3[2],
                             double b2vPrime2[2][2]) {
    double bEnd[2], dbEnddt[2], d2bEnddt2[2], d3bEnddt3[2];
    double H2OTr    = t/H2OTc;
    double CO2Tr    = t/CO2Tc;
    double k1       = 0.0;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;
    double dk1dt    = 0.0;
    double d2k1dt2  = 0.0;
    double d3k1dt3  = 0.0;

    if (useLowPcoeff) {
        bEnd[H2O] = H2OLa1 + H2OLa2/H2OTr/H2OTr + H2OLa3/H2OTr/H2OTr/H2OTr;
        bEnd[CO2] = CO2La1 + CO2La2/CO2Tr/CO2Tr + CO2La3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 3.131 - 5.0624e-3*t + 1.8641e-6*t*t - 31.409/t;
        dbEnddt[H2O] = - 2.0*H2OLa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dbEnddt[CO2] = - 2.0*CO2La2*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La3*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk1dt = - 5.0624e-3 + 2.0*1.8641e-6*t + 31.409/t/t;
        d2bEnddt2[H2O] = 6.0*H2OLa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2bEnddt2[CO2] = 6.0*CO2La2*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La3*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k1dt2 = 2.0*1.8641e-6 - 2.0*31.409/t/t/t;
        d3bEnddt3[H2O] = - 24.0*H2OLa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3bEnddt3[CO2] = - 24.0*CO2La2*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La3*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k1dt3 =  6.0*31.409/t/t/t/t;
    } else {
        bEnd[H2O] = H2OHa1 + H2OHa2/H2OTr/H2OTr + H2OHa3/H2OTr/H2OTr/H2OTr;
        bEnd[CO2] = CO2Ha1 + CO2Ha2/CO2Tr/CO2Tr + CO2Ha3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 9.034 - 7.9212e-3*t + 2.3285e-6*t*t - 2.4221e3/t;
        dbEnddt[H2O] = - 2.0*H2OHa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dbEnddt[CO2] = - 2.0*CO2Ha2*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha3*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk1dt = - 7.9212e-3 + 2.0*2.3285e-6*t + 2.4221e3/t/t;
        d2bEnddt2[H2O] = 6.0*H2OHa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2bEnddt2[CO2] = 6.0*CO2Ha2*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha3*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k1dt2 = 2.0*2.3285e-6 - 2.0*2.4221e3/t/t/t;
        d3bEnddt3[H2O] = - 24.0*H2OHa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3bEnddt3[CO2] = - 24.0*CO2Ha2*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha3*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k1dt3 = 6.0*2.4221e3/t/t/t/t;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *bv = 0.0;
        bvPrime[H2O] = 0.0;
        bvPrime[CO2] = 0.0;

        b2vPrime2[H2O][H2O] = 0.0;
        b2vPrime2[H2O][CO2] = 0.0;
        b2vPrime2[CO2][CO2] = 0.0;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = 0.0;
        *d2bvdt2 = 0.0;
        *d3bvdt3 = 0.0;

        dbvPrimedt[H2O] = 0.0;
        dbvPrimedt[CO2] = 0.0;

        d2bvPrimedt2[H2O] = 0.0;
        d2bvPrimedt2[CO2] = 0.0;

        d3bvPrimedt3[H2O] = 0.0;
        d3bvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *bv = bEnd[H2O]*H2OVc;
        bvPrime[H2O] = 2.0*bEnd[H2O]*H2OVc;
        bvPrime[CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = dbEnddt[H2O]*H2OVc;
        *d2bvdt2 = d2bEnddt2[H2O]*H2OVc;
        *d3bvdt3 = d3bEnddt3[H2O]*H2OVc;

        dbvPrimedt[H2O] = 2.0*dbEnddt[H2O]*H2OVc;
        dbvPrimedt[CO2] = 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                               + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                               )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        d2bvPrimedt2[H2O] = 2.0*d2bEnddt2[H2O]*H2OVc;
        d2bvPrimedt2[CO2] = 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                 + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

        d3bvPrimedt3[H2O] = 2.0*d3bEnddt3[H2O]*H2OVc;
        d3bvPrimedt3[CO2] = 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                 + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                 + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *bv = bEnd[CO2]*CO2Vc;
        bvPrime[H2O] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        bvPrime[CO2] = 2.0*bEnd[CO2]*CO2Vc;

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt   = dbEnddt[CO2]*CO2Vc;
        *d2bvdt2 = d2bEnddt2[CO2]*CO2Vc;
        *d3bvdt3 = d3bEnddt3[CO2]*CO2Vc;

        dbvPrimedt[H2O] = 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                               + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                               )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        dbvPrimedt[CO2] = 2.0*dbEnddt[CO2]*CO2Vc;

        d2bvPrimedt2[H2O] = 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                 + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        d2bvPrimedt2[CO2] = 2.0*d2bEnddt2[CO2]*CO2Vc;

        d3bvPrimedt3[H2O] = 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                 + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                 + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                 +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                 )*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        d3bvPrimedt3[CO2] = 2.0*d3bEnddt3[CO2]*CO2Vc;

    } else {
        *bv  = 0.0;
        *bv += bEnd[H2O]*H2OVc*x[H2O]*x[H2O];
        *bv += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *bv += bEnd[CO2]*CO2Vc*x[CO2]*x[CO2];

        bvPrime[H2O]  = 0.0;
        bvPrime[H2O] += 2.0*bEnd[H2O]*H2OVc*x[H2O];
        bvPrime[H2O] += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        bvPrime[CO2]  = 0.0;
        bvPrime[CO2] += 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        bvPrime[CO2] += 2.0*bEnd[CO2]*CO2Vc*x[CO2];

        b2vPrime2[H2O][H2O] = 2.0*bEnd[H2O]*H2OVc;
        b2vPrime2[H2O][CO2] = 2.0*powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*k1*powSum(H2OVc, 1.0, CO2Vc, 1.0);
        b2vPrime2[CO2][CO2] = 2.0*bEnd[CO2]*CO2Vc;
        b2vPrime2[CO2][H2O] = b2vPrime2[H2O][CO2];

        *dbvdt  = 0.0;
        *dbvdt += dbEnddt[H2O]*H2OVc*x[H2O]*x[H2O];
        *dbvdt += 2.0*(DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1 + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                       )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *dbvdt += dbEnddt[CO2]*CO2Vc*x[CO2]*x[CO2];

        *d2bvdt2  = 0.0;
        *d2bvdt2 += d2bEnddt2[H2O]*H2OVc*x[H2O]*x[H2O];
        *d2bvdt2 += 2.0*(  D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                         + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                         + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                         )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *d2bvdt2 += d2bEnddt2[CO2]*CO2Vc*x[CO2]*x[CO2];

        *d3bvdt3  = 0.0;
        *d3bvdt3 += d3bEnddt3[H2O]*H2OVc*x[H2O]*x[H2O];
        *d3bvdt3 += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                         + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                         + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                         +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                         )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O]*x[CO2];
        *d3bvdt3 += d3bEnddt3[CO2]*CO2Vc*x[CO2]*x[CO2];

        dbvPrimedt[H2O]  = 0.0;
        dbvPrimedt[H2O] += 2.0*dbEnddt[H2O]*H2OVc*x[H2O];
        dbvPrimedt[H2O] += 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                                + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                                )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        dbvPrimedt[CO2]  = 0.0;
        dbvPrimedt[CO2] += 2.0*(  DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*k1
                                + powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*dk1dt
                                )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        dbvPrimedt[CO2] += 2.0*dbEnddt[CO2]*CO2Vc*x[CO2];

        d2bvPrimedt2[H2O]  = 0.0;
        d2bvPrimedt2[H2O] += 2.0*d2bEnddt2[H2O]*H2OVc*x[H2O];
        d2bvPrimedt2[H2O] += 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                  + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        d2bvPrimedt2[CO2]  = 0.0;
        d2bvPrimedt2[CO2] += 2.0*(      D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*k1
                                  + 2.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*dk1dt
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d2k1dt2
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        d2bvPrimedt2[CO2] += 2.0*d2bEnddt2[CO2]*CO2Vc*x[CO2];

        d3bvPrimedt3[H2O]  = 0.0;
        d3bvPrimedt3[H2O] += 2.0*d3bEnddt3[H2O]*H2OVc*x[H2O];
        d3bvPrimedt3[H2O] += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                  + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                  + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[CO2];

        d3bvPrimedt3[CO2]  = 0.0;
        d3bvPrimedt3[CO2] += 2.0*(      D3powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2], d3bEnddt3[H2O], d3bEnddt3[CO2])*k1
                                  + 3.0*D2powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2], d2bEnddt2[H2O], d2bEnddt2[CO2])*dk1dt
                                  + 3.0*DpowSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0, dbEnddt[H2O], dbEnddt[CO2])*d2k1dt2
                                  +     powSum(bEnd[H2O], 1.0, bEnd[CO2], 1.0)*d3k1dt3
                                  )*powSum(H2OVc, 1.0, CO2Vc, 1.0)*x[H2O];
        d3bvPrimedt3[CO2] += 2.0*d3bEnddt3[CO2]*CO2Vc*x[CO2];

    }

    return;
}

static void CVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *cv, double cvPrime[2], double *dcvdt, double *d2cvdt2, double *d3cvdt3, double dcvPrimedt[2], double d2cvPrimedt2[2], double d3cvPrimedt3[2],
                             double c2vPrime2[2][2]) {
    double cEnd[2], dcEnddt[2], d2cEnddt2[2], d3cEnddt3[2];
    double H2OTr    = t/H2OTc;
    double CO2Tr    = t/CO2Tc;
    double k2       = 0.0;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;
    double dk2dt    = 0.0;
    double d2k2dt2  = 0.0;
    double d3k2dt3  = 0.0;

    if (useLowPcoeff) {
        cEnd[H2O] = H2OLa4 + H2OLa5/H2OTr/H2OTr + H2OLa6/H2OTr/H2OTr/H2OTr;
        cEnd[CO2] = CO2La4 + CO2La5/CO2Tr/CO2Tr + CO2La6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -46.646 + 4.2877e-2*t - 1.0892e-5*t*t + 1.5782e4/t;
        dcEnddt[H2O] = - 2.0*H2OLa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dcEnddt[CO2] = - 2.0*CO2La5*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La6*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk2dt = 4.2877e-2 - 2.0*1.0892e-5*t - 1.5782e4/t/t;
        d2cEnddt2[H2O] = 6.0*H2OLa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2cEnddt2[CO2] = 6.0*CO2La5*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La6*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k2dt2 = - 2.0*1.0892e-5 + 2.0*1.5782e4/t/t/t;
        d3cEnddt3[H2O] = -24.0*H2OLa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3cEnddt3[CO2] = -24.0*CO2La5*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La6*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k2dt3 = - 3.0*2.0*1.5782e4/t/t/t/t;
    } else {
        cEnd[H2O] = H2OHa4 + H2OHa5/H2OTr/H2OTr + H2OHa6/H2OTr/H2OTr/H2OTr;
        cEnd[CO2] = CO2Ha4 + CO2Ha5/CO2Tr/CO2Tr + CO2Ha6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -1.068 + 1.8756e-3*t - 4.9371e-7*t*t + 6.6180e2/t;
        dcEnddt[H2O] = - 2.0*H2OHa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dcEnddt[CO2] = - 2.0*CO2Ha5*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha6*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        dk2dt = 1.8756e-3 - 2.0*4.9371e-7*t - 6.6180e2/t/t;
        d2cEnddt2[H2O] = 6.0*H2OHa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2cEnddt2[CO2] = 6.0*CO2Ha5*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha6*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2k2dt2 = - 2.0*4.9371e-7 + 2.0*6.6180e2/t/t/t;
        d3cEnddt3[H2O] = -24.0*H2OHa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3cEnddt3[CO2] = -24.0*CO2Ha5*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha6*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3k2dt3 = - 3.0*2.0*6.6180e2/t/t/t/t;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *cv = 0.0;
        cvPrime[H2O] = 0.0;
        cvPrime[CO2] = 0.0;

        c2vPrime2[H2O][H2O] = 0.0;
        c2vPrime2[H2O][CO2] = 0.0;
        c2vPrime2[CO2][CO2] = 0.0;
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = 0.0;
        *d2cvdt2 = 0.0;
        *d3cvdt3 = 0.0;

        dcvPrimedt[H2O] = 0.0;
        dcvPrimedt[CO2] = 0.0;

        d2cvPrimedt2[H2O] = 0.0;
        d2cvPrimedt2[CO2] = 0.0;

        d3cvPrimedt3[H2O] = 0.0;
        d3cvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *cv = cEnd[H2O]*H2OVc*H2OVc;
        cvPrime[H2O] = 3.0*cEnd[H2O]*H2OVc*H2OVc;
        cvPrime[CO2] = 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        c2vPrime2[H2O][H2O] = 6.0*cEnd[H2O]*H2OVc*H2OVc;
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        c2vPrime2[CO2][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = dcEnddt[H2O]*H2OVc*H2OVc;
        *d2cvdt2 = d2cEnddt2[H2O]*H2OVc*H2OVc;
        *d3cvdt3 = d3cEnddt3[H2O]*H2OVc*H2OVc;

        dcvPrimedt[H2O] = 3.0*dcEnddt[H2O]*H2OVc*H2OVc;
        dcvPrimedt[CO2] = 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        d2cvPrimedt2[H2O] = 3.0*d2cEnddt2[H2O]*H2OVc*H2OVc;
        d2cvPrimedt2[CO2] = 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        d3cvPrimedt3[H2O] = 3.0*d3cEnddt3[H2O]*H2OVc*H2OVc;
        d3cvPrimedt3[CO2] = 3.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 9.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 9.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *cv = cEnd[CO2]*CO2Vc*CO2Vc;
        cvPrime[H2O] = 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        cvPrime[CO2] = 3.0*cEnd[CO2]*CO2Vc*CO2Vc;

        c2vPrime2[H2O][H2O] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        c2vPrime2[CO2][CO2] = 6.0*cEnd[CO2]*CO2Vc*CO2Vc;
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt   = dcEnddt[CO2]*CO2Vc*CO2Vc;
        *d2cvdt2 = d2cEnddt2[CO2]*CO2Vc*CO2Vc;
        *d3cvdt3 = d3cEnddt3[CO2]*CO2Vc*CO2Vc;

        dcvPrimedt[H2O] = 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        dcvPrimedt[CO2] = 3.0*dcEnddt[CO2]*CO2Vc*CO2Vc;

        d2cvPrimedt2[H2O] = 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        d2cvPrimedt2[CO2] = 3.0*d2cEnddt2[CO2]*CO2Vc*CO2Vc;

        d3cvPrimedt3[H2O] = 3.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 9.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 9.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        d3cvPrimedt3[CO2] = 3.0*d3cEnddt3[CO2]*CO2Vc*CO2Vc;

    } else {
        *cv  = 0.0;
        *cv += cEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *cv += 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *cv += 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *cv += cEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        cvPrime[H2O]  = 0.0;
        cvPrime[H2O] += 3.0*cEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        cvPrime[H2O] += 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        cvPrime[H2O] += 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        cvPrime[CO2]  = 0.0;
        cvPrime[CO2] += 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        cvPrime[CO2] += 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        cvPrime[CO2] += 3.0*cEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        c2vPrime2[H2O][H2O] = 6.0*cEnd[H2O]*H2OVc*H2OVc*x[H2O]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[CO2];
        c2vPrime2[H2O][CO2] = 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2];
        c2vPrime2[CO2][CO2] = 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]
        + 6.0*cEnd[CO2]*CO2Vc*CO2Vc*x[CO2];
        c2vPrime2[CO2][H2O] = c2vPrime2[H2O][CO2];

        *dcvdt  = 0.0;
        *dcvdt += dcEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *dcvdt += 3.0*(DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2 + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt
                       )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *dcvdt += 3.0*(DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2 + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt
                       )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *dcvdt += dcEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        *d2cvdt2  = 0.0;
        *d2cvdt2 += d2cEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *d2cvdt2 += 3.0*(  D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2
                         + 2.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt
                         + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2
                         )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *d2cvdt2 += 3.0*(  D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2
                         + 2.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt
                         + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2
                         )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *d2cvdt2 += d2cEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        *d3cvdt3  = 0.0;
        *d3cvdt3 += d3cEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *d3cvdt3 += 3.0*(  D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2
                         + 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt
                         + 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2
                         + powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3
                         )*pow( powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *d3cvdt3 += 3.0*(  D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2
                         + 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt
                         + 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2
                         + powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3
                         )*pow( powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *d3cvdt3 += d3cEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        dcvPrimedt[H2O]  = 0.0;
        dcvPrimedt[H2O] += 3.0*dcEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        dcvPrimedt[H2O] += 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        dcvPrimedt[H2O] += 3.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        dcvPrimedt[CO2]  = 0.0;
        dcvPrimedt[CO2] += 3.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        dcvPrimedt[CO2] += 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        dcvPrimedt[CO2] += 3.0*dcEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        d2cvPrimedt2[H2O]  = 0.0;
        d2cvPrimedt2[H2O] += 3.0*d2cEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        d2cvPrimedt2[H2O] += 6.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 12.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        d2cvPrimedt2[H2O] += 3.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 6.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        d2cvPrimedt2[CO2]  = 0.0;
        d2cvPrimedt2[CO2] += 3.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 6.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        d2cvPrimedt2[CO2] += 6.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 12.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        d2cvPrimedt2[CO2] += 3.0*d2cEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        d3cvPrimedt3[H2O]  = 0.0;
        d3cvPrimedt3[H2O] += 3.0*d3cEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        d3cvPrimedt3[H2O] += 6.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        d3cvPrimedt3[H2O] += 3.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 9.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 9.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2]
        + 3.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        d3cvPrimedt3[CO2]  = 0.0;
        d3cvPrimedt3[CO2] += 3.0*D3powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 9.0*D2powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 9.0*DpowSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]
        + 3.0*powSum(cEnd[H2O], 2.0, cEnd[CO2], 1.0)*d3k2dt3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        d3cvPrimedt3[CO2] += 6.0*D3powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2], d3cEnddt3[H2O], d3cEnddt3[CO2])*k2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*D2powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2], d2cEnddt2[H2O], d2cEnddt2[CO2])*dk2dt*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 18.0*DpowSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0, dcEnddt[H2O], dcEnddt[CO2])*d2k2dt2*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]
        + 6.0*powSum(cEnd[H2O], 1.0, cEnd[CO2], 2.0)*d3k2dt3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        d3cvPrimedt3[CO2] += 3.0*d3cEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

    }

    return;
}

static void DVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *dv, double dvPrime[2], double *ddvdt, double *d2dvdt2, double *d3dvdt3, double ddvPrimedt[2], double d2dvPrimedt2[2], double d3dvPrimedt3[2],
                             double d2vPrime2[2][2]) {
    double dEnd[2], ddEnddt[2], d2dEnddt2[2], d3dEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        dEnd[H2O] = H2OLa7 + H2OLa8/H2OTr/H2OTr + H2OLa9/H2OTr/H2OTr/H2OTr;
        dEnd[CO2] = CO2La7 + CO2La8/CO2Tr/CO2Tr + CO2La9/CO2Tr/CO2Tr/CO2Tr;
        ddEnddt[H2O] = - 2.0*H2OLa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        ddEnddt[CO2] = - 2.0*CO2La8*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La9*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2dEnddt2[H2O] = 6.0*H2OLa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2dEnddt2[CO2] = 6.0*CO2La8*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La9*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3dEnddt3[H2O] = - 24.0*H2OLa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3dEnddt3[CO2] = - 24.0*CO2La8*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La9*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        dEnd[H2O] = H2OHa7 + H2OHa8/H2OTr/H2OTr + H2OHa9/H2OTr/H2OTr/H2OTr;
        dEnd[CO2] = CO2Ha7 + CO2Ha8/CO2Tr/CO2Tr + CO2Ha9/CO2Tr/CO2Tr/CO2Tr;
        ddEnddt[H2O] = - 2.0*H2OHa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        ddEnddt[CO2] = - 2.0*CO2Ha8*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha9*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2dEnddt2[H2O] = 6.0*H2OHa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2dEnddt2[CO2] = 6.0*CO2Ha8*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha9*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3dEnddt3[H2O] = - 24.0*H2OHa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3dEnddt3[CO2] = - 24.0*CO2Ha8*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha9*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *dv = 0.0;
        dvPrime[H2O] = 0.0;
        dvPrime[CO2] = 0.0;

        d2vPrime2[H2O][H2O] = 0.0;
        d2vPrime2[H2O][CO2] = 0.0;
        d2vPrime2[CO2][CO2] = 0.0;
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = 0.0;
        *d2dvdt2 = 0.0;
        *d3dvdt3 = 0.0;

        ddvPrimedt[H2O] = 0.0;
        ddvPrimedt[CO2] = 0.0;

        d2dvPrimedt2[H2O] = 0.0;
        d2dvPrimedt2[CO2] = 0.0;

        d3dvPrimedt3[H2O] = 0.0;
        d3dvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *dv = dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        dvPrime[H2O] = 5.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        dvPrime[CO2] = 5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d2vPrime2[H2O][H2O] = 20.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d2vPrime2[CO2][CO2] = 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0);
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        *d2dvdt2 = d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        *d3dvdt3 = d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;

        ddvPrimedt[H2O] = 5.0*ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        ddvPrimedt[CO2] = 5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d2dvPrimedt2[H2O] = 5.0*d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d2dvPrimedt2[CO2] = 5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

        d3dvPrimedt3[H2O] = 5.0*d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc;
        d3dvPrimedt3[CO2] = 5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *dv = dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        dvPrime[H2O] = 5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        dvPrime[CO2] = 5.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2vPrime2[H2O][H2O] = 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0);
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0);
        d2vPrime2[CO2][CO2] = 20.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt   = ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d2dvdt2 = d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d3dvdt3 = d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        ddvPrimedt[H2O] = 5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        ddvPrimedt[CO2] = 5.0*ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2dvPrimedt2[H2O] = 5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d2dvPrimedt2[CO2] = 5.0*d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d3dvPrimedt3[H2O] = 5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0);
        d3dvPrimedt3[CO2] = 5.0*d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

    } else {
        *dv  = 0.0;
        *dv += dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *dv +=  5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *dv += 10.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *dv += 10.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *dv +=  5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *dv += dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        dvPrime[H2O]  =  0.0;
        dvPrime[H2O] +=  5.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        dvPrime[H2O] += 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        dvPrime[H2O] += 30.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        dvPrime[H2O] += 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        dvPrime[H2O] +=  5.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        dvPrime[CO2]  =  0.0;
        dvPrime[CO2] +=  5.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        dvPrime[CO2] += 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        dvPrime[CO2] += 30.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        dvPrime[CO2] += 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        dvPrime[CO2] += 5.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2vPrime2[H2O][H2O] = 20.0*dEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[H2O][CO2] = 20.0*powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[CO2][CO2] = 20.0*powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]
        + 60.0*powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]
        + 60.0*powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]
        + 20.0*dEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];
        d2vPrime2[CO2][H2O] = d2vPrime2[H2O][CO2];

        *ddvdt  = 0.0;
        *ddvdt += ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *ddvdt +=  5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *ddvdt += 10.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *ddvdt += 10.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *ddvdt +=  5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ddvdt += ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d2dvdt2  = 0.0;
        *d2dvdt2 += d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d2dvdt2 +=  5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d2dvdt2 += 10.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d2dvdt2 += 10.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d2dvdt2 +=  5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2dvdt2 += d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d3dvdt3  = 0.0;
        *d3dvdt3 += d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d3dvdt3 +=  5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d3dvdt3 += 10.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d3dvdt3 += 10.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d3dvdt3 +=  5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3dvdt3 += d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        ddvPrimedt[H2O]  =  0.0;
        ddvPrimedt[H2O] +=  5.0*ddEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        ddvPrimedt[H2O] += 20.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        ddvPrimedt[H2O] += 30.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        ddvPrimedt[H2O] += 20.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        ddvPrimedt[H2O] +=  5.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        ddvPrimedt[CO2]  =  0.0;
        ddvPrimedt[CO2] +=  5.0*DpowSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        ddvPrimedt[CO2] += 20.0*DpowSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        ddvPrimedt[CO2] += 30.0*DpowSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        ddvPrimedt[CO2] += 20.0*DpowSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        ddvPrimedt[CO2] += 5.0*ddEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2dvPrimedt2[H2O]  =  0.0;
        d2dvPrimedt2[H2O] +=  5.0*d2dEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2dvPrimedt2[H2O] += 20.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2dvPrimedt2[H2O] += 30.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2dvPrimedt2[H2O] += 20.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2dvPrimedt2[H2O] +=  5.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2dvPrimedt2[CO2]  =  0.0;
        d2dvPrimedt2[CO2] +=  5.0*D2powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2dvPrimedt2[CO2] += 20.0*D2powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2dvPrimedt2[CO2] += 30.0*D2powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2dvPrimedt2[CO2] += 20.0*D2powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2dvPrimedt2[CO2] += 5.0*d2dEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3dvPrimedt3[H2O]  =  0.0;
        d3dvPrimedt3[H2O] +=  5.0*d3dEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3dvPrimedt3[H2O] += 20.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3dvPrimedt3[H2O] += 30.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3dvPrimedt3[H2O] += 20.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3dvPrimedt3[H2O] +=  5.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3dvPrimedt3[CO2]  =  0.0;
        d3dvPrimedt3[CO2] +=  5.0*D3powSum(dEnd[H2O], 4.0, dEnd[CO2], 1.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 1.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3dvPrimedt3[CO2] += 20.0*D3powSum(dEnd[H2O], 3.0, dEnd[CO2], 2.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 2.0), 4.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3dvPrimedt3[CO2] += 30.0*D3powSum(dEnd[H2O], 2.0, dEnd[CO2], 3.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 3.0), 4.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3dvPrimedt3[CO2] += 20.0*D3powSum(dEnd[H2O], 1.0, dEnd[CO2], 4.0, ddEnddt[H2O], ddEnddt[CO2], d2dEnddt2[H2O], d2dEnddt2[CO2], d3dEnddt3[H2O], d3dEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 4.0), 4.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3dvPrimedt3[CO2] += 5.0*d3dEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];

    }

    return;
}

static void EVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *ev, double evPrime[2], double *devdt, double *d2evdt2, double *d3evdt3, double devPrimedt[2], double d2evPrimedt2[2], double d3evPrimedt3[2],
                             double e2vPrime2[2][2]) {
    double eEnd[2], deEnddt[2], d2eEnddt2[2], d3eEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        eEnd[H2O] = H2OLa10 + H2OLa11/H2OTr/H2OTr + H2OLa12/H2OTr/H2OTr/H2OTr;
        eEnd[CO2] = CO2La10 + CO2La11/CO2Tr/CO2Tr + CO2La12/CO2Tr/CO2Tr/CO2Tr;
        deEnddt[H2O] = - 2.0*H2OLa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OLa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        deEnddt[CO2] = - 2.0*CO2La11*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2La12*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2eEnddt2[H2O] = 6.0*H2OLa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OLa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2eEnddt2[CO2] = 6.0*CO2La11*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2La12*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3eEnddt3[H2O] = - 24.0*H2OLa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OLa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3eEnddt3[CO2] = - 24.0*CO2La11*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2La12*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        eEnd[H2O] = H2OHa10 + H2OHa11/H2OTr/H2OTr + H2OHa12/H2OTr/H2OTr/H2OTr;
        eEnd[CO2] = CO2Ha10 + CO2Ha11/CO2Tr/CO2Tr + CO2Ha12/CO2Tr/CO2Tr/CO2Tr;
        deEnddt[H2O] = - 2.0*H2OHa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2OHa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        deEnddt[CO2] = - 2.0*CO2Ha11*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr - 3.0*CO2Ha12*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2eEnddt2[H2O] = 6.0*H2OHa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2OHa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2eEnddt2[CO2] = 6.0*CO2Ha11*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr + 12.0*CO2Ha12*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3eEnddt3[H2O] = - 24.0*H2OHa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2OHa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3eEnddt3[CO2] = - 24.0*CO2Ha11*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr - 60.0*CO2Ha12*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *ev = 0.0;
        evPrime[H2O] = 0.0;
        evPrime[CO2] = 0.0;

        e2vPrime2[H2O][H2O] = 0.0;
        e2vPrime2[H2O][CO2] = 0.0;
        e2vPrime2[CO2][CO2] = 0.0;
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt   = 0.0;
        *d2evdt2 = 0.0;
        *d3evdt3 = 0.0;

        devPrimedt[H2O] = 0.0;
        devPrimedt[CO2] = 0.0;

        d2evPrimedt2[H2O] = 0.0;
        d2evPrimedt2[CO2] = 0.0;

        d3evPrimedt3[H2O] = 0.0;
        d3evPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *ev = eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        evPrime[H2O] = 6.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        evPrime[CO2] = 6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        e2vPrime2[H2O][H2O] =  30.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        e2vPrime2[H2O][CO2] =  30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);
        e2vPrime2[CO2][CO2] =  30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0);
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt = deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        *d2evdt2 = d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        *d3evdt3 = d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;

        devPrimedt[H2O] = 6.0*deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        devPrimedt[CO2] = 6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        d2evPrimedt2[H2O] = 6.0*d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        d2evPrimedt2[CO2] = 6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

        d3evPrimedt3[H2O] = 6.0*d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
        d3evPrimedt3[CO2] = 6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *ev = eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        evPrime[H2O] = 6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        evPrime[CO2] = 6.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        e2vPrime2[H2O][H2O] = 30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0);
        e2vPrime2[H2O][CO2] = 30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        e2vPrime2[CO2][CO2] = 30.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt = deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d2evdt2 = d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;
        *d3evdt3 = d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        devPrimedt[H2O] = 6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        devPrimedt[CO2] = 6.0*deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d2evPrimedt2[H2O] = 6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        d2evPrimedt2[CO2] = 6.0*d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

        d3evPrimedt3[H2O] = 6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0);
        d3evPrimedt3[CO2] = 6.0*d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc;

    } else {
        *ev  = 0.0;
        *ev += eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *ev +=  6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *ev += 15.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *ev += 20.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *ev += 15.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ev +=  6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *ev += eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        evPrime[H2O]  =  0.0;
        evPrime[H2O] +=  6.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        evPrime[H2O] += 30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        evPrime[H2O] += 60.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        evPrime[H2O] += 60.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        evPrime[H2O] += 30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        evPrime[H2O] +=  6.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        evPrime[CO2]  =  0.0;
        evPrime[CO2] +=  6.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        evPrime[CO2] += 30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        evPrime[CO2] += 60.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        evPrime[CO2] += 60.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        evPrime[CO2] += 30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        evPrime[CO2] +=  6.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        e2vPrime2[H2O][H2O] =  30.0*eEnd[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[H2O][CO2] =  30.0*powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0)*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[CO2][CO2] =  30.0*powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0)*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]
        + 120.0*powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0)*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]
        + 180.0*powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0)*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]
        + 120.0*powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]
        +  30.0*eEnd[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        e2vPrime2[CO2][H2O] = e2vPrime2[H2O][CO2];

        *devdt  = 0.0;
        *devdt += deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *devdt +=  6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *devdt += 15.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *devdt += 20.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *devdt += 15.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *devdt +=  6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *devdt += deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d2evdt2  = 0.0;
        *d2evdt2 += d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d2evdt2 +=  6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d2evdt2 += 15.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d2evdt2 += 20.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 += 15.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 +=  6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d2evdt2 += d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        *d3evdt3  = 0.0;
        *d3evdt3 += d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        *d3evdt3 +=  6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        *d3evdt3 += 15.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        *d3evdt3 += 20.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 += 15.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 +=  6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        *d3evdt3 += d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        devPrimedt[H2O]  =  0.0;
        devPrimedt[H2O] +=  6.0*deEnddt[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        devPrimedt[H2O] += 30.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        devPrimedt[H2O] += 60.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        devPrimedt[H2O] += 60.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[H2O] += 30.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[H2O] +=  6.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        devPrimedt[CO2]  =  0.0;
        devPrimedt[CO2] +=  6.0*DpowSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        devPrimedt[CO2] += 30.0*DpowSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        devPrimedt[CO2] += 60.0*DpowSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        devPrimedt[CO2] += 60.0*DpowSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[CO2] += 30.0*DpowSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        devPrimedt[CO2] +=  6.0*deEnddt[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d2evPrimedt2[H2O]  =  0.0;
        d2evPrimedt2[H2O] +=  6.0*d2eEnddt2[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2evPrimedt2[H2O] += 30.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2evPrimedt2[H2O] += 60.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] += 60.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] += 30.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[H2O] +=  6.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        d2evPrimedt2[CO2]  =  0.0;
        d2evPrimedt2[CO2] +=  6.0*D2powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d2evPrimedt2[CO2] += 30.0*D2powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d2evPrimedt2[CO2] += 60.0*D2powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] += 60.0*D2powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] += 30.0*D2powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d2evPrimedt2[CO2] +=  6.0*d2eEnddt2[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

        d3evPrimedt3[H2O]  =  0.0;
        d3evPrimedt3[H2O] +=  6.0*d3eEnddt3[H2O]*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3evPrimedt3[H2O] += 30.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3evPrimedt3[H2O] += 60.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] += 60.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] += 30.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[H2O] +=  6.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];;

        d3evPrimedt3[CO2]  =  0.0;
        d3evPrimedt3[CO2] +=  6.0*D3powSum(eEnd[H2O], 5.0, eEnd[CO2], 1.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 5.0, CO2Vc, 1.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[H2O];
        d3evPrimedt3[CO2] += 30.0*D3powSum(eEnd[H2O], 4.0, eEnd[CO2], 2.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 4.0, CO2Vc, 2.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[H2O]*x[CO2];
        d3evPrimedt3[CO2] += 60.0*D3powSum(eEnd[H2O], 3.0, eEnd[CO2], 3.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 3.0, CO2Vc, 3.0), 5.0)*x[H2O]*x[H2O]*x[H2O]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] += 60.0*D3powSum(eEnd[H2O], 2.0, eEnd[CO2], 4.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 2.0, CO2Vc, 4.0), 5.0)*x[H2O]*x[H2O]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] += 30.0*D3powSum(eEnd[H2O], 1.0, eEnd[CO2], 5.0, deEnddt[H2O], deEnddt[CO2], d2eEnddt2[H2O], d2eEnddt2[CO2], d3eEnddt3[H2O], d3eEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 5.0), 5.0)*x[H2O]*x[CO2]*x[CO2]*x[CO2]*x[CO2];
        d3evPrimedt3[CO2] +=  6.0*d3eEnddt3[CO2]*CO2Vc*CO2Vc*CO2Vc*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2]*x[CO2]*x[CO2];

    }

    return;
}

static void FVcAndDerivative(int useLowPcoeff, double t, double x[2],
                             double *fv, double fvPrime[2], double *dfvdt, double *d2fvdt2, double *d3fvdt3, double dfvPrimedt[2], double d2fvPrimedt2[2], double d3fvPrimedt3[2],
                             double f2vPrime2[2][2]) {
    double fEnd[2], dfEnddt[2], d2fEnddt2[2], d3fEnddt3[2];
    double H2OTr = t/H2OTc;
    double CO2Tr = t/CO2Tc;
    double dH2OTrdt = 1.0/H2OTc;
    double dCO2Trdt = 1.0/CO2Tc;

    if (useLowPcoeff) {
        fEnd[H2O] = H2OLa/H2OTr/H2OTr/H2OTr;
        fEnd[CO2] = CO2La/CO2Tr/CO2Tr/CO2Tr;
        dfEnddt[H2O] = - 3.0*H2OLa*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dfEnddt[CO2] = - 3.0*CO2La*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2fEnddt2[H2O] = 12.0*H2OLa*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2fEnddt2[CO2] = 12.0*CO2La*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3fEnddt3[H2O] = - 60.0*H2OLa*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3fEnddt3[CO2] = - 60.0*CO2La*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    } else {
        fEnd[H2O] = H2OHa/H2OTr/H2OTr/H2OTr;
        fEnd[CO2] = CO2Ha/CO2Tr/CO2Tr/CO2Tr;
        dfEnddt[H2O] = - 3.0*H2OHa*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
        dfEnddt[CO2] = - 3.0*CO2Ha*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d2fEnddt2[H2O] = 12.0*H2OHa*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d2fEnddt2[CO2] = 12.0*CO2Ha*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
        d3fEnddt3[H2O] = - 60.0*H2OHa*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
        d3fEnddt3[CO2] = - 60.0*CO2Ha*dCO2Trdt*dCO2Trdt*dCO2Trdt/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr/CO2Tr;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *fv = 0.0;
        fvPrime[H2O] = 0.0;
        fvPrime[CO2] = 0.0;

        f2vPrime2[H2O][H2O] = 0.0;
        f2vPrime2[H2O][CO2] = 0.0;
        f2vPrime2[CO2][CO2] = 0.0;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = 0.0;
        *d2fvdt2 = 0.0;
        *d3fvdt3 = 0.0;

        dfvPrimedt[H2O] = 0.0;
        dfvPrimedt[CO2] = 0.0;

        d2fvPrimedt2[H2O] = 0.0;
        d2fvPrimedt2[CO2] = 0.0;

        d3fvPrimedt3[H2O] = 0.0;
        d3fvPrimedt3[CO2] = 0.0;

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *fv = fEnd[H2O]*H2OVc*H2OVc;
        fvPrime[H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        fvPrime[CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = dfEnddt[H2O]*H2OVc*H2OVc;
        *d2fvdt2 = d2fEnddt2[H2O]*H2OVc*H2OVc;
        *d3fvdt3 = d3fEnddt3[H2O]*H2OVc*H2OVc;

        dfvPrimedt[H2O] = 2.0*dfEnddt[H2O]*H2OVc*H2OVc;
        dfvPrimedt[CO2] = 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        d2fvPrimedt2[H2O] = 2.0*d2fEnddt2[H2O]*H2OVc*H2OVc;
        d2fvPrimedt2[CO2] = 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

        d3fvPrimedt3[H2O] = 2.0*d3fEnddt3[H2O]*H2OVc*H2OVc;
        d3fvPrimedt3[CO2] = 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *fv = fEnd[CO2]*CO2Vc*CO2Vc;
        fvPrime[H2O] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        fvPrime[CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt   = dfEnddt[CO2]*CO2Vc*CO2Vc;
        *d2fvdt2 = d2fEnddt2[CO2]*CO2Vc*CO2Vc;
        *d3fvdt3 = d3fEnddt3[CO2]*CO2Vc*CO2Vc;

        dfvPrimedt[H2O] = 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        dfvPrimedt[CO2] = 2.0*dfEnddt[CO2]*CO2Vc*CO2Vc;

        d2fvPrimedt2[H2O] = 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        d2fvPrimedt2[CO2] = 2.0*d2fEnddt2[CO2]*CO2Vc*CO2Vc;

        d3fvPrimedt3[H2O] = 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        d3fvPrimedt3[CO2] = 2.0*d3fEnddt3[CO2]*CO2Vc*CO2Vc;

    } else {
        *fv  = 0.0;
        *fv += fEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *fv += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *fv += fEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        fvPrime[H2O]  = 0.0;
        fvPrime[H2O] += 2.0*fEnd[H2O]*H2OVc*H2OVc*x[H2O];
        fvPrime[H2O] += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        fvPrime[CO2]  = 0.0;
        fvPrime[CO2] += 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        fvPrime[CO2] += 2.0*fEnd[CO2]*CO2Vc*CO2Vc*x[CO2];

        f2vPrime2[H2O][H2O] = 2.0*fEnd[H2O]*H2OVc*H2OVc;
        f2vPrime2[H2O][CO2] = 2.0*powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0)*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0);
        f2vPrime2[CO2][CO2] = 2.0*fEnd[CO2]*CO2Vc*CO2Vc;
        f2vPrime2[CO2][H2O] = f2vPrime2[H2O][CO2];

        *dfvdt  = 0.0;
        *dfvdt += dfEnddt[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *dfvdt += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *dfvdt += dfEnddt[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        *d2fvdt2  = 0.0;
        *d2fvdt2 += d2fEnddt2[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *d2fvdt2 += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *d2fvdt2 += d2fEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        *d3fvdt3  = 0.0;
        *d3fvdt3 += d3fEnddt3[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        *d3fvdt3 += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        *d3fvdt3 += d3fEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        dfvPrimedt[H2O]  = 0.0;
        dfvPrimedt[H2O] += 2.0*dfEnddt[H2O]*H2OVc*H2OVc*x[H2O];
        dfvPrimedt[H2O] += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        dfvPrimedt[CO2]  = 0.0;
        dfvPrimedt[CO2] += 2.0*DpowSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        dfvPrimedt[CO2] += 2.0*dfEnddt[CO2]*CO2Vc*CO2Vc*x[CO2];

        d2fvPrimedt2[H2O]  = 0.0;
        d2fvPrimedt2[H2O] += 2.0*d2fEnddt2[H2O]*H2OVc*H2OVc*x[H2O];
        d2fvPrimedt2[H2O] += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        d2fvPrimedt2[CO2]  = 0.0;
        d2fvPrimedt2[CO2] += 2.0*D2powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        d2fvPrimedt2[CO2] += 2.0*d2fEnddt2[CO2]*CO2Vc*CO2Vc*x[CO2];

        d3fvPrimedt3[H2O]  = 0.0;
        d3fvPrimedt3[H2O] += 2.0*d3fEnddt3[H2O]*H2OVc*H2OVc*x[H2O];
        d3fvPrimedt3[H2O] += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[CO2];

        d3fvPrimedt3[CO2]  = 0.0;
        d3fvPrimedt3[CO2] += 2.0*D3powSum(fEnd[H2O], 1.0, fEnd[CO2], 1.0, dfEnddt[H2O], dfEnddt[CO2], d2fEnddt2[H2O], d2fEnddt2[CO2], d3fEnddt3[H2O], d3fEnddt3[CO2])*pow(powSum(H2OVc, 1.0, CO2Vc, 1.0), 2.0)*x[H2O];
        d3fvPrimedt3[CO2] += 2.0*d3fEnddt3[CO2]*CO2Vc*CO2Vc*x[CO2];

    }

    return;
}

static void BetaAndDerivative(int useLowPcoeff, double t, double x[2], double *beta, double betaPrime[2]) {
    double betaEnd[2];

    if (useLowPcoeff) {
        betaEnd[H2O] = H2OLb;
        betaEnd[CO2] = CO2Lb;
    } else {
        betaEnd[H2O] = H2OHb;
        betaEnd[CO2] = CO2Hb;
    }

    *beta  = betaEnd[H2O]*x[H2O] + betaEnd[CO2]*x[CO2];

    betaPrime[H2O]  = betaEnd[H2O];
    betaPrime[CO2]  = betaEnd[CO2];

    return;
}

static void GammaVcAndDerivative(int useLowPcoeff, double t, double x[2], double *gammav, double gammavPrime[2], double gamma2vPrime2[2][2]) {
    double gammaEnd[2];
    double k3 = 0.0;

    if (useLowPcoeff) {
        gammaEnd[H2O] = H2OLc;
        gammaEnd[CO2] = CO2Lc;
        k3 = 0.9;
    } else {
        gammaEnd[H2O] = H2OHc;
        gammaEnd[CO2] = CO2Hc;
        k3 = 1.0;
    }

    if ((x[H2O] == 0.0) && (x[CO2] == 0.0)) {
        *gammav = 0.0;
        gammavPrime[H2O] = 0.0;
        gammavPrime[CO2] = 0.0;

        gamma2vPrime2[H2O][H2O] = 0.0;
        gamma2vPrime2[H2O][CO2] = 0.0;
        gamma2vPrime2[CO2][CO2] = 0.0;
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else if ((x[H2O] == 1.0) && (x[CO2] == 0.0)) {
        *gammav = gammaEnd[H2O]*H2OVc*H2OVc;
        gammavPrime[H2O] = 3.0*gammaEnd[H2O]*H2OVc*H2OVc;
        gammavPrime[CO2] = 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);

        gamma2vPrime2[H2O][H2O] = 6.0*gammaEnd[H2O]*H2OVc*H2OVc;
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        gamma2vPrime2[CO2][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else if ((x[H2O] == 0.0) && (x[CO2] == 1.0)) {
        *gammav = gammaEnd[CO2]*CO2Vc*CO2Vc;
        gammavPrime[H2O] = 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gammavPrime[CO2] = 3.0*gammaEnd[CO2]*CO2Vc*CO2Vc;

        gamma2vPrime2[H2O][H2O] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0);
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0);
        gamma2vPrime2[CO2][CO2] = 6.0*gammaEnd[CO2]*CO2Vc*CO2Vc;
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    } else {
        *gammav  = 0.0;
        *gammav += gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O]*x[H2O];
        *gammav += 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O]*x[CO2];
        *gammav += 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2]*x[CO2];
        *gammav += gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2]*x[CO2];

        gammavPrime[H2O]  = 0.0;
        gammavPrime[H2O] += 3.0*gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]*x[H2O];
        gammavPrime[H2O] += 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[CO2];
        gammavPrime[H2O] += 3.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2]*x[CO2];

        gammavPrime[CO2]  = 0.0;
        gammavPrime[CO2] += 3.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]*x[H2O];
        gammavPrime[CO2] += 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]*x[CO2];
        gammavPrime[CO2] += 3.0*gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2]*x[CO2];

        gamma2vPrime2[H2O][H2O] = 6.0*gammaEnd[H2O]*H2OVc*H2OVc*x[H2O]
        + 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[CO2];
        gamma2vPrime2[H2O][CO2] = 6.0*powSum(gammaEnd[H2O], 2.0, gammaEnd[CO2], 1.0)*k3*pow(powSum(H2OVc, 2.0, CO2Vc, 1.0), 2.0)*x[H2O]
        + 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[CO2];
        gamma2vPrime2[CO2][CO2] = 6.0*powSum(gammaEnd[H2O], 1.0, gammaEnd[CO2], 2.0)*k3*pow(powSum(H2OVc, 1.0, CO2Vc, 2.0), 2.0)*x[H2O]
        + 6.0*gammaEnd[CO2]*CO2Vc*CO2Vc*x[CO2];
        gamma2vPrime2[CO2][H2O] = gamma2vPrime2[H2O][CO2];

    }

    return;
}

static void idealGasH2O(double t, double *cp, double *s0, double *h0, double *dcpdt) {
    int i;

    for (i=0, *cp=0.0; i<7; i++) *cp += idealCoeff[i][H2O]*pow(t/1000.0, (double) i);
    for (i=7; i<13; i++)         *cp += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-6));

    for (i=1, *dcpdt=0.0; i<7; i++) *dcpdt +=  ((double) i)  *idealCoeff[i][H2O]*pow(t/1000.0, (double) i-1);
    for (i=7; i<13; i++)            *dcpdt += -((double) i-6)*idealCoeff[i][H2O]/pow(t/1000.0, (double) (i+1-6));

    for (i=0, *h0=0.0; i<7; i++) *h0 += idealCoeff[i][H2O]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += idealCoeff[7][H2O]*log(t/1000.0);
    for (i=8; i<13; i++)         *h0 += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-7))/((double) (7-i));

    *s0  = idealCoeff[0][H2O]*log(t/1000.0);
    for (i=1; i<7; i++)	       *s0 += idealCoeff[i][H2O]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++)         *s0 += idealCoeff[i][H2O]/pow(t/1000.0, (double) (i-6))/((double) (6-i));

    *cp    *= 8.31451;
    *h0    *= 8.31451*1000.0;
    *s0    *= 8.31451;
    *dcpdt *= 8.31451/1000.0;

    *h0 += - 355665.4136;
    *s0 +=   359.6505;
}

static void duanH2ODriver(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                          double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {
    double bv, cv, dv, ev, fv, beta, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2;
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3;
    double bvPrime[2], cvPrime[2], dvPrime[2], evPrime[2], fvPrime[2], betaPrime[2], gammavPrime[2];
    double b2vPrime2[2][2], c2vPrime2[2][2], d2vPrime2[2][2], e2vPrime2[2][2], f2vPrime2[2][2], gamma2vPrime2[2][2];
    double dbvPrimedt[2], d2bvPrimedt2[2], d3bvPrimedt3[2], dcvPrimedt[2], d2cvPrimedt2[2], d3cvPrimedt3[2],
    ddvPrimedt[2], d2dvPrimedt2[2], d3dvPrimedt3[2], devPrimedt[2], d2evPrimedt2[2], d3evPrimedt3[2],
    dfvPrimedt[2], d2fvPrimedt2[2], d3fvPrimedt3[2];
    double lnPhiH2O, dlnPhiH2Odv, dlnPhiH2Odt, d2lnPhiH2Odv2, d2lnPhiH2Odt2, d2lnPhiH2Odvdt;
    double x[2] = { 1.0, 0.0 }; /* H2O */

    BVcAndDerivative    (useLowPcoeff, t, x, &bv,     bvPrime,    &dbvdt,    &d2bvdt2,  &d3bvdt3, dbvPrimedt, d2bvPrimedt2, d3bvPrimedt3, b2vPrime2);
    CVcAndDerivative    (useLowPcoeff, t, x, &cv,     cvPrime,    &dcvdt,    &d2cvdt2,  &d3cvdt3, dcvPrimedt, d2cvPrimedt2, d3cvPrimedt3, c2vPrime2);
    DVcAndDerivative    (useLowPcoeff, t, x, &dv,     dvPrime,	&ddvdt,    &d2dvdt2,  &d3dvdt3, ddvPrimedt, d2dvPrimedt2, d3dvPrimedt3, d2vPrime2);
    EVcAndDerivative    (useLowPcoeff, t, x, &ev,     evPrime,	&devdt,    &d2evdt2,  &d3evdt3, devPrimedt, d2evPrimedt2, d3evPrimedt3, e2vPrime2);
    FVcAndDerivative    (useLowPcoeff, t, x, &fv,     fvPrime,	&dfvdt,    &d2fvdt2,  &d3fvdt3, dfvPrimedt, d2fvPrimedt2, d3fvPrimedt3, f2vPrime2);
    BetaAndDerivative   (useLowPcoeff, t, x, &beta,   betaPrime);
    GammaVcAndDerivative(useLowPcoeff, t, x, &gammav, gammavPrime, gamma2vPrime2);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0, vLowest = 1.7, vHighest = 1.0e6; // 8/21/18 - MSG

        v = 8.314467*t/p;
        // if block added 8/21/18 - MSG
        if ((p > 1.0) && (p < 220.0) && (t > 273.15) && (t < 646.15)) {
            static const double psat_coeff[5] = { 1.44021565e+00, -2.75944904e-02,  3.50602876e-04, -2.44834016e-06,
                                                  1.57085668e-08 };
            static const double vsat_coeff[6] = { 1.84252342e+01, -3.06586710e-02,  5.65750627e-04, -2.69937313e-06,
                                                  4.67555414e-09, -8.89632469e+00 };
            double x = t - 273.15;
            double psat = psat_coeff[0] + psat_coeff[1]*x + psat_coeff[2]*x*x + psat_coeff[3]*x*x*x + psat_coeff[4]*x*x*x*x;
            if (p > psat) vHighest = (vsat_coeff[0] + vsat_coeff[1]*x + vsat_coeff[2]*x*x + vsat_coeff[3]*x*x*x
                + vsat_coeff[4]*x*x*x*x)/10.0 + vsat_coeff[5]/(x-374.0);
        }
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            if (v < 1.0) v = vLowest;
            if (v > vHighest) v = vHighest; // 8/21/18 - MSG
            iter++;
            if (iter == 200) {
                printf("duanH2ODriver(a): t = %g, p = %g, z = %g, v = %g, vPrev = %g, delv = %g, delvPrev = %g, iter = %d\n", t, p, z, v, vPrevious, delv, delvPrevious, iter);
                if (fabs(v-vLowest) <= 100.0*DBL_EPSILON) {
                    vLowest *= 0.90;
                    v = 8.314467*t/p;
                    iter = 0;
                }
            }
        }
        if (fabs(delv) > v*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                v = rtb + (dx *= 0.5);
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) printf("duanH2ODriver(b): t = %g, p = %g, z = %g, v = %g, delv = %g, iter = %d\n", t, p, z, v, delv, iter);
        }
    }

    lnPhiH2O  = 0.0;
    lnPhiH2O += -log(z);
    lnPhiH2O += bvPrime[H2O]/v;
    lnPhiH2O += cvPrime[H2O]/2.0/v/v;
    lnPhiH2O += dvPrime[H2O]/4.0/v/v/v/v;
    lnPhiH2O += evPrime[H2O]/5.0/v/v/v/v/v;
    lnPhiH2O += ((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiH2O += ((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)
    *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiH2O += ((gammavPrime[H2O]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    dzdv = - bv/v/v - 2.0*cv/v/v/v + - 4.0*dv/v/v/v/v/v - 5.0*ev/v/v/v/v/v/v
    - 2.0*(fv/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(fv/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(fv/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    dzdt = dbvdt/v + dcvdt/v/v + ddvdt/v/v/v/v + devdt/v/v/v/v/v + (dfvdt/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    dlnPhiH2Odv  = 0.0;
    dlnPhiH2Odv += -dzdv/z;
    dlnPhiH2Odv += -bvPrime[H2O]/v/v;
    dlnPhiH2Odv += -cvPrime[H2O]/v/v/v;
    dlnPhiH2Odv += -dvPrime[H2O]/v/v/v/v/v;
    dlnPhiH2Odv += -evPrime[H2O]/v/v/v/v/v/v;
    dlnPhiH2Odv += -((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += -((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += +((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiH2Odt  = 0.0;
    dlnPhiH2Odt += -dzdt/z;
    dlnPhiH2Odt += dbvPrimedt[H2O]/v;
    dlnPhiH2Odt += dcvPrimedt[H2O]/2.0/v/v;
    dlnPhiH2Odt += ddvPrimedt[H2O]/4.0/v/v/v/v;
    dlnPhiH2Odt += devPrimedt[H2O]/5.0/v/v/v/v/v;
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiH2Odt += ((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    dlnPhiH2Odt += ((gammavPrime[H2O]-gammav)*dfvdt/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    d2zdv2 =  2.0*bv/v/v/v + 6.0*cv/v/v/v/v + 20.0*dv/v/v/v/v/v/v + 30.0*ev/v/v/v/v/v/v/v
    +  6.0*(fv/v/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    + 14.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    - 14.0*(fv/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    -  8.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    +  4.0*(fv/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
    d2zdvdt = - dbvdt/v/v - 2.0*dcvdt/v/v/v + - 4.0*ddvdt/v/v/v/v/v - 5.0*devdt/v/v/v/v/v/v
    - 2.0*(dfvdt/v/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v)
    - 2.0*(dfvdt/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
    + 2.0*(dfvdt/v/v) * (beta + gammav/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    d2zdt2 = d2bvdt2/v + d2cvdt2/v/v + d2dvdt2/v/v/v/v + d2evdt2/v/v/v/v/v + (d2fvdt2/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);

    d2lnPhiH2Odv2  = 0.0;
    d2lnPhiH2Odv2 += dzdv*dzdv/z/z - d2zdv2/z;
    d2lnPhiH2Odv2 +=  2.0*bvPrime[H2O]/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*cvPrime[H2O]/v/v/v/v;
    d2lnPhiH2Odv2 +=  5.0*dvPrime[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  6.0*evPrime[H2O]/v/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*((fvPrime[H2O]*beta + betaPrime[H2O]*fv)/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*beta*gammav + betaPrime[H2O]*fv*gammav)/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  5.0*((fvPrime[H2O]*gammav+gammavPrime[H2O]*fv-fv*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -2.0*((fvPrime[H2O]*gammav*gammav+gammavPrime[H2O]*fv*gammav-fv*beta*gammav*(gammavPrime[H2O]-gammav))/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 += -7.0*((gammavPrime[H2O]-gammav)*fv*gammav/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  2.0*((gammavPrime[H2O]-gammav)*fv*gammav*gammav/v/v/v/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odvdt  = 0.0;
    d2lnPhiH2Odvdt += dzdv*dzdt/z/z -d2zdvdt/z;
    d2lnPhiH2Odvdt += -dbvPrimedt[H2O]/v/v;
    d2lnPhiH2Odvdt += -dcvPrimedt[H2O]/v/v/v;
    d2lnPhiH2Odvdt += -ddvPrimedt[H2O]/v/v/v/v/v;
    d2lnPhiH2Odvdt += -devPrimedt[H2O]/v/v/v/v/v/v;
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*beta + betaPrime[H2O]*dfvdt)/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt += -((dfvPrimedt[H2O]*gammav+gammavPrime[H2O]*dfvdt-dfvdt*beta*(gammavPrime[H2O]-gammav))/v/v/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt +=  ((gammavPrime[H2O]-gammav)*dfvdt*gammav/v/v/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odt2  = 0.0;
    d2lnPhiH2Odt2 += dzdt*dzdt/z/z - d2zdt2/z;
    d2lnPhiH2Odt2 += d2bvPrimedt2[H2O]/v;
    d2lnPhiH2Odt2 += d2cvPrimedt2[H2O]/2.0/v/v;
    d2lnPhiH2Odt2 += d2dvPrimedt2[H2O]/4.0/v/v/v/v;
    d2lnPhiH2Odt2 += d2evPrimedt2[H2O]/5.0/v/v/v/v/v;
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*beta + betaPrime[H2O]*d2fvdt2)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((d2fvPrimedt2[H2O]*gammav+gammavPrime[H2O]*d2fvdt2-d2fvdt2*beta*(gammavPrime[H2O]-gammav))/2.0/gammav/gammav)*(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    d2lnPhiH2Odt2 += ((gammavPrime[H2O]-gammav)*d2fvdt2/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    *vPt = v;
    *zPt = z;
    *phi = exp(lnPhiH2O);
    *dvdp = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2  = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2 *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2  = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    *dlnphidt    = dlnPhiH2Odv*(*dvdt) + dlnPhiH2Odt;
    *dlnphidp    = dlnPhiH2Odv*(*dvdp);
    *d2lnphidt2  = d2lnPhiH2Odv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiH2Odvdt*(*dvdt) + dlnPhiH2Odv*(*d2vdt2) + d2lnPhiH2Odt2;
    *d2lnphidtdp = d2lnPhiH2Odv2*(*dvdt)*(*dvdp) + dlnPhiH2Odv*(*d2vdtdp) + d2lnPhiH2Odvdt*(*dvdp);
    *d2lnphidp2  = d2lnPhiH2Odv2*(*dvdp)*(*dvdp) + dlnPhiH2Odv*(*d2vdp2);
}

static void duanH2O(int useLowPcoeff, double t, double p, double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                    double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {

    duanH2ODriver(useLowPcoeff, t, p, vPt, zPt, phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2);

    if (!useLowPcoeff) {
        double vRef, zRef, phiRef, dvdtRef, dvdpRef, d2vdt2Ref, d2vdtdpRef, d2vdp2Ref, dlnphidtRef, dlnphidpRef, d2lnphidt2Ref, d2lnphidtdpRef, d2lnphidp2Ref;

        duanH2ODriver(1, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         *= phiRef;
        *dlnphidt    += dlnphidtRef;
        *d2lnphidt2  += d2lnphidt2Ref;

        duanH2ODriver(0, t, 2000.0, &vRef, &zRef, &phiRef, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, &dlnphidtRef, &dlnphidpRef, &d2lnphidt2Ref, &d2lnphidtdpRef, &d2lnphidp2Ref);
        *phi         /= phiRef;
        *dlnphidt    -= dlnphidtRef;
        *d2lnphidt2  -= d2lnphidt2Ref;
    }

}

static void propertiesOfPureH2O(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt,
                                double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    double z, phi, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2;

    idealGasH2O(t, cp, s, h, dcpdt);
    duanH2O((p <= 2000.0) ? 1 : 0, t, p, v, &z, &phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, &dlnphidt, &dlnphidp, &d2lnphidt2, &d2lnphidtdp, &d2lnphidp2);

    *g      = *h - t*(*s) + R*t*log(phi*p);
    *s     += - (R*log(phi*p) + R*t*dlnphidt);
    *h     += R*t*log(phi*p) - t*(R*log(phi*p) + R*t*dlnphidt);
    *cp    += -t*(2.0*R*dlnphidt + R*t*d2lnphidt2);
    {
        double zTemp, phiTemp, dlnphidtTemp, dlnphidpTemp, d2lnphidt2Temp, d2lnphidtdpTemp, d2lnphidp2Temp, vTemp, dvdtTemp, dvdpTemp, d2vdt2Temp, d2vdtdpTemp, d2vdp2Temp, d3lnphidt3;

        duanH2O((p <= 2000.0) ? 1 : 0, t*(1.0+sqrt(DBL_EPSILON)), p, &vTemp, &zTemp, &phiTemp, &dvdtTemp, &dvdpTemp, &d2vdt2Temp, &d2vdtdpTemp, &d2vdp2Temp,
                &dlnphidtTemp, &dlnphidpTemp, &d2lnphidt2Temp, &d2lnphidtdpTemp, &d2lnphidp2Temp);

        d3lnphidt3 = (d2lnphidt2Temp - d2lnphidt2)/t/sqrt(DBL_EPSILON);
        *dcpdt += -(2.0*R*dlnphidt + R*t*d2lnphidt2) -t*(3.0*R*d2lnphidt2 + R*t*d3lnphidt3);
    }
}

static const unsigned int returnValueOfG       =  1;
static const unsigned int returnValueOfH       =  2;
static const unsigned int returnValueOfS       =  3;
static const unsigned int returnValueOfCP      =  4;
static const unsigned int returnValueOfDCPDT   =  5;
static const unsigned int returnValueOfV       =  6;
static const unsigned int returnValueOfdVdT    =  7;
static const unsigned int returnValueOfdVdP    =  8;
static const unsigned int returnValueOfd2VdT2  =  9;
static const unsigned int returnValueOfd2VdTdP = 10;
static const unsigned int returnValueOfd2VdP2  = 11;

static double calculate(double t, double p, int returnMode) {
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
	double result = 0.0;

    propertiesOfPureH2O(t, p, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

	switch (returnMode) {
		case 1:
			result = gH2O;
			break;
		case 2:
			result = hH2O;
			break;
		case 3:
			result = sH2O;
			break;
		case 4:
			result = cpH2O;
			break;
		case 5:
			result = dcpdtH2O;
			break;
		case 6:
			result = vH2O;
			break;
		case 7:
			result = dvdtH2O;
			break;
		case 8:
			result = dvdpH2O;
			break;
		case 9:
			result = d2vdt2H2O;
			break;
		case 10:
			result = d2vdtdpH2O;
			break;
		case 11:
			result = d2vdp2H2O;
			break;
		default:
			break;
	}
	return result;
}

double DuanAndZhang2006_getGibbsFreeEnergy(double t, double p) {
	return calculate(t, p, returnValueOfG);
}

double DuanAndZhang2006_getEnthalpy(double t, double p) {
	return calculate(t, p, returnValueOfH);
}

double DuanAndZhang2006_getEntropy(double t, double p) {
	return calculate(t, p, returnValueOfS);
}

double DuanAndZhang2006_getHeatCapacity(double t, double p) {
	return calculate(t, p, returnValueOfCP);
}

double DuanAndZhang2006_getDcpDt(double t, double p) {
	return calculate(t, p, returnValueOfDCPDT);
}

double DuanAndZhang2006_getVolume(double t, double p) {
	return  calculate(t, p, returnValueOfV);
}

double DuanAndZhang2006_getDvDt(double t, double p) {
	return calculate(t, p, returnValueOfdVdT);
}

double DuanAndZhang2006_getDvDp(double t, double p) {
	return  calculate(t, p, returnValueOfdVdP);
}

double DuanAndZhang2006_getD2vDt2(double t, double p) {
	return calculate(t, p, returnValueOfd2VdT2);
}

double DuanAndZhang2006_getD2vDtDp(double t, double p) {
	return calculate(t, p, returnValueOfd2VdTdP);
}

double DuanAndZhang2006_getD2vDp2(double t, double p) {
	return calculate(t, p, returnValueOfd2VdP2);
}
