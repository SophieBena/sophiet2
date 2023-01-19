#include "holten.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
  holten.c
  ThermoFit

	Program implementing the equation of state in the article
	"Equation of state for supercooled water at pressures up to 400 MPa",
	Journal of Physical and Chemical Reference Data,
	Vincent Holten, Jan V. Sengers, and Mikhail A. Anisimov,
	Institute for Physical Science and Technology and
	Department of Chemical and Biomolecular Engineering,
	University of Maryland, College Park, Maryland 20742, U.S.A.

*/

static const double Tc     =  228.2;
static const double Pc     =    0;
static const double rho0   = 1081.6482;
static const double R      =  461.523087;
static const double omega0 =    0.52122690;
static const double L0     =    0.76317954;
static const double k0     =    0.072158686;
static const double k1     =  - 0.31569232;
static const double k2     =    5.2992608;

static const double MolesPerKg = 55.508435;

static const int n = 20;

static const double c[n] = {-8.1570681381655, 1.2875032e+000, 7.0901673598012,
    -3.2779161e-002, 7.3703949e-001, -2.1628622e-001, -5.1782479e+000,
    4.2293517e-004, 2.3592109e-002, 4.3773754e+000, -2.9967770e-003,
    -9.6558018e-001, 3.7595286e+000, 1.2632441e+000, 2.8542697e-001,
    -8.5994947e-001, -3.2916153e-001, 9.0019616e-002, 8.1149726e-002,
    -3.2788213e+000};
static const double a[n] = {0, 0, 1, -0.2555, 1.5762, 1.64, 3.6385, -0.3828,
    1.6219, 4.3287, 3.4763, 5.1556, -0.3593, 5.0361, 2.9786, 6.2373,
    4.046, 5.3558, 9.0157, 1.2194};
static const double b[n] = {0, 1, 0, 2.1051, 1.1422, 0.951, 0, 3.6402,
    2.076, -0.0016, 2.2769, 0.0008, 0.3706, -0.3975, 2.973, -0.318,
    2.9805, 2.9265, 0.4456, 0.1298};
static const double d[n] = {0, 0, 0, -0.0016, 0.6894, 0.013, 0.0002, 0.0435,
    0.05, 0.0004, 0.0528, 0.0147, 0.8584, 0.9924, 1.0041, 1.0961,
    1.0228, 1.0303, 1.618, 0.5213};

static double B(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

static double Bp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau, a[i]) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) * exp(-d[i]*pi);
    return sum;
}

static double Bt(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

static double Bpp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * pow(tau,a[i]) * pow(pi, b[i]-2) * (pow(d[i]*pi - b[i], 2) - b[i]) * exp(-d[i]*pi);
    return sum;
}

static double Btp(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * pow(tau,a[i]-1) * pow(pi, b[i]-1) * (b[i]-d[i]*pi) *  exp(-d[i]*pi);
    return sum;
}

static double Btt(double tau, double pi) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += c[i] * a[i] * (a[i]-1) * pow(tau,a[i]-2) * pow(pi, b[i]) * exp(-d[i]*pi);
    return sum;
}

static double xefun(double x, double L, double W) {
    return L + log(x/(1-x)) + W*(1-2*x);
}

typedef enum { false, true } bool;

static double findxe(double L, double W) {
    double x0, x1;
    bool flip = false;

    if (L < 0) {
        flip = true;
        L = -L;
    }

    /* Find starting values for x */
    if (W < 1.1111111*(2.944439 - L)) { /* xe = 0.05 isoline, W = (10/9) * (ln(19) - L) */
        x0 = 0.049;
        x1 = 0.5;
    } else if (W < 1.0204081*(4.595119 - L)) { /* xe = 0.01 isoline, W = (50/49) * (ln(99) - L) */
        x0 = 0.0099;
        x1 = 0.051;
    } else {
        x0 = 0.99 * exp(-1.0204081 * L - W);
        x1 = 1.01 * 1.087 * exp(-L - W);
        if (x1 > 0.0101) x1 = 0.0101;
    }

    double y0 = xefun(x0, L, W),
    y1 = xefun(x1, L, W);
    if (y0*y1 < 0) printf("Error in findxe(): starting values for x are incorrect.");

    /*
       Bisection algorithm
       This could be replaced with a root-finding function from a numerical library
    */
    int N = 0;
    while (fabs(x1 - x0) > 10.0*DBL_EPSILON) {
        double x = (x0 + x1) / 2.0;
        double y = xefun(x, L, W);
        if (y0 * y >= 0) {
            x0 = x;
            y0 = y;
        } else {
            x1 = x;
        }
        if (N++ < 50) printf("Error in findxe(): bisection does not converge.");
    }

    double x = (x0 + x1) / 2.0;
    return (flip) ? 1.0 - x : x;
}

typedef enum { ValueOfG, ValueOfS, ValueOfRho, ValueOfKap, ValueOfAlp, ValueOfCp, ValueOfCv, ValueOfU } ReturnValue;

static double evalAtTinK(double T, double P, ReturnValue returnValue) {
    static double tOld = -1.0, pOld = -1.0;
    static double G, S, rho, Kap, Alp, CP, CV, U;
    static double EPS = -1.0, EPS2 = -1.0;
    if (EPS == -1.0) {
        EPS = sqrt(DBL_EPSILON);
        EPS2 = sqrt(sqrt(DBL_EPSILON));
    }
    if ((fabs(T - tOld) < 10.0*DBL_EPSILON) && (fabs(P - pOld) < 10.0*DBL_EPSILON)) {
        if      (returnValue == ValueOfG)       return G;
        else if (returnValue == ValueOfS)       return S;
        else if (returnValue == ValueOfRho)     return rho;
        else if (returnValue == ValueOfKap)     return Kap;
        else if (returnValue == ValueOfAlp)     return Alp;
        else if (returnValue == ValueOfCp)      return CP;
        else if (returnValue == ValueOfCv)      return CV;
        else if (returnValue == ValueOfU)       return U;
        else return 0.0;
    }

    const double P0 = -300e6;

    /* Dimensionless temperature and pressure */
    double t   =  (T - Tc)/Tc;
    double p   =  (P - Pc)/(rho0*R*Tc);
    double tau = T/Tc;
    double pi  = (P - P0)/(rho0*R*Tc);

    /* Field L and its derivatives */
    double K1  = sqrt(pow(1.0+k0*k2+k1*(p-k2*t), 2.0) - 4.0*k0*k1*k2*(p-k2*t));
    double K3  = pow(K1, 3.0);
    double K2  = sqrt(1.0 + k2*k2);
    double L   = L0 * K2 * (1.0 - K1 + k0*k2 + k1*(p + k2*t)) / (2.0*k1*k2);
    double Lt  = L0 * 0.5 * K2 * (1 + (1 - k0*k2 + k1*(p - k2*t))/K1);
    double Lp  = L0 * K2 * (K1 + k0*k2 - k1*p + k1*k2*t - 1) / (2.0*k2*K1);
    double Ltt = -2.0*L0*K2*k0*k1*k2*k2 / K3;
    double Ltp =  2.0*L0*K2*k0*k1*k2 / K3;
    double Lpp = -2.0*L0*K2*k0*k1 / K3;

    /* Interaction parameter omega */
    double omega = 2.0 + omega0 * p;

    /* Calculate equilibrium fraction xe */
    double x = findxe(L, omega);

    /* Order parameter f and susceptibility chi */
    double f   = 2.0 * x - 1.0;
    double f2  = f*f;
    double chi = 1.0 / (2.0 / (1.0 - f2) - omega);

    /* Dimensionless properties */
    double g0  = x*L + x*log(x) + (1.0-x)*log(1.0-x) + omega*x*(1-x); /* g0 = (g - gA)/tau */
    double g   = B(tau,pi) + tau*g0;
    double s   = -0.5*(f+1)*Lt*tau - g0 - Bt(tau,pi);
    double v   = 0.5 * tau * (omega0/2.0*(1.0-f2) + Lp*(f+1)) + Bp(tau,pi);
    double kap = (1.0/v) * (tau/2.0 * (chi * pow(Lp - omega0*f, 2.0) - (f+1.0)*Lpp) - Bpp(tau,pi));
    double alp = (1.0/v) * (Ltp/2.0 * tau*(f+1.0) + (omega0/2.0*(1.0-f2) + Lp*(f+1.0))/2.0 - tau*Lt/2.0 * chi*(Lp - omega0*f) + Btp(tau,pi));
    double cp  = tau * ( -Lt * (f+1) + tau*(Lt*Lt*chi - Ltt*(f+1)) / 2 - Btt(tau,pi));

    /* Properties in SI units */
    G   = R * Tc * g;					  /* Specific Gibbs energy       J/kg   */
    S   = R * s;						  /* Specific entropy            J/kg-K */
    rho = rho0 / v;						  /* Density                     kg/m^3 */
    Kap = kap / (rho0 * R * Tc);		  /* Isothermal compressibility  Pa     */
    Alp = alp / Tc;						  /* Expansion coefficient       1/K    */
    CP  = R * cp;						  /* Isobaric heat capacity      J/kg-K */
    CV  = CP - T*Alp*Alp / (rho * Kap);   /* Isochoric heat capacity     J/kg-K */
    U   = 1/sqrt(rho*Kap - T*Alp*Alp/CP); /* Speed of sound              m/s    */

    tOld = T;
    pOld = P;

    if      (returnValue == ValueOfG)       return G;
    else if (returnValue == ValueOfS)       return S;
    else if (returnValue == ValueOfRho)     return rho;
    else if (returnValue == ValueOfKap)     return Kap;
    else if (returnValue == ValueOfAlp)     return Alp;
    else if (returnValue == ValueOfCp)      return CP;
    else if (returnValue == ValueOfCv)      return CV;
    else if (returnValue == ValueOfU)       return U;
    else return 0.0;
}

double HoltenEtAl2014_homogeneousIceNucleationTemperatureForPressureInBars(double p) {
    double p1 = p/10.0; // MPa
    double result = 172.82 + 0.03718*p1 + 3.403e-5*p1*p1 - 1.573e-8*p1*p1*p1;
    return (result > 0.0) ? result : 1.0;
}

double HoltenEtAl2014_homogeneousIceNucleationPressureForTemperatureInK(double t) {
    double p0 = 0.1;    // MPa
    double t0 = 235.15; // K
    double theta = t/t0;
    double result = 1.0 + 2282.7*(1.0-pow(theta, 6.243)) + 157.24*(1.0-pow(theta, 79.81));
    return result*p0*10.0; //bars
}

double HoltenEtAl2014_getGibbsFreeEnergy(double t, double p) {
    double G = evalAtTinK(t, p*1.0e5, ValueOfG);
    return G/MolesPerKg; /* J/mol */
}

double HoltenEtAl2014_getEnthalpy(double t, double p) {
    double G = evalAtTinK(t, p*1.0e5, ValueOfG);
    double S = evalAtTinK(t, p*1.0e5, ValueOfS);
    return (G+t*S)/MolesPerKg; // J/mol
}

double HoltenEtAl2014_getEntropy(double t, double p) {
    double S = evalAtTinK(t, p*1.0e5, ValueOfS);
    return S/MolesPerKg; // J/mol-K
}

double HoltenEtAl2014_getHeatCapacity(double t, double p) {
    double CP = evalAtTinK(t, p*1.0e5, ValueOfCp);
    return CP/MolesPerKg; // J/mol-K
}

double HoltenEtAl2014_getDcpDt(double t, double p) {
    static double EPS = -1.0;
    if (EPS == -1.0) EPS = sqrt(DBL_EPSILON);
    double CP = evalAtTinK(t*(1.0+EPS), p*1.0e5, ValueOfCp);
    double result = CP;
    CP = evalAtTinK(t*(1.0-EPS), p*1.0e5, ValueOfCp);
    result -= CP;
    return result/2.0/t/EPS/MolesPerKg;
}

double HoltenEtAl2014_getVolume(double t, double p) {
    double rho = evalAtTinK(t, p*1.0e5, ValueOfRho);
    return 1.0e5/rho/MolesPerKg; /* J/mol-bar */
}

double HoltenEtAl2014_getDvDt(double t, double p) {
    double rho = evalAtTinK(t, p*1.0e5, ValueOfRho);
    double Alp = evalAtTinK(t, p*1.0e5, ValueOfAlp);
    return Alp*1.0e5/rho/MolesPerKg; /* J/mol-bar-K */
}

double HoltenEtAl2014_getDvDp(double t, double p) {
    double rho = evalAtTinK(t, p*1.0e5, ValueOfRho);
    double Kap = evalAtTinK(t, p*1.0e5, ValueOfKap);
    return  -1.0e10*Kap/rho/MolesPerKg; // J/mol-bar^2
}

double HoltenEtAl2014_getD2vDt2(double t, double p) {
    static double EPS2 = -1.0;
    if (EPS2 == -1.0) EPS2 = sqrt(sqrt(DBL_EPSILON));
    double rho = evalAtTinK(t*(1.0+EPS2), p*1.0e5, ValueOfRho);
    double Alp = evalAtTinK(t*(1.0+EPS2), p*1.0e5, ValueOfAlp);
    double result = Alp/rho;
    rho = evalAtTinK(t*(1.0-EPS2), p*1.0e5, ValueOfRho);
    Alp = evalAtTinK(t*(1.0-EPS2), p*1.0e5, ValueOfAlp);
    result -= Alp/rho;
    return result*1.0e5/2.0/t/EPS2/MolesPerKg;
}

double HoltenEtAl2014_getD2vDtDp(double t, double p) {
    static double EPS2 = -1.0;
    if (EPS2 == -1.0) EPS2 = sqrt(sqrt(DBL_EPSILON));
    double rho = evalAtTinK(t, p*1.0e5*(1.0-2.0*EPS2), ValueOfRho);
    double Alp = evalAtTinK(t, p*1.0e5*(1.0-2.0*EPS2), ValueOfAlp);
    double result = Alp/rho/12.0;
    rho = evalAtTinK(t, p*1.0e5*(1.0-EPS2), ValueOfRho);
    Alp = evalAtTinK(t, p*1.0e5*(1.0-EPS2), ValueOfAlp);
    result -= 2.0*Alp/rho/3.0;
    rho = evalAtTinK(t, p*1.0e5*(1.0+EPS2), ValueOfRho);
    Alp = evalAtTinK(t, p*1.0e5*(1.0+EPS2), ValueOfAlp);
    result += 2.0*Alp/rho/3.0;
    rho = evalAtTinK(t, p*1.0e5*(1.0+2.0*EPS2), ValueOfRho);
    Alp = evalAtTinK(t, p*1.0e5*(1.0+2.0*EPS2), ValueOfAlp);
    result -= Alp/rho/12.0;
    return result*1.0e5/p/EPS2/MolesPerKg; // *1.0e5
}

double HoltenEtAl2014_getD2vDp2(double t, double p) {
    static double EPS2 = -1.0;
    if (EPS2 == -1.0) EPS2 = sqrt(sqrt(DBL_EPSILON));
    double rho = evalAtTinK(t, p*1.0e5*(1.0+EPS2), ValueOfRho);
    double Kap = evalAtTinK(t, p*1.0e5*(1.0+EPS2), ValueOfKap);
    double result = Kap/rho;
    rho = evalAtTinK(t, p*1.0e5*(1.0-EPS2), ValueOfRho);
    Kap = evalAtTinK(t, p*1.0e5*(1.0-EPS2), ValueOfKap);
    result -= Kap/rho;
    return -result*1.0e10/2.0/p/EPS2/MolesPerKg;
}
