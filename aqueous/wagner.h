#ifndef INCLUDE_WAGNER
#define INCLUDE_WAGNER

double WagnerEtAl2002_getGibbsFreeEnergy(double t, double p);
double WagnerEtAl2002_getEnthalpy(double t, double p);
double WagnerEtAl2002_getEntropy(double t, double p);
double WagnerEtAl2002_getHeatCapacity(double t, double p);
double WagnerEtAl2002_getDcpDt(double t, double p);
double WagnerEtAl2002_getVolume(double t, double p);
double WagnerEtAl2002_getDvDt(double t, double p);
double WagnerEtAl2002_getDvDp(double t, double p);
double WagnerEtAl2002_getD2vDt2(double t, double p);
double WagnerEtAl2002_getD2vDtDp(double t, double p);
double WagnerEtAl2002_getD2vDp2(double t, double p);

#endif
