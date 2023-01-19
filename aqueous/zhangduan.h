#ifndef INCLUDE_ZHANGDUAN
#define INCLUDE_ZHANGDUAN

double ZhangAndDuan2005_getGibbsFreeEnergy(double t, double p);
double ZhangAndDuan2005_getEnthalpy(double t, double p);
double ZhangAndDuan2005_getEntropy(double t, double p);
double ZhangAndDuan2005_getHeatCapacity(double t, double p);
double ZhangAndDuan2005_getDcpDt(double t, double p);
double ZhangAndDuan2005_getVolume(double t, double p);
double ZhangAndDuan2005_getDvDt(double t, double p);
double ZhangAndDuan2005_getDvDp(double t, double p);
double ZhangAndDuan2005_getD2vDt2(double t, double p);
double ZhangAndDuan2005_getD2vDtDp(double t, double p);
double ZhangAndDuan2005_getD2vDp2(double t, double p);

#endif
