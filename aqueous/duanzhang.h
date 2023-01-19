#ifndef INCLUDE_DUANZHANG
#define INCLUDE_DUANZHANG

extern double DuanAndZhang2006_getGibbsFreeEnergy(double t, double p);
extern double DuanAndZhang2006_getEnthalpy(double t, double p);
extern double DuanAndZhang2006_getEntropy(double t, double p);
extern double DuanAndZhang2006_getHeatCapacity(double t, double p);
extern double DuanAndZhang2006_getDcpDt(double t, double p);
extern double DuanAndZhang2006_getVolume(double t, double p);
extern double DuanAndZhang2006_getDvDt(double t, double p);
extern double DuanAndZhang2006_getDvDp(double t, double p);
extern double DuanAndZhang2006_getD2vDt2(double t, double p);
extern double DuanAndZhang2006_getD2vDtDp(double t, double p);
extern double DuanAndZhang2006_getD2vDp2(double t, double p);

#endif
