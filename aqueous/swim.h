#ifndef INCLUDE_SWIM
#define INCLUDE_SWIM

typedef enum {
    NotApplicable,
    DuanAndZhang2006,
    ZhangAndDuan2005,
    HoltenEtAl2014,
    WagnerEtAl2002
} SWIM_RegionType;

double SWIM_getGibbsFreeEnergy(double t, double p, SWIM_RegionType region);
double SWIM_getEnthalpy(double t, double p, SWIM_RegionType region);
double SWIM_getEntropy(double t, double p, SWIM_RegionType region);
double SWIM_getHeatCapacity(double t, double p, SWIM_RegionType region);
double SWIM_getDcpDt(double t, double p, SWIM_RegionType region);
double SWIM_getVolume(double t, double p, SWIM_RegionType region);
double SWIM_getDvDt(double t, double p, SWIM_RegionType region);
double SWIM_getDvDp(double t, double p, SWIM_RegionType region);
double SWIM_getD2vDt2(double t, double p, SWIM_RegionType region);
double SWIM_getD2vDtDp(double t, double p, SWIM_RegionType region);
double SWIM_getD2vDp2(double t, double p, SWIM_RegionType region);

const char *eosForRegion(double t, double p);
double lowerTemperatureLimitAtPinBars(double p);

double tTransHoltenWagner(void);
double tTransWagnerDZ2006(void);
double tTransHoltenZD2005(void);
double pTransWagnerZD2005(void);
double pTransDZ2006ZD2005(void);

double SWIM_getRho(double t, double p);
double SWIM_getDrhoDt(double t, double p);
double SWIM_getDrhoDp(double t, double p);
double SWIM_getD2rhoDt2(double t, double p);
double SWIM_getD2rhoDtDp(double t, double p);
double SWIM_getD2rhoDp2(double t, double p);

#endif
