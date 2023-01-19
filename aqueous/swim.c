
#include "swim.h"
#include "wagner.h"
#include "zhangduan.h"
#include "duanzhang.h"
#include "holten.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double widthOfWagnerToDZ2006SmoothInterval  =  50.0;  // K
static const double centerOfWagnerToDZ2006SmoothInterval = 673.15; // K
static const double widthOfWagnerToZD2005SmoothInterval  =  100.0; // bars
static const double centerOfWagnerToZD2005SmoothInterval = 1000.0; // bars
static const double widthOfDZ2006ToZD2005SmoothInterval  =  100.0; // bars
static const double lowTofZD2005AdjustmentRegion         = 298.15; // K
static const double highTofZD2005AdjustmentRegion        = 398.15; // K
static const double centerOfDZ2006ToZD2005SmoothInterval = 1000.0; // bars

static const double HliqReference = -285830.0;      // J/mol, 298.15 K, 1 bar  NIST, CODATA
static const double SliqReference =  69.950;        // J/K-mol, 298.15, 1 bar  NIST, CODATA
static const double GliqReference = -56687.0*4.184; // J/mol, 298.15 K, 1 bar Helgeson and Kirkham, 1974

typedef enum { valueOfG, valueOfH, valueOfS, valueOfCP, valueOfDCPDT, valueOfV, valueOfdVdT,
                     valueOfdVdP, valueOfd2VdT2, valueOfd2VdTdP, valueOfd2VdP2 } SWIM_Property;

static const char *region0 = "No EOS is applicable at the specified temperature and pressure";
static const char *region1 = "MELTS H2O-CO2 from Duan and Zhang (2006)";
static const char *region2 = "DEW H2O from Zhang and Duan (2005)";
static const char *region3 = "Supercooled H2O from Holten et al. (2014)";
static const char *region4 = "Steam Properties from Wagner et al. (2002)";

SWIM_RegionType classifyRegion(double t, double p) {
    SWIM_RegionType region = NotApplicable;
    if      ((t > 673.15) && (p < 1000.0)) region = DuanAndZhang2006;
    else if ((t > 298.15) && (p > 1000.0)) region = ZhangAndDuan2005;
    else if  (t < 298.15)                  region = HoltenEtAl2014;
    else if  (t < 673.15)                  region = WagnerEtAl2002;
    else                                   region = ZhangAndDuan2005;
    return region;
}

double tTransHoltenWagner(void) { return 298.15; }
double tTransWagnerDZ2006(void) { return 673.15; }
double tTransHoltenZD2005(void) { return 298.15; }
double pTransWagnerZD2005(void) { return 1000.0; }
double pTransDZ2006ZD2005(void) { return 1000.0; }

const char *eosForRegion(double t, double p) {
    SWIM_RegionType region = classifyRegion(t, p);
    const char *result = region0;
    switch (region) {
        case DuanAndZhang2006:
            result = region1;
            break;
        case ZhangAndDuan2005:
            result = region2;
            break;
        case HoltenEtAl2014:
            result = region3;
            break;
        case WagnerEtAl2002:
            result = region4;
            break;
        default:
            break;
    }
    return result;
}

double lowerTemperatureLimitAtPinBars(double p) {
    double TK = 0.0;
    if (p >= 2000.0) TK = HoltenEtAl2014_homogeneousIceNucleationTemperatureForPressureInBars(p);
    else {
        double TKold = 200.0;
        double pOld = HoltenEtAl2014_homogeneousIceNucleationPressureForTemperatureInK(TKold);
        TK = 220.0;
        double pGuess = HoltenEtAl2014_homogeneousIceNucleationPressureForTemperatureInK(TK);
        int iter = 0;
        while ((fabs(pGuess-p) > 0.1) && (iter < 200)) {
            double deriv = (pGuess - pOld)/(TK-TKold);
            TKold = TK;
            pOld = pGuess;
            TK = TK + (p-pGuess)/deriv;
            if (TK < 175.0 ) TK = 175.0;
            if (TK > 235.15) TK = 235.15;
            pGuess = HoltenEtAl2014_homogeneousIceNucleationPressureForTemperatureInK(TK);
            iter++;
        }
    }
    return TK;
}

static double getProperty(SWIM_Property valueOfProperty, SWIM_RegionType region, double t, double p) {
    static double HliqReferenceWagner2002  = -1.0;
    static double SliqReferenceWagner2002  = -1.0;
    static double CpliqReferenceWagner2002 = -1.0;
    double result = 0.0;
    if ( HliqReferenceWagner2002 == -1.0)  HliqReferenceWagner2002 = WagnerEtAl2002_getEnthalpy(298.15, 1.0);
    if ( SliqReferenceWagner2002 == -1.0)  SliqReferenceWagner2002 = WagnerEtAl2002_getEntropy(298.15, 1.0);
    if (CpliqReferenceWagner2002 == -1.0) CpliqReferenceWagner2002 = WagnerEtAl2002_getHeatCapacity(298.15, 1.0);
    switch (valueOfProperty) {
        case valueOfG:
        {
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getGibbsFreeEnergy(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getGibbsFreeEnergy(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getGibbsFreeEnergy(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getGibbsFreeEnergy(t, p);
            else                                 result = 0.0;
            if ((region == HoltenEtAl2014) || (region == WagnerEtAl2002)) {
                result +=     HliqReference - HliqReferenceWagner2002;
                result += -t*(SliqReference - SliqReferenceWagner2002);
            }
            result += GliqReference - (HliqReference - 298.15*SliqReference);
            break;
        }
        case valueOfH: {
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getEnthalpy(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getEnthalpy(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getEnthalpy(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getEnthalpy(t, p);
            else                                 result = 0.0;
            if ((region == HoltenEtAl2014) || (region == WagnerEtAl2002)) {
                result += HliqReference - HliqReferenceWagner2002;
            }
            break;
        }
        case valueOfS: {
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getEntropy(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getEntropy(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getEntropy(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getEntropy(t, p);
            else                                 result = 0.0;
            if ((region == HoltenEtAl2014) || (region == WagnerEtAl2002)) {
                result += SliqReference - SliqReferenceWagner2002;
            }
            break;
        }
        case valueOfCP:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getHeatCapacity(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getHeatCapacity(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getHeatCapacity(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getHeatCapacity(t, p);
            else                                 result = 0.0;
            break;
        case valueOfDCPDT:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getDcpDt(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getDcpDt(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getDcpDt(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getDcpDt(t, p);
            else                                 result = 0.0;
            break;
        case valueOfV:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getVolume(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getVolume(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getVolume(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getVolume(t, p);
            else                                 result = 0.0;
            break;
        case valueOfdVdT:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getDvDt(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getDvDt(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getDvDt(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getDvDt(t, p);
            else                                 result = 0.0;
            break;
        case valueOfdVdP:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getDvDp(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getDvDp(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getDvDp(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getDvDp(t, p);
            else                                 result = 0.0;
            break;
        case valueOfd2VdT2:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getD2vDt2(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getD2vDt2(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getD2vDt2(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getD2vDt2(t, p);
            else                                 result = 0.0;
            break;
        case valueOfd2VdTdP:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getD2vDtDp(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getD2vDtDp(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getD2vDtDp(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getD2vDtDp(t, p);
            else                                 result = 0.0;
            break;
        case valueOfd2VdP2:
            if      (region == DuanAndZhang2006) result = DuanAndZhang2006_getD2vDp2(t, p);
            else if (region == ZhangAndDuan2005) result = ZhangAndDuan2005_getD2vDp2(t, p);
            else if (region == HoltenEtAl2014)   result = HoltenEtAl2014_getD2vDp2(t, p);
            else if (region == WagnerEtAl2002)   result = WagnerEtAl2002_getD2vDp2(t, p);
            else                                 result = 0.0;
            break;
        default:
            result = 0.0;
            break;
    }
    return result;
}

static double smoothedProperty(double t, double p, SWIM_Property valueOfProperty, SWIM_RegionType region) {
    static double offsetsForZD2005AdjustmentRegion[11] = { -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 };
    static double pOfZD2005AdjustmentRegion = -1.0;
    double result = 0.0;
    if (region != NotApplicable) {
        switch (region) {
            case DuanAndZhang2006: {
                result = getProperty(valueOfProperty, DuanAndZhang2006, t, p);
                break;
            }
            case ZhangAndDuan2005: {
                result = getProperty(valueOfProperty, ZhangAndDuan2005, t, p);
                break;
            }
            case HoltenEtAl2014: {
                result = getProperty(valueOfProperty, HoltenEtAl2014, t, p);
                break;
            }
            case WagnerEtAl2002: {
                result = getProperty(valueOfProperty, WagnerEtAl2002, t, p);
                break;
            }
            default:
                break;
        }
        return  result;
    }
    region = classifyRegion(t, p);
    switch (region) {
        case DuanAndZhang2006: {
            double weightT = 1.0/2.0 - tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
            double weightP = 1.0/2.0 + tanh((p-centerOfDZ2006ToZD2005SmoothInterval)/widthOfDZ2006ToZD2005SmoothInterval)/2.0;
            if (weightP <= 0.001) weightP = 0.0;
            if (weightT <= 0.001) weightT = 0.0;
            if (weightP  > 0.001) result += weightP*getProperty(valueOfProperty, ZhangAndDuan2005, t, p);
            if (weightT  > 0.001) result += weightT*(1.0-weightP)*getProperty(valueOfProperty, WagnerEtAl2002, t, p);
            result += (1.0-weightP-weightT*(1.0-weightP))*getProperty(valueOfProperty, DuanAndZhang2006, t, p);
            break;
        }
        case ZhangAndDuan2005: {
            if (t < centerOfWagnerToDZ2006SmoothInterval) {
                double weightT = 1.0/2.0 + tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
                double weightP = 1.0/2.0 - tanh((p-centerOfWagnerToZD2005SmoothInterval)/widthOfWagnerToZD2005SmoothInterval)/2.0;
                if (weightP <= 0.001) weightP = 0.0;
                if (weightT <= 0.001) weightT = 0.0;
                if (weightP  > 0.001) result += weightP*(1.0-weightT)*getProperty(valueOfProperty, WagnerEtAl2002, t, p);
                if (weightT  > 0.001) result += weightP*weightT*getProperty(valueOfProperty, DuanAndZhang2006, t, p);
                if (t < highTofZD2005AdjustmentRegion) {
                    if (p != pOfZD2005AdjustmentRegion) {
                        double holtenAtLowT = getProperty(valueOfProperty, HoltenEtAl2014, lowTofZD2005AdjustmentRegion, p);
                        double ZD2005AtLowT = getProperty(valueOfProperty, ZhangAndDuan2005, lowTofZD2005AdjustmentRegion, p);
                        offsetsForZD2005AdjustmentRegion[valueOfProperty] = holtenAtLowT - ZD2005AtLowT;
                        pOfZD2005AdjustmentRegion = p;
                    }
                    double adjust = -offsetsForZD2005AdjustmentRegion[valueOfProperty]*(t-lowTofZD2005AdjustmentRegion)
                                     /(highTofZD2005AdjustmentRegion-lowTofZD2005AdjustmentRegion)
                                  + offsetsForZD2005AdjustmentRegion[valueOfProperty];
                    result += (1.0-weightP)*adjust;
                }
                result += (1.0-weightP)*getProperty(valueOfProperty, ZhangAndDuan2005, t, p);
            } else {
                double weightT = 1.0/2.0 - tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
                double weightP = 1.0/2.0 - tanh((p-centerOfDZ2006ToZD2005SmoothInterval)/widthOfDZ2006ToZD2005SmoothInterval)/2.0;
                if (weightP <= 0.001) weightP = 0.0;
                if (weightT <= 0.001) weightT = 0.0;
                if (weightP  > 0.001) result += weightP*(1.0-weightT)*getProperty(valueOfProperty, DuanAndZhang2006, t, p);
                if (weightT  > 0.001) result += weightP*weightT*getProperty(valueOfProperty, WagnerEtAl2002, t, p);
                result += (1.0-weightP)*getProperty(valueOfProperty, ZhangAndDuan2005, t, p);
            }
            break;
        }
        case HoltenEtAl2014: {
            result = getProperty(valueOfProperty, HoltenEtAl2014, t, p);
            break;
        }
        case WagnerEtAl2002: {
            double weightT = 1.0/2.0 + tanh((t-centerOfWagnerToDZ2006SmoothInterval)/widthOfWagnerToDZ2006SmoothInterval)/2.0;
            double weightP = 1.0/2.0 + tanh((p-centerOfWagnerToZD2005SmoothInterval)/widthOfWagnerToZD2005SmoothInterval)/2.0;
            if (weightP <= 0.001) weightP = 0.0;
            if (weightT <= 0.001) weightT = 0.0;
            if (weightP  > 0.001) {
                if (t < highTofZD2005AdjustmentRegion) {
                    if (p != pOfZD2005AdjustmentRegion) {
                        double holtenAtLowT = getProperty(valueOfProperty, HoltenEtAl2014, lowTofZD2005AdjustmentRegion, p);
                        double ZD2005AtLowT = getProperty(valueOfProperty, ZhangAndDuan2005, lowTofZD2005AdjustmentRegion, p);
                        offsetsForZD2005AdjustmentRegion[valueOfProperty] = holtenAtLowT - ZD2005AtLowT;
                        pOfZD2005AdjustmentRegion = p;
                    }
                    double adjust = -offsetsForZD2005AdjustmentRegion[valueOfProperty]*(t-lowTofZD2005AdjustmentRegion)
                                     /(highTofZD2005AdjustmentRegion-lowTofZD2005AdjustmentRegion)
                                  + offsetsForZD2005AdjustmentRegion[valueOfProperty];
                    result += weightP*adjust;
                }
                result += weightP*getProperty(valueOfProperty, ZhangAndDuan2005, t, p);
            }
            if (weightT  > 0.001) result += weightT*(1.0-weightP)*getProperty(valueOfProperty, DuanAndZhang2006, t, p);
            result += (1.0-weightP-weightT*(1.0-weightP))*getProperty(valueOfProperty, WagnerEtAl2002, t, p);
            break;
        }
        default:
            break;
    }
    return result;
}

double SWIM_getGibbsFreeEnergy(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfG, region);
}

double SWIM_getEnthalpy(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfH, region);
}

double SWIM_getEntropy(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfS, region);
}

double SWIM_getHeatCapacity(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfCP, region);
}

double SWIM_getDcpDt(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfDCPDT, region);
}

double SWIM_getVolume(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfV, region);
}

double SWIM_getDvDt(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfdVdT, region);
}

double SWIM_getDvDp(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfdVdP, region);
}

double SWIM_getD2vDt2(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfd2VdT2, region);
}

double SWIM_getD2vDtDp(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfd2VdTdP, region);
}

double SWIM_getD2vDp2(double t, double p, SWIM_RegionType region) {
    return smoothedProperty(t, p, valueOfd2VdP2, region);
}

/*
 * Convenience functions for density and density derivatives
*/

static const double MW = 18.01528;

double SWIM_getRho(double t, double p) {
    double v = 10.0*SWIM_getVolume(t, p, NotApplicable);
    return MW/v;
}

double SWIM_getDrhoDt(double t, double p) {
    double v    = 10.0*SWIM_getVolume(t, p, NotApplicable);
    double dvdt = 10.0*SWIM_getDvDt(t, p, NotApplicable);
    return -MW*dvdt/v/v;
}

double SWIM_getDrhoDp(double t, double p) {
    double v    = 10.0*SWIM_getVolume(t, p, NotApplicable);
    double dvdp = 10.0*SWIM_getDvDp(t, p, NotApplicable);
    return -MW*dvdp/v/v;
}

double SWIM_getD2rhoDt2(double t, double p) {
    double v      = 10.0*SWIM_getVolume(t, p, NotApplicable);
    double dvdt   = 10.0*SWIM_getDvDt(t, p, NotApplicable);
    double d2vdt2 = 10.0*SWIM_getD2vDt2(t, p, NotApplicable);
    return MW*(2.0*dvdt*dvdt/v/v/v - d2vdt2/v/v);
}

double SWIM_getD2rhoDtDp(double t, double p) {
    double v       = 10.0*SWIM_getVolume(t, p, NotApplicable);
    double dvdt    = 10.0*SWIM_getDvDt(t, p, NotApplicable);
    double dvdp    = 10.0*SWIM_getDvDp(t, p, NotApplicable);
    double d2vdtdp = 10.0*SWIM_getD2vDtDp(t, p, NotApplicable);
    return MW*(2.0*dvdt*dvdp/v/v/v - d2vdtdp/v/v);
}

double SWIM_getD2rhoDp2(double t, double p) {
    double v      = 10.0*SWIM_getVolume(t, p, NotApplicable);
    double dvdp   = 10.0*SWIM_getDvDp(t, p, NotApplicable);
    double d2vdp2 = 10.0*SWIM_getD2vDp2(t, p, NotApplicable);
    return MW*(2.0*dvdp*dvdp/v/v/v - d2vdp2/v/v);
}
