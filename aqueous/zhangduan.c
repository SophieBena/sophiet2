#include "zhangduan.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct _zAndDcorrections {
    double t;
    double g;
    double h;
    double s;
    double cp;
    double dcpdt;
} ZandDcorrections;

const static ZandDcorrections zAndDcorrections[] = {
    { 270, 45039.85053, 163354.1531, 438.2011291, -1486.30666, 17.92271899 },
    { 280, 40919.27332, 149343.9978, 387.2311428, -1319.863346, 15.43848334 },
    { 290, 37271.18718, 136880.9981, 343.482104, -1176.211447, 13.35210746 },
    { 300, 34029.611, 125756.01, 305.75467, -1051.7193, 11.592283 },
    { 310, 31139.313, 115792.62, 273.0752, -943.43561, 10.10518 },
    { 320, 28553.945, 106841.62, 244.64898, -848.87035, 8.8417461 },
    { 330, 26234.341, 98776.334, 219.82422, -765.98133, 7.7642901 },
    { 340, 24147.258, 91488.771, 198.06327, -693.06794, 6.8419677 },
    { 350, 22264.357, 84886.483, 178.92036, -628.70984, 6.049418 },
    { 360, 20561.369, 78890.045, 162.0241, -571.71678, 5.3658192 },
    { 370, 19017.428, 73430.952, 147.06358, -521.08773, 4.7740205 },
    { 380, 17614.521, 68449.908, 133.77733, -475.97773, 4.2598701 },
    { 390, 16337.04, 63895.398, 121.94451, -435.67093, 3.8116101 },
    { 400, 15171.413, 59722.517, 111.37776, -399.55869, 3.4194793 },
    { 410, 14105.799, 55891.982, 101.91752, -367.12165, 3.0753259 },
    { 420, 13129.832, 52369.323, 93.42736, -337.91513, 2.772325 },
    { 430, 12234.41, 49124.191, 85.790187, -311.5571, 2.5047538 },
    { 440, 11411.522, 46129.788, 78.90515, -287.71828, 2.2677528 },
    { 450, 10654.092, 43362.381, 72.685087, -266.1139, 2.0572376 },
    { 460, 9955.8569, 40800.889, 67.054417, -246.49696, 1.8697293 },
    { 470, 9311.2591, 38426.536, 61.947398, -228.65252, 1.7022843 },
    { 480, 8715.3554, 36222.558, 57.306672, -212.39295, 1.5523746 },
    { 490, 8163.7394, 34173.944, 53.082051, -197.55398, 1.4178481 },
    { 500, 7652.4748, 32267.226, 49.229502, -183.99128, 1.2968587 },
    { 510, 7178.0387, 30490.29, 45.710296, -171.57766, 1.1876457 },
    { 520, 6737.272, 28832.219, 42.490283, -160.20059, 1.0891408 },
    { 530, 6327.3367, 27283.157, 39.539285, -149.76013, 1.0000356 },
    { 540, 5945.679, 25834.193, 36.830582, -140.16711, 0.9193478 },
    { 550, 5589.9967, 24477.258, 34.340475, -131.34153, 0.846215 },
    { 560, 5258.2112, 23205.046, 32.047919, -123.2111, 0.7799196 },
    { 570, 4948.4424, 22010.939, 29.934204, -115.70995, 0.7198632 },
    { 580, 4658.9866, 20888.952, 27.982699, -108.77733, 0.6655649 },
    { 590, 4388.2972, 19833.687, 26.178627, -102.35631, 0.6166528 },
    { 600, 4134.9665, 18840.302, 24.508894, -96.392376, 0.5728691 },
    { 610, 3897.7102, 17904.495, 22.961942, -90.831803, 0.5340626 },
    { 620, 3675.3522, 17022.505, 21.527665, -85.619549, 0.5001807 },
    { 630, 3466.7324, 16190.873, 20.197049, -80.759613, 0.4555729 },
    { 640, 3271.0205, 15405.433, 18.960019, -76.398056, 0.413445 },
    { 650, 3087.2512, 14662.05, 17.807383, -72.308197, 0.391237 },
    { 660, 2914.6103, 13958.72, 16.733499, -68.386476, 0.369032 },
    { 670, 2752.3353, 13293.689, 15.733364, -64.654187, 0.3447 },
    { 680, 2599.7147, 12664.9, 14.801743, -61.141838, 0.320853 },
    { 690, 2425.3585, 12126.76, 14.060002, -55.205467, 0.299928 },
    { 700, 2288.6708, 11589.025, 13.28622, -52.363481, 0.282965 },
    { 710, 2159.4661, 11079.091, 12.562853, -49.641338, 0.269871 },
    { 720, 2037.2555, 10595.872, 11.886968, -47.017116, 0.259967 },
    { 730, 1921.5769, 10138.486, 11.25604, -44.471907, 0.252337 },
    { 740, 1811.9922, 9706.2153, 10.667869, -41.992614, 0.245864 },
    { 750, 1708.0837, 9298.4262, 10.120457, -39.576127, 0.239112 },
    { 760, 1609.4537, 8914.4454, 9.611831, -37.234529, 0.230222 },
    { 770, 1515.725, 8553.3814, 9.139814, -34.999948, 0.2170172 },
    { 780, 1426.5442, 8213.9122, 8.701754, -32.926219, 0.197453 },
    { 790, 1341.5878, 7894.0871, 8.294302, -31.083562, 0.1704893 },
    { 800, 1260.5698, 7591.2297, 7.913325, -29.543588, 0.1370493 },
    { 810, 1183.2489, 7302.0335, 7.554055, -28.356667, 0.1004524 },
    { 820, 1109.4329, 7022.8885, 7.211531, -27.53029, 0.06572863 },
    { 830, 1038.9778, 6750.3691, 6.881195, -27.020039, 0.0378869 },
    { 840, 971.78056, 6481.7207, 6.559452, -26.739252, 0.0201114 },
    { 850, 907.76791, 6215.1715, 6.244004, -26.582479, 0.012818 },
    { 860, 846.88296, 5949.9905, 5.933846, -26.450714, 0.014598 },
    { 870, 789.07326, 5686.3292, 5.62903, -26.268364, 0.022487 },
    { 880, 734.2821, 5424.9497, 5.330304, -25.988983, 0.033602 },
    { 890, 682.44324, 5166.942, 5.038762, -25.592517, 0.045609 },
    { 900, 633.479, 4913.4924, 4.755571, -25.078549, 0.056923 },
    { 910, 587.30046, 4665.7236, 4.481784, -24.459007, 0.066639 },
    { 920, 543.80911, 4424.6035, 4.218255, -23.752141, 0.074368 },
    { 930, 502.8991, 4190.9041, 3.965597, -22.978269, 0.0800603 },
    { 940, 464.45965, 3965.1955, 3.724187, -22.157128, 0.0838613 },
    { 950, 428.37729, 3747.8597, 3.494192, -21.306448, 0.0860158 },
    { 960, 394.53779, 3539.1143, 3.275601, -20.441341, 0.0867958 },
    { 970, 362.82777, 3339.0394, 3.068259, -19.574185, 0.0864701 },
    { 980, 333.13596, 3147.6044, 2.871907, -18.714796, 0.0852809 },
    { 990, 305.35414, 2964.6922, 2.6862, -17.870725, 0.0834389 },
    { 1000, 279.37779, 2790.1199, 2.510742, -17.047597, 0.0811196 },
    { 1010, 255.10658, 2623.6569, 2.345099, -16.249444, 0.0784656 },
    { 1020, 232.4446, 2465.0385, 2.188818, -15.47902, 0.0755912 },
    { 1030, 211.30057, 2313.9782, 2.041435, -14.738058, 0.0725864 },
    { 1040, 191.58779, 2170.179, 1.902488, -14.027497, 0.0695209 },
    { 1050, 173.22419, 2033.3258, 1.771525, -13.347662, 0.0664483 },
    { 1060, 156.13215, 1903.1207, 1.648102, -12.698419, 0.0634084 },
    { 1070, 140.23843, 1779.257, 1.531793, -12.079288, 0.0604301 },
    { 1080, 125.47394, 1661.437, 1.422188, -11.48954, 0.0575346 },
    { 1090, 111.77362, 1549.3712, 1.318896, -10.928271, 0.0547363 },
    { 1100, 99.076198, 1442.78, 1.221549, -10.394457, 0.0520446 },
    { 1110, 87.324008, 1341.3942, 1.129793, -9.887003, 0.0494652 },
    { 1120, 76.462805, 1244.9559, 1.043298, -9.404771, 0.0470005 },
    { 1130, 66.441569, 1153.2186, 0.961749, -8.94661, 0.0446509 },
    { 1140, 57.212304, 1065.9473, 0.884856, -8.511372, 0.0424154 },
    { 1150, 48.729875, 982.91848, 0.812337, -8.09793, 0.04029147 },
    { 1160, 40.951824, 903.91971, 0.743938, -7.705183, 0.03827585 },
    { 1170, 33.838209, 828.7494, 0.679411, -7.332065, 0.03636472 },
    { 1180, 27.351451, 757.2164, 0.61853, -6.977553, 0.03455395 },
    { 1190, 21.456183, 689.1396, 0.561079, -6.640665, 0.03283914 },
    { 1200, 16.119117, 624.34747, 0.506857, -6.320466, 0.03121563 },
    { 1210, 11.308907, 562.67763, 0.455676, -6.016062, 0.02967904 },
    { 1220, 6.996036, 503.97639, 0.407361, -5.726609, 0.02822489 },
    { 1230, 3.152692, 448.09829, 0.361744, -5.451304, 0.02684882 },
    { 1240, -0.24733, 394.90568, 0.318672, -5.189387, 0.02554652 },
    { 1250, -3.228737, 344.26831, 0.277997, -4.940141, 0.02431395 },
    { 1260, -5.814822, 296.06289, 0.239586, -4.702888, 0.02314733 },
    { 1270, -8.027559, 250.17271, 0.203308, -4.476988, 0.02204285 },
    { 1280, -9.88768, 206.48731, 0.169043, -4.261836, 0.02099694 },
    { 1290, -11.414754, 164.9027, 0.13668, -4.056864, 0.02000638 },
    { 1300, -12.627258, 125.3179, 0.106112, -3.861534, 0.019067832 },
    { 1310, -13.54264, 87.640922, 0.07724, -3.675342, 0.018178488 },
    { 1320, -14.177385, 51.782183, 0.049969, -3.497811, 0.017335461 },
    { 1330, -14.547075, 17.657352, 0.024214, -3.328488, 0.016536028 },
    { 1340, -14.666435, -14.813535, -0.000109, -3.166953, 0.015777733 },
    { 1350, -14.549392, -45.706326, -0.023079, -3.012804, 0.015058208 },
    { 1360, -14.209121, -75.092985, -0.044768, -2.865666, 0.01437538 },
    { 1370, -13.658082, -103.04182, -0.065244, -2.725182, 0.01372699 },
    { 1380, -12.908073, -129.61769, -0.084572, -2.591018, 0.013111092 },
    { 1390, -11.970259, -154.88219, -0.102814, -2.462857, 0.01252598 },
    { 1400, -10.855215, -178.89385, -0.120027, -2.340401, 0.01196996 },
    { 1410, -9.572955, -201.70829, -0.136266, -2.223368, 0.01144119 },
    { 1420, -8.132969, -223.37839, -0.151581, -2.111491, 0.01093835 },
    { 1430, -6.544248, -243.95445, -0.166022, -2.004519, 0.01045995 },
    { 1440, -4.815313, -263.48433, -0.179632, -1.902215, 0.01000466 },
    { 1450, -2.954245, -282.01355, -0.192455, -1.804353, 0.00957121 },
    { 1460, -0.968705, -299.58549, -0.204532, -1.710721, 0.00915842 },
    { 1470, 1.134037, -316.24142, -0.215902, -1.621119, 0.0087652 },
    { 1480, 3.34709, -332.0267, -0.226599, -1.535356, 0.00839048 },
    { 1490, 5.663911, -346.9673, -0.236661, -1.453252, 0.00803328 },
    { 1500, 8.078291, -361.09733, -0.246117, -1.374635, 0.00769269 },
    { 1510, 10.584339, -374.46452, -0.254999, -1.299346, 0.00736782 },
    { 1520, 13.17646, -387.09481, -0.263337, -1.227229, 0.00705786 },
    { 1530, 15.849342, -399.0192, -0.271156, -1.158141, 0.00676204 },
    { 1540, 18.597939, -410.26727, -0.278484, -1.091944, 0.00647963 },
    { 1550, 21.417463, -420.86727, -0.285345, -1.028506, 0.00620994 },
    { 1560, 24.303361, -430.84617, -0.291762, -0.967704, 0.00595233 },
    { 1570, 27.251309, -440.22975, -0.297759, -0.909421, 0.0057062 },
    { 1580, 30.257198, -449.04261, -0.303355, -0.853544, 0.00547096 },
    { 1590, 33.317127, -457.30829, -0.30857, -0.799967, 0.00524608 },
    { 1600, 36.427382, -465.04929, -0.313423, -0.74859, 0.00503103 },
    { 1610, 39.584439, -472.2871, -0.317932, -0.699316, 0.00482534 },
    { 1620, 42.784946, -479.0423, -0.322116, -0.652053, 0.00462855 },
    { 1630, 46.025716, -485.33458, -0.325988, -0.606716, 0.00444023 },
    { 1640, 49.30372, -491.18277, -0.329565, -0.563222, 0.004259965 },
    { 1650, 52.616079, -496.60489, -0.332862, -0.521492, 0.004087377 },
    { 1660, 55.960054, -501.61822, -0.335891, -0.48145, 0.003922098 },
    { 1670, 59.333043, -506.23928, -0.338666, -0.443026, 0.003763779 },
    { 1680, 62.732571, -510.48391, -0.3412, -0.406152, 0.0036121 },
    { 1690, 66.156284, -514.36727, -0.343505, -0.370763, 0.003466744 },
    { 1700, 69.601945, -517.90391, -0.345592, -0.336797, 0.003327412 },
    { 1710, 73.067426, -521.10776, -0.34747, -0.304195, 0.003193832 },
    { 1720, 76.550705, -523.99218, -0.349153, -0.272902, 0.003065745 },
    { 1730, 80.049858, -526.56998, -0.350648, -0.242863, 0.002942895 },
    { 1740, 83.563057, -528.85346, -0.351964, -0.214027, 0.002825042 },
    { 1750, 87.088564, -530.85438, -0.353111, -0.186346, 0.00271196 },
    { 1760, 90.624725, -532.58408, -0.354096, -0.159773, 0.00260345 },
    { 1770, 94.169972, -534.05339, -0.354928, -0.134263, 0.00249929 },
    { 1780, 97.722809, -535.27273, -0.355615, -0.109773, 0.00239929 },
    { 1790, 101.28182, -536.25212, -0.356164, -0.086264, 0.00230327 },
    { 1800, 104.84565, -537.00114, -0.356582, -0.063695, 0.00221107 },
    { 1810, 108.41303, -537.52903, -0.356874, -0.04203, 0.00212249 },
    { 1820, 111.98274, -537.84464, -0.357048, -0.021234, 0.0020374 },
    { 1830, 115.55362, -537.95649, -0.357109, -0.001271, 0.00195563 },
    { 1840, 119.12457, -537.87274, -0.357064, 0.01789, 0.00187706 },
    { 1850, 122.69456, -537.60126, -0.356917, 0.03628, 0.00180154 },
    { 1860, 126.26258, -537.14961, -0.356673, 0.05393, 0.00172894 },
    { 1870, 129.82772, -536.52503, -0.356338, 0.070868, 0.00165914 },
    { 1880, 133.38906, -535.73452, -0.355917, 0.087121, 0.00159202 },
    { 1890, 136.94578, -534.78479, -0.355413, 0.102717, 0.00152746 },
    { 1900, 140.49706, -533.68229, -0.354831, 0.117679, 0.00146538 },
    { 1910, 144.04216, -532.43324, -0.354176, 0.132032, 0.00140564 },
    { 1920, 147.58034, -531.04359, -0.35345, 0.1458, 0.00134818 },
    { 1930, 151.11094, -529.51912, -0.352658, 0.159003, 0.00129288 },
    { 1940, 154.63329, -527.86533, -0.351803, 0.171664, 0.00123967 },
    { 1950, 158.14681, -526.08757, -0.35089, 0.183804, 0.00118846 },
    { 1960, 161.65089, -524.1994, -0.349919, 0.19544, 0.00113915 },
    { 1970, 165.14502, -522.1838, -0.348896, 0.206592, 0.00109169 },
    { 1980, 168.62865, -520.06, -0.347822, 0.217279, 0.00104599 },
    { 1990, 172.10132, -517.83628, -0.346702, 0.227518, 0.00100199 },
    { 2000, 175.56255, -515.51172, -0.345537, 0.237325, 0.00095961 }
};

const static int nZandDcorrections = (sizeof(zAndDcorrections)/sizeof(struct _zAndDcorrections));

static const double R = 8.3143;

/* pure EOS terms for 0.1 to 10 GPa */
static const double H2Oa1  =  3.49824207e-01;
static const double H2Oa2  = -2.91046273e+00;
static const double H2Oa3  =  2.00914688e+00;
static const double H2Oa4  =  1.12819964e-01;
static const double H2Oa5  =  7.48997714e-01;
static const double H2Oa6  = -8.73207040e-01;
static const double H2Oa7  =  1.70609505e-02;
static const double H2Oa8  = -1.46355822e-02;
static const double H2Oa9  =  5.79768283e-02;
static const double H2Oa10 = -8.41246372e-04;
static const double H2Oa11 =  4.95186474e-03;
static const double H2Oa12 = -9.16248538e-03;
static const double H2Oa13 = -1.00358152e-01;
static const double H2Oa14 = -1.82674744e-03;
static const double H2Ogam =  1.05999998e-02;

/* H2O, critical constants, K, bars, J/bar */

static const double H2OTc = 647.25;
static const double H2OVc =  55.9480373/10.0;
#define H2OVPc (8.314467*H2OTc/H2OVc)

static const double idealCoeff[13] = {
      3.10409601236035e+01,
     -3.91422080460869e+01,
      3.79695277233575e+01,
     -2.18374910952284e+01,
      7.42251494566339e+00,
     -1.38178929609470e+00,
      1.08807067571454e-01,
     -1.20771176848589e+01,
      3.39105078851732e+00,
     -5.84520979955060e-01,
      5.89930846488082e-02,
     -3.12970001415882e-03,
      6.57460740981757e-05
};

static void BVcAndDerivative(double t, double *bv, double *dbvdt, double *d2bvdt2, double *d3bvdt3) {
    double H2OTr    = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double bEnd = H2Oa1 + H2Oa2/H2OTr/H2OTr + H2Oa3/H2OTr/H2OTr/H2OTr;
    double dbEnddt = - 2.0*H2Oa2*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa3*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2bEnddt2 = 6.0*H2Oa2*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa3*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3bEnddt3 = - 24.0*H2Oa2*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr
                       - 60.0*H2Oa3*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *bv = bEnd*H2OVc;
    *dbvdt   = dbEnddt*H2OVc;
    *d2bvdt2 = d2bEnddt2*H2OVc;
    *d3bvdt3 = d3bEnddt3*H2OVc;

    return;
}

static void CVcAndDerivative(double t, double *cv, double *dcvdt, double *d2cvdt2, double *d3cvdt3) {
    double H2OTr    = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double cEnd = H2Oa4 + H2Oa5/H2OTr/H2OTr + H2Oa6/H2OTr/H2OTr/H2OTr;
    double dcEnddt = - 2.0*H2Oa5*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa6*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2cEnddt2 = 6.0*H2Oa5*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa6*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3cEnddt3 = -24.0*H2Oa5*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2Oa6*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *cv = cEnd*H2OVc*H2OVc;
    *dcvdt   = dcEnddt*H2OVc*H2OVc;
    *d2cvdt2 = d2cEnddt2*H2OVc*H2OVc;
    *d3cvdt3 = d3cEnddt3*H2OVc*H2OVc;

    return;
}

static void DVcAndDerivative(double t, double *dv, double *ddvdt, double *d2dvdt2, double *d3dvdt3) {
    double H2OTr = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double dEnd = H2Oa7 + H2Oa8/H2OTr/H2OTr + H2Oa9/H2OTr/H2OTr/H2OTr;
    double ddEnddt = - 2.0*H2Oa8*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa9*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2dEnddt2 = 6.0*H2Oa8*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa9*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3dEnddt3 = - 24.0*H2Oa8*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr - 60.0*H2Oa9*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *dv = dEnd*H2OVc*H2OVc*H2OVc*H2OVc;
    *ddvdt   = ddEnddt*H2OVc*H2OVc*H2OVc*H2OVc;
    *d2dvdt2 = d2dEnddt2*H2OVc*H2OVc*H2OVc*H2OVc;
    *d3dvdt3 = d3dEnddt3*H2OVc*H2OVc*H2OVc*H2OVc;

    return;
}

static void EVcAndDerivative(double t, double *ev, double *devdt, double *d2evdt2, double *d3evdt3) {
    double H2OTr = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double eEnd = H2Oa10 + H2Oa11/H2OTr/H2OTr + H2Oa12/H2OTr/H2OTr/H2OTr;
    double deEnddt = - 2.0*H2Oa11*dH2OTrdt/H2OTr/H2OTr/H2OTr - 3.0*H2Oa12*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;
    double d2eEnddt2 = 6.0*H2Oa11*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr + 12.0*H2Oa12*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;
    double d3eEnddt3 = - 24.0*H2Oa11*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr
                     - 60.0*H2Oa12*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr/H2OTr;

    *ev = eEnd*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
    *devdt = deEnddt*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
    *d2evdt2 = d2eEnddt2*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;
    *d3evdt3 = d3eEnddt3*H2OVc*H2OVc*H2OVc*H2OVc*H2OVc;

    return;
}

static void FVcAndDerivative(double t, double *fv, double *dfvdt, double *d2fvdt2, double *d3fvdt3) {
    double H2OTr = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double fEnd      = H2Oa13/H2OTr;
    double dfEnddt   = -     H2Oa13*dH2OTrdt/H2OTr/H2OTr;
    double d2fEnddt2 =   2.0*H2Oa13*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr;
    double d3fEnddt3 = - 6.0*H2Oa13*dH2OTrdt*dH2OTrdt*dH2OTrdt/H2OTr/H2OTr/H2OTr/H2OTr;

    *fv = fEnd*H2OVc*H2OVc;
    *dfvdt   = dfEnddt*H2OVc*H2OVc;
    *d2fvdt2 = d2fEnddt2*H2OVc*H2OVc;
    *d3fvdt3 = d3fEnddt3*H2OVc*H2OVc;

    return;
}

static void GVcAndDerivative(double t, double *gv, double *dgvdt, double *d2gvdt2, double *d3gvdt3) {
    double H2OTr = t/H2OTc;
    double dH2OTrdt = 1.0/H2OTc;

    double gEnd      = H2Oa14*H2OTr;
    double dgEnddt   = H2Oa14*dH2OTrdt;
    double d2gEnddt2 = 0.0;
    double d3gEnddt3 = 0.0;

    *gv = gEnd*H2OVc*H2OVc*H2OVc*H2OVc;
    *dgvdt   = dgEnddt*H2OVc*H2OVc*H2OVc*H2OVc;
    *d2gvdt2 = d2gEnddt2*H2OVc*H2OVc*H2OVc*H2OVc;
    *d3gvdt3 = d3gEnddt3*H2OVc*H2OVc*H2OVc*H2OVc;

    return;
}

static void GammaVcAndDerivative(double t, double *gammav) {
    *gammav = H2Ogam*H2OVc*H2OVc;

    return;
}

static void idealGasH2O(double t, double *cp, double *s0, double *h0, double *dcpdt) {
    int i;

    for (i=0, *cp=0.0; i<7; i++) *cp += idealCoeff[i]*pow(t/1000.0, (double) i);
    for (i=7; i<13; i++)         *cp += idealCoeff[i]/pow(t/1000.0, (double) (i-6));

    for (i=1, *dcpdt=0.0; i<7; i++) *dcpdt +=  ((double) i)  *idealCoeff[i]*pow(t/1000.0, (double) i-1);
    for (i=7; i<13; i++)            *dcpdt += -((double) i-6)*idealCoeff[i]/pow(t/1000.0, (double) (i+1-6));

    for (i=0, *h0=0.0; i<7; i++) *h0 += idealCoeff[i]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += idealCoeff[7]*log(t/1000.0);
    for (i=8; i<13; i++)         *h0 += idealCoeff[i]/pow(t/1000.0, (double) (i-7))/((double) (7-i));

    *s0  = idealCoeff[0]*log(t/1000.0);
    for (i=1; i<7; i++)	       *s0 += idealCoeff[i]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++)         *s0 += idealCoeff[i]/pow(t/1000.0, (double) (i-6))/((double) (6-i));

    *cp    *= 8.31451;
    *h0    *= 8.31451*1000.0;
    *s0    *= 8.31451;
    *dcpdt *= 8.31451/1000.0;

    *h0 += - 355665.4136;
    *s0 +=   359.6505;
}

static void dewH2ODriver(double t, double p,
                          double *vPt, double *zPt, double *phi, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2,
                          double *dlnphidt, double *dlnphidp, double *d2lnphidt2, double *d2lnphidtdp, double *d2lnphidp2) {
    double bv, cv, dv, ev, fv, gv, gammav, v, z = 1.0, dzdv, dzdt, d2zdv2, d2zdvdt, d2zdt2;
    double dbvdt, dcvdt, ddvdt, devdt, dfvdt, dgvdt;
    double d2bvdt2, d2cvdt2, d2dvdt2, d2evdt2, d2fvdt2, d2gvdt2;
    double d3bvdt3, d3cvdt3, d3dvdt3, d3evdt3, d3fvdt3, d3gvdt3;
    double lnPhiH2O, dlnPhiH2Odv, dlnPhiH2Odt, d2lnPhiH2Odv2, d2lnPhiH2Odt2, d2lnPhiH2Odvdt;

    BVcAndDerivative    (t, &bv,     &dbvdt,    &d2bvdt2,  &d3bvdt3);
    CVcAndDerivative    (t, &cv,     &dcvdt,    &d2cvdt2,  &d3cvdt3);
    DVcAndDerivative    (t, &dv,     &ddvdt,    &d2dvdt2,  &d3dvdt3);
    EVcAndDerivative    (t, &ev,     &devdt,    &d2evdt2,  &d3evdt3);
    FVcAndDerivative    (t, &fv,     &dfvdt,    &d2fvdt2,  &d3fvdt3);
    GVcAndDerivative    (t, &gv,     &dgvdt,    &d2gvdt2,  &d3gvdt3);
    GammaVcAndDerivative(t, &gammav);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0, vLowest = 1.7;

        v = 8.314467*t/p;
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v + gv/v/v/v/v) * exp(-gammav/v/v);
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            if (v < 1.0) v = vLowest;
            if ((t < 473.15) && (v > 2.0)) v = 2.0;  /* MSG 8-17-18 */
            iter++;
            if (iter == 200) {
                printf("dewH2ODriver(a): t = %g, p = %g, z = %g, v = %g, vPrev = %g, delv = %g, delvPrev = %g, iter = %d\n", t, p, z, v, vPrevious, delv, delvPrevious, iter);
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
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v + gv/v/v/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
            if ( (iter == 200) || (fabs(dx) > 100.0*DBL_EPSILON) ) printf("dewH2ODriver(b): t = %g, p = %g, z = %g, v = %g, delv = %g, iter = %d\n", t, p, z, v, delv, iter);
        }
    }

    lnPhiH2O  = 0.0;
    lnPhiH2O += -log(z) + z - 1.0;
    lnPhiH2O += bv/v;
    lnPhiH2O += cv/2.0/v/v;
    lnPhiH2O += dv/4.0/v/v/v/v;
    lnPhiH2O += ev/5.0/v/v/v/v/v;
    lnPhiH2O += (fv/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiH2O += (gv/2.0/gammav/gammav)*(1.0 - (gammav/v/v+1.0)*exp(-gammav/v/v));

    // z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v + gv/v/v/v/v) * exp(-gammav/v/v);
    dzdv = - bv/v/v - 2.0*cv/v/v/v + - 4.0*dv/v/v/v/v/v - 5.0*ev/v/v/v/v/v/v
         - 2.0*(fv/v/v/v) * exp(-gammav/v/v) + 2.0*(fv/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
         - 4.0*(gv/v/v/v/v/v) * exp(-gammav/v/v) + 2.0*(gv/v/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    dzdt = dbvdt/v + dcvdt/v/v + ddvdt/v/v/v/v + devdt/v/v/v/v/v + (dfvdt/v/v + dgvdt/v/v/v/v) * exp(-gammav/v/v);

    dlnPhiH2Odv  = 0.0;
    dlnPhiH2Odv += -dzdv/z + dzdv;
    dlnPhiH2Odv += -bv/v/v;
    dlnPhiH2Odv += -cv/v/v/v;
    dlnPhiH2Odv += -dv/v/v/v/v/v;
    dlnPhiH2Odv += -ev/v/v/v/v/v/v;
    dlnPhiH2Odv += -(fv/v/v/v)*exp(-gammav/v/v);
    dlnPhiH2Odv += -(gv/v/v/v/v/v)*exp(-gammav/v/v);

    dlnPhiH2Odt  = 0.0;
    dlnPhiH2Odt += -dzdt/z + dzdt;
    dlnPhiH2Odt += dbvdt/v;
    dlnPhiH2Odt += dcvdt/2.0/v/v;
    dlnPhiH2Odt += ddvdt/4.0/v/v/v/v;
    dlnPhiH2Odt += devdt/5.0/v/v/v/v/v;
    dlnPhiH2Odt += (dfvdt/2.0/gammav)*(1.0-exp(-gammav/v/v));
    dlnPhiH2Odt += (dgvdt/2.0/gammav/gammav)*(1.0 - (gammav/v/v+1.0)*exp(-gammav/v/v));

    d2zdv2 =  2.0*bv/v/v/v + 6.0*cv/v/v/v/v + 20.0*dv/v/v/v/v/v/v + 30.0*ev/v/v/v/v/v/v/v
           + 6.0*(fv/v/v/v/v) * exp(-gammav/v/v) - 4.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
           - 4.0*(fv/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
           - 6.0*(fv/v/v) * (gammav/v/v/v/v) * exp(-gammav/v/v)
           + 4.0*(fv/v/v) * (gammav*gammav/v/v/v/v/v/v) * exp(-gammav/v/v)
           + 20.0*(gv/v/v/v/v/v/v) * exp(-gammav/v/v) - 8.0*(gv/v/v/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
           - 8.0*(gv/v/v/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
           - 6.0*(gv/v/v/v/v) * (gammav/v/v/v/v) * exp(-gammav/v/v)
           + 4.0*(gv/v/v/v/v) * (gammav/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);

    d2zdvdt = - dbvdt/v/v - 2.0*dcvdt/v/v/v + - 4.0*ddvdt/v/v/v/v/v - 5.0*devdt/v/v/v/v/v/v
            - 2.0*(dfvdt/v/v/v) * exp(-gammav/v/v) + 2.0*(dfvdt/v/v) * (gammav/v/v/v) * exp(-gammav/v/v)
            - 4.0*(dgvdt/v/v/v/v/v) * exp(-gammav/v/v) + 2.0*(dgvdt/v/v/v/v) * (gammav/v/v/v) * exp(-gammav/v/v);
    d2zdt2 = d2bvdt2/v + d2cvdt2/v/v + d2dvdt2/v/v/v/v + d2evdt2/v/v/v/v/v + (d2fvdt2/v/v + d2gvdt2/v/v/v/v) * exp(-gammav/v/v);

    d2lnPhiH2Odv2  = 0.0;
    d2lnPhiH2Odv2 += dzdv*dzdv/z/z - d2zdv2/z + d2zdv2;
    d2lnPhiH2Odv2 +=  2.0*bv/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*cv/v/v/v/v;
    d2lnPhiH2Odv2 +=  5.0*dv/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  6.0*ev/v/v/v/v/v/v/v;
    d2lnPhiH2Odv2 +=  3.0*(fv/v/v/v/v)*exp(-gammav/v/v) - 2.0*(fv/v/v/v)*(gammav/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odv2 +=  5.0*(gv/v/v/v/v/v/v)*exp(-gammav/v/v) - 2.0*(gv/v/v/v/v/v)*(gammav/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odvdt  = 0.0;
    d2lnPhiH2Odvdt += dzdv*dzdt/z/z -d2zdvdt/z + d2zdvdt;
    d2lnPhiH2Odvdt += -dbvdt/v/v;
    d2lnPhiH2Odvdt += -dcvdt/v/v/v;
    d2lnPhiH2Odvdt += -ddvdt/v/v/v/v/v;
    d2lnPhiH2Odvdt += -devdt/v/v/v/v/v/v;
    d2lnPhiH2Odvdt += -(dfvdt/v/v/v)*exp(-gammav/v/v);
    d2lnPhiH2Odvdt += -(dgvdt/v/v/v/v/v)*exp(-gammav/v/v);

    d2lnPhiH2Odt2  = 0.0;
    d2lnPhiH2Odt2 += dzdt*dzdt/z/z - d2zdt2/z + d2zdt2;
    d2lnPhiH2Odt2 += d2bvdt2/v;
    d2lnPhiH2Odt2 += d2cvdt2/2.0/v/v;
    d2lnPhiH2Odt2 += d2dvdt2/4.0/v/v/v/v;
    d2lnPhiH2Odt2 += d2evdt2/5.0/v/v/v/v/v;
    d2lnPhiH2Odt2 += (d2fvdt2/2.0/gammav)*(1.0-exp(-gammav/v/v));
    d2lnPhiH2Odt2 += (d2gvdt2/2.0/gammav/gammav)*(1.0 - (gammav/v/v+1.0)*exp(-gammav/v/v));

    *vPt      = v;
    *zPt      = z;
    *phi      = exp(lnPhiH2O);
    *dvdp     = 1.0/( p*(dzdv/z - 1.0/v) );
    *dvdt     = (1.0/t + dzdt/z)/(1.0/v - dzdv/z);
    *d2vdp2   = p*(1.0/v-dzdv/z)/v - dzdv/(*dvdp)/z + 1.0/(*dvdp)/(*dvdp)/p + p*d2zdv2/z;
    *d2vdp2  *= -(*dvdp)*(*dvdp)*(*dvdp);
    *d2vdtdp  = -(*dvdp)*(1.0/t + dzdt/z) - p*(*d2vdp2)*(1.0/t + dzdt/z) - p*(*dvdp)*(*dvdp)*(-dzdv*dzdt/z/z + d2zdvdt/z);
    *d2vdt2   = -p*(*d2vdtdp)*(1.0/t + dzdt/z) + p*(*dvdp)/t/t + p*(*dvdp)*dzdt*(dzdv*(*dvdt) + dzdt)/z/z
              - p*(*dvdp)*(*dvdt)*d2zdvdt/z - p*(*dvdp)*d2zdt2/z;

    *dlnphidt    = dlnPhiH2Odv*(*dvdt) + dlnPhiH2Odt;
    *dlnphidp    = dlnPhiH2Odv*(*dvdp);
    *d2lnphidt2  = d2lnPhiH2Odv2*(*dvdt)*(*dvdt) + 2.0*d2lnPhiH2Odvdt*(*dvdt) + dlnPhiH2Odv*(*d2vdt2) + d2lnPhiH2Odt2;
    *d2lnphidtdp = d2lnPhiH2Odv2*(*dvdt)*(*dvdp) + dlnPhiH2Odv*(*d2vdtdp) + d2lnPhiH2Odvdt*(*dvdp);
    *d2lnphidp2  = d2lnPhiH2Odv2*(*dvdp)*(*dvdp) + dlnPhiH2Odv*(*d2vdp2);
}

static void propertiesOfPureH2O(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt,
                                double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    double z, phi, dlnphidt, dlnphidp, d2lnphidt2, d2lnphidtdp, d2lnphidp2;

    idealGasH2O(t, cp, s, h, dcpdt);
    dewH2ODriver(t, p, v, &z, &phi, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, &dlnphidt, &dlnphidp, &d2lnphidt2, &d2lnphidtdp, &d2lnphidp2);

    *g      = *h - t*(*s) + R*t*log(phi*p);
    *s     += - (R*log(phi*p) + R*t*dlnphidt);
    *h     += R*t*log(phi*p) - t*(R*log(phi*p) + R*t*dlnphidt);
    *cp    += -t*(2.0*R*dlnphidt + R*t*d2lnphidt2);
    {
        double zTemp, phiTemp, dlnphidtTemp, dlnphidpTemp, d2lnphidt2Temp, d2lnphidtdpTemp, d2lnphidp2Temp,
        vTemp, dvdtTemp, dvdpTemp, d2vdt2Temp, d2vdtdpTemp, d2vdp2Temp, d3lnphidt3;

        dewH2ODriver(t*(1.0+sqrt(DBL_EPSILON)), p, &vTemp, &zTemp, &phiTemp, &dvdtTemp, &dvdpTemp, &d2vdt2Temp, &d2vdtdpTemp, &d2vdp2Temp,
                &dlnphidtTemp, &dlnphidpTemp, &d2lnphidt2Temp, &d2lnphidtdpTemp, &d2lnphidp2Temp);

        d3lnphidt3 = (d2lnphidt2Temp - d2lnphidt2)/t/sqrt(DBL_EPSILON);
        *dcpdt += -(2.0*R*dlnphidt + R*t*d2lnphidt2) -t*(3.0*R*d2lnphidt2 + R*t*d3lnphidt3);
    }
}

static void spline(double x[], double y[], int n, double yp1, double ypn, double y2[], double u[]) {
    int i,k;
    double p, qn, sig, un;

    if (yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }

    for (i=1; i<(n-1); i++) {
        sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p     = sig*y2[i-1] + 2.0;
        y2[i] = (sig-1.0)/p;
        u[i]  = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]  = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    if (ypn > 0.99e30)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

    for (k=n-2;k>=0;k--) y2[k] = y2[k]*y2[k+1]+u[k];
}

static double splint(double xa[], double ya[], double y2a[], int n, double x) {
    int klo, khi, k;
    double h, b, a;

    klo = 0;
    khi = n-1;
    while (khi-klo > 1) {
        k = (khi+klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo = k;
    }

    h = xa[khi] - xa[klo];
    if (h == 0.0) printf("Internal error in spline function.");

    a = (xa[khi]-x)/h;
    b = (x-xa[klo])/h;

    return a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(h*h)/6.0;
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
    static double xSpline[nZandDcorrections];
    static double uSpline[nZandDcorrections];
    static double ySplineForG[nZandDcorrections];
    static double ySplineForH[nZandDcorrections];
    static double ySplineForS[nZandDcorrections];
    static double ySplineForCp[nZandDcorrections];
    static double ySplineForDCpDt[nZandDcorrections];
    static double y2SplineForG[nZandDcorrections];
    static double y2SplineForH[nZandDcorrections];
    static double y2SplineForS[nZandDcorrections];
    static double y2SplineForCp[nZandDcorrections];
    static double y2SplineForDCpDt[nZandDcorrections];
    static int FIRST = 0;
    double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
    if (!FIRST) {
        for (int i=0; i<nZandDcorrections; i++) {
            xSpline[i]         = zAndDcorrections[i].t;
            ySplineForG[i]     = zAndDcorrections[i].g;
            ySplineForH[i]     = zAndDcorrections[i].h;
            ySplineForS[i]     = zAndDcorrections[i].s;
            ySplineForCp[i]    = zAndDcorrections[i].cp;
            ySplineForDCpDt[i] = zAndDcorrections[i].dcpdt;
        }
        spline(xSpline, ySplineForG,     nZandDcorrections, 1.0e30, 1.0e30, y2SplineForG,     uSpline);
        spline(xSpline, ySplineForH,     nZandDcorrections, 1.0e30, 1.0e30, y2SplineForH,     uSpline);
        spline(xSpline, ySplineForS,     nZandDcorrections, 1.0e30, 1.0e30, y2SplineForS,     uSpline);
        spline(xSpline, ySplineForCp,    nZandDcorrections, 1.0e30, 1.0e30, y2SplineForCp,    uSpline);
        spline(xSpline, ySplineForDCpDt, nZandDcorrections, 1.0e30, 1.0e30, y2SplineForDCpDt, uSpline);
        FIRST = 1;
    }
    double result = 0.0;

    propertiesOfPureH2O(t, p, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

    switch (returnMode) {
        case 1:
            result = gH2O - splint(xSpline, ySplineForG, y2SplineForG, nZandDcorrections, t);
            break;
        case 2:
            result = hH2O - splint(xSpline, ySplineForH, y2SplineForH, nZandDcorrections, t);
            break;
        case 3:
            result = sH2O - splint(xSpline, ySplineForS, y2SplineForS, nZandDcorrections, t);
            break;
        case 4:
            result = cpH2O - splint(xSpline, ySplineForCp, y2SplineForCp, nZandDcorrections, t);
            break;
        case 5:
            result = dcpdtH2O  - splint(xSpline, ySplineForDCpDt, y2SplineForDCpDt, nZandDcorrections, t);
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

double ZhangAndDuan2005_getGibbsFreeEnergy(double t, double p) {
    return calculate(t, p, returnValueOfG);
}

double ZhangAndDuan2005_getEnthalpy(double t, double p) {
    return calculate(t, p, returnValueOfH);
}

double ZhangAndDuan2005_getEntropy(double t, double p) {
    return calculate(t, p, returnValueOfS);
}

double ZhangAndDuan2005_getHeatCapacity(double t, double p) {
    return calculate(t, p, returnValueOfCP);
}

double ZhangAndDuan2005_getDcpDt(double t, double p) {
    return calculate(t, p, returnValueOfDCPDT);
}

double ZhangAndDuan2005_getVolume(double t, double p) {
    return  calculate(t, p, returnValueOfV);
}

double ZhangAndDuan2005_getDvDt(double t, double p) {
    return calculate(t, p, returnValueOfdVdT);
}

double ZhangAndDuan2005_getDvDp(double t, double p) {
    return  calculate(t, p, returnValueOfdVdP);
}

double ZhangAndDuan2005_getD2vDt2(double t, double p) {
    return calculate(t, p, returnValueOfd2VdT2);
}

double ZhangAndDuan2005_getD2vDtDp(double t, double p) {
    return calculate(t, p, returnValueOfd2VdTdP);
}

double ZhangAndDuan2005_getD2vDp2(double t, double p) {
    return calculate(t, p, returnValueOfd2VdP2);
}
