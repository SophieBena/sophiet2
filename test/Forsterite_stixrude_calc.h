
const char *Forsterite_stixrude_identifier(void);
const char *Forsterite_stixrude_name(void);
const char *Forsterite_stixrude_formula(void);
const double Forsterite_stixrude_mw(void);
const double *Forsterite_stixrude_elements(void);

double Forsterite_stixrude_g(double T, double P);
double Forsterite_stixrude_dgdt(double T, double P);
double Forsterite_stixrude_dgdp(double T, double P);
double Forsterite_stixrude_d2gdt2(double T, double P);
double Forsterite_stixrude_d2gdtdp(double T, double P);
double Forsterite_stixrude_d2gdp2(double T, double P);
double Forsterite_stixrude_d3gdt3(double T, double P);
double Forsterite_stixrude_d3gdt2dp(double T, double P);
double Forsterite_stixrude_d3gdtdp2(double T, double P);
double Forsterite_stixrude_d3gdp3(double T, double P);

double Forsterite_stixrude_s(double T, double P);
double Forsterite_stixrude_v(double T, double P);
double Forsterite_stixrude_cv(double T, double P);
double Forsterite_stixrude_cp(double T, double P);
double Forsterite_stixrude_dcpdt(double T, double P);
double Forsterite_stixrude_alpha(double T, double P);
double Forsterite_stixrude_beta(double T, double P);
double Forsterite_stixrude_K(double T, double P);
double Forsterite_stixrude_Kp(double T, double P);

