#ifndef TORQUE2_RING_H
#define TORQUE2_RING_H

#define Pi 3.14159265358979


double Ab(double a, double e, double r1x, double r1y, double r1z, double b);
double Bce(double a, double e, double r1x, double r1y, double r1z, double b);
double Bse(double a, double e, double r1x, double r1y, double r1z, double b);
double B(double a, double e, double r1x, double r1y, double r1z, double b);
double CC(double a, double e);
double Q(double a, double e, double r1x, double r1y, double r1z, double b);
double R(double a, double e, double r1x, double r1y, double r1z, double b);
double theta(double a, double e, double r1x, double r1y, double r1z, double b);
double* lambda(double a, double e, double r1x, double r1y, double r1z, double b);
void Qmat(double a, double e, double r1x, double r1y, double r1z, double b, double QQ[3][3]);
double* UU(double a, double e, double r1x, double r1y, double r1z, double b);
double* VV(double a, double e, double r1x, double r1y, double r1z, double b);
void Fmat(double a, double e, double r1x, double r1y, double r1z, double b, double FF[3][3]);
double* Fu(double a, double e, double r1x, double r1y, double r1z, double b);
double* Fv(double a, double e, double r1x, double r1y, double r1z, double b);
double* force(double a, double e, double r1x, double r1y, double r1z, double b);
double* torque(double a, double e, double r1x, double r1y, double r1z, double b);

#endif
