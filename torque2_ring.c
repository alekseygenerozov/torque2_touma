#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_ellint.h>
#include "torque2_ring.h"

#define Pi 3.14159265358979


double Ab(double a, double e, double r1x, double r1y, double r1z, double b){
	//Now with compact softening!!!
	double nu1=atan2(r1y, r1x);
	double x[3];
	x[0] = a*(1.0-e*e)*cos(nu1)/(1.0+e*cos(nu1));
	x[1] = a*(1.0-e*e)*sin(nu1)/(1.0+e*cos(nu1));
	x[2] = 0;

 	double dist = pow((r1x-x[0])*(r1x-x[0])+(r1y-x[1])*(r1y-x[1])+r1z*r1z, 0.5);
 	double b1=(dist<b)?b:0;

	return r1x*r1x+r1y*r1y+r1z*r1z+a*a+b1*b1+2.0*a*e*r1x;

}
double Bce(double a, double e, double r1x, double r1y, double r1z, double b){
	return a*r1x+a*a*e;
}

double Bse(double a, double e, double r1x, double r1y, double r1z, double b){
	return a*pow(1.0-e*e,0.5)*r1y;
}

double B(double a, double e, double r1x, double r1y, double r1z, double b){
	double b1=Bse(a, e, r1x, r1y, r1z, b);
	double b2=Bce(a, e, r1x, r1y, r1z, b);
	return pow(b1*b1+b2*b2, 0.5);
}


double CC(double a, double e){
	return a*a*e*e;
}

double Q(double a, double e, double r1x, double r1y, double r1z, double b){
	double Ab1=Ab(a, e, r1x, r1y, r1z, b);
	double B1=B(a, e, r1x, r1y, r1z, b);
	double CC1=CC(a, e);
	return (1.0/9.0)*(CC1-Ab1)*(CC1-Ab1)-(1.0/3.0)*(B1*B1-Ab1*CC1);
}

double R(double a, double e, double r1x, double r1y, double r1z, double b){
	double Ab1=Ab(a, e, r1x, r1y, r1z, b);
	double B1=B(a, e, r1x, r1y, r1z, b);
	double CC1=CC(a, e);
	double Bse1=Bse(a, e, r1x, r1y, r1z, b);

	return (1.0/27.0)*pow(CC1-Ab1, 3)-(1.0/6.0)*(CC1-Ab1)*(B1*B1-Ab1*CC1)+0.5*CC1*Bse1*Bse1;
}

double theta(double a, double e, double r1x, double r1y, double r1z, double b){
	return acos(R(a, e, r1x, r1y, r1z, b)/pow(Q(a, e, r1x, r1y, r1z, b), 1.5));
}

double* lambda(double a, double e, double r1x, double r1y, double r1z, double b){
	static double lam[3];

	double Ab1=Ab(a, e, r1x, r1y, r1z, b);
	double CC1=CC(a,e);
	double th=theta(a, e, r1x, r1y, r1z, b);
	double Q1=Q(a, e, r1x, r1y, r1z, b);

	lam[0]=-2.0*pow(Q1,0.5)*cos((1.0/3.0)*th+(2.0/3.0)*Pi)-(1.0/3.0)*(CC1-Ab1);
	lam[1]=-2.0*pow(Q1,0.5)*cos((1.0/3.0)*th-(2.0/3.0)*Pi)-(1.0/3.0)*(CC1-Ab1);
	lam[2]=-2.0*pow(Q1,0.5)*cos((1.0/3.0)*th)-(1.0/3.0)*(CC1-Ab1);

	return lam;
}

void Qmat(double a, double e, double r1x, double r1y, double r1z, double b, double QQ[3][3]){
	double* lams=lambda(a, e, r1x, r1y, r1z, b);
	double lam0=lams[0];
	double lam1=lams[1];
	double lam2=lams[2];
	double CC1=CC(a, e);
	double Bse1=Bse(a, e, r1x, r1y, r1z, b);
	double Bce1=Bce(a, e, r1x, r1y, r1z, b);

	QQ[0][0]=sqrt((lam0*(CC1 + lam0))/((lam0 - lam1)*(lam0 - lam2)));
	QQ[0][1]=sqrt((lam1*(CC1 + lam1))/((lam0 - lam1)*(lam1 - lam2)));
	QQ[0][2]=sqrt((CC1 + lam2)/((lam0 - lam2)*(lam1 - lam2)))*sqrt(fabs(lam2));
	QQ[1][0]=Bse1*sqrt((CC1 + lam0)/(lam0*(lam0 - lam1)*(lam0 - lam2)));
	QQ[1][1]=Bse1*sqrt((CC1 + lam1)/((lam0 - lam1)*lam1*(lam1 - lam2)));
	QQ[1][2]=-((Bse1*sqrt((CC1 + lam2)/((lam0 - lam2)*(lam1 - lam2))))/sqrt(fabs(lam2)));
	QQ[2][0]=Bce1*sqrt(lam0/((CC1+ lam0)*(lam0 - lam1)*(lam0 - lam2)));
	QQ[2][1]=Bce1*sqrt(lam1/((lam0 - lam1)*(CC1 + lam1)*(lam1 - lam2)));
	QQ[2][2]=Bce1*sqrt(1/((lam0 - lam2)*(lam1 - lam2)*(CC1 + lam2)))*sqrt(fabs(lam2));


}

double* UU(double a, double e, double r1x, double r1y, double r1z, double b){
	static double UU1[3];
	double QQ[3][3];
	Qmat(a, e, r1x, r1y, r1z, b, QQ);
	int i=0;
	for (i=0; i<3; i++){
		UU1[i]=QQ[0][0]*QQ[i][0]-e*QQ[2][0]*QQ[i][0]+QQ[0][2]*QQ[i][2]-e*QQ[2][2]*QQ[i][2];
	}

	return UU1;
}

double* VV(double a, double e, double r1x, double r1y, double r1z, double b){
	static double VV1[3];
	double QQ[3][3];
	Qmat(a, e, r1x, r1y, r1z, b, QQ);
	int i=0;
	for (i=0; i<3; i++){
		VV1[i]=QQ[0][1]*QQ[i][1]-e*QQ[2][1]*QQ[i][1]-QQ[0][2]*QQ[i][2]+e*QQ[2][2]*QQ[i][2];
	}

	return VV1;
}

void Fmat(double a, double e, double r1x, double r1y, double r1z, double b, double FF[3][3]){
	FF[0][0]=-r1x-a*e;
	FF[0][1]=-r1y;
	FF[0][2]=-r1z;
	FF[1][0]=0;
	FF[1][1]=a*pow(1.0-e*e, 0.5);
	FF[1][2]=0;
	FF[2][0]=a;
	FF[2][1]=0;
	FF[2][2]=0;

}


double* Fu(double a, double e, double r1x, double r1y, double r1z, double b){
	static double Fu1[3];
	double FF[3][3];
	double* UU1=UU(a, e, r1x, r1y, r1z, b);
	Fmat(a, e, r1x, r1y, r1z, b, FF);
	int i=0;
	int j=0;

	for (i=0; i<3; i++){
		Fu1[i]=0;
		for (j=0; j<3; j++){
			Fu1[i]+=FF[j][i]*UU1[j];
		}
	}
	return Fu1;
}


double* Fv(double a, double e, double r1x, double r1y, double r1z, double b){
	static double Fv1[3];
	double FF[3][3];
	double* VV1=VV(a, e, r1x, r1y, r1z, b);
	Fmat(a, e, r1x, r1y, r1z, b, FF);
	int i=0;
	int j=0;

	for (i=0; i<3; i++){
		Fv1[i]=0;
		for (j=0; j<3; j++){
			Fv1[i]+=FF[j][i]*VV1[j];
		}
	}
	return Fv1;
}

double* force(double a, double e, double r1x, double r1y, double r1z, double b){
	static double force1[3];
	double* Fu1=Fu(a, e, r1x, r1y, r1z, b);
	double* Fv1=Fv(a, e, r1x, r1y, r1z, b);
	double* lams=lambda(a, e, r1x, r1y, r1z, b);
	double lam0=lams[0];
	double lam1=lams[1];
	double lam2=lams[2];
	double k=sqrt((lam1-lam2)/(lam0-lam2));
	int i=0;
	for (i=0; i<3; i++){
		force1[i]=(2.0/Pi)*sqrt(lam0-lam2)/((lam0-lam1)*(lam1-lam2))*((k*k*Fu1[i]+Fv1[i])*gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE)-(1.0-k*k)*Fv1[i]*gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE));
	}
	return force1;
}

double* torque(double a, double e, double r1x, double r1y, double r1z, double b){
	static double torque1[3];
	double* force1=force(a, e, r1x, r1y, r1z, b);
	
	torque1[0]=r1y*force1[2]-r1z*force1[1];
	torque1[1]=-(r1x*force1[2]-r1z*force1[0]);
	torque1[2]=r1x*force1[1]-r1y*force1[0];
	return torque1;
}

// int main(int argc, char* argv[]){

// 	static double QQ[3][3];
// 	double* test=torque(1.0, 0.7, 0.1, 0.2, 0.3, 1e-5);
// 	printf("%f %f %f\n", test[0], test[1], test[2]);

// 	return 0;
// }