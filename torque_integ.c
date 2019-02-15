
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <system.h>
// #include <fenv.h>

#include "torque2_ring.h"
#include "rebound.h"




double sig(double a, double a1, double a2, double p){
	double No=(1.0-p)/(pow((a2/a1), 1.0-p)-1.0);
	return No*pow((a/a1), -p);
	// return a2;
}

double* prec_tot(double a1, double a2, double p, double e0, double q, double e_test, double a_test, double omega_test, double b, int N, int flag){
	// Central body
    static double sol[2];
    struct reb_particle star = {0};
    struct reb_simulation* r = reb_create_simulation();
    star.m = 1;
    reb_add(r, star);
    //Number 
    int N0=100;
    int i=0;
    int j=0;
    struct reb_particle test;
    double deltaMM=2.0*Pi/N;
    double MM=1.1e-4+flag*deltaMM*0.5;
    double deltaA=(a2-a1)/N0;

   	double* torque1;
   	double* force1;
   	double torqueTot=0;
   
   	double ex=e_test*cos(omega_test);
   	double ey=e_test*sin(omega_test);
    double jz=pow(a_test*(1.0-e_test*e_test), 0.5);
    
    double a=a1;
    double e=e0;
    double f=0;
    double edotx=0;
    double edoty=0;
    double edot=0;
    double mm=1.0e-30;
    //FILE* forceFile;
    //forceFile=fopen("forces", "w");

    for (j=0; j<N0; j++){
	    for (i=0; i<N; i++){
	    	f=reb_tools_M_to_f(e_test, MM);
	    	//test=reb_tools_orbit2d_to_particle(1.0, star, 2.5e-7, a_test, e_test, omega_test, f);
	    	test=reb_tools_orbit_to_particle(1.0, star, mm, a_test, e_test, 0.0, omega_test, 0.0, f);
	    	force1=force(a, e, test.x, test.y, test.z, b);
	    	//fprintf(forceFile, "%lf %lf %lf %lf %lf %lf\n", test.x, test.y, test.z, force1[0], force1[1], force1[2]);

	    	torque1=torque(a, e, test.x, test.y, test.z, b);

	    	MM+=deltaMM;
	    	edotx=force1[1]*jz;//+test.vy*torque1[2]-test.vz*torque1[1];
	    	edoty=-force1[0]*jz;//-test.vx*torque1[2]-test.vz*torque1[0];
	    	edot+=(ex*edoty-ey*edotx)/(e_test*e_test)*sig(a, a1, a2, p)*deltaA*deltaMM/2.0/Pi;
	    	torqueTot+=sig(a, a1, a2, p)*torque1[2]*deltaA*deltaMM/2.0/Pi;
    	}
    //printf("%lf %lf %lf\n", star.x, star.y, star.z);
    a+=deltaA;
    e=pow((a/a1), q)*e0;

	}
	sol[0]=torqueTot;
	sol[1]=edot;
	return sol;
}



int main(int argc, char* argv[]){
    // feenableexcept(FE_INVALID | FE_OVERFLOW);

	char* tag="a";
	if (atoi(argv[6])){
		tag="b";
	}

	double m = 2.5e-7;
	double norm=1000.0*pow(m,2.0);
	double* sol=prec_tot(1.0, 2.0, 1.01, 0.7, 0.0, atof(argv[1]), atof(argv[2]), atof(argv[3])*Pi/180.0, atof(argv[4]), atoi(argv[5]), atoi(argv[6]) );
    printf("%d \n", atoi(argv[6]));

	char fname[80]="";
    strcat(fname, "tau");
    strcat(fname,"_");
    strcat(fname, tag);
    int i=1;
    for (i=1; i<4; i++){
        strcat(fname, "_");
        strcat(fname, argv[i]);
    }
    printf("%s\n",fname);
    FILE *f=fopen(fname, "w");
	fprintf(f, "%0.10e\n", 
		norm*sol[0]	
	);
	fclose(f);

	char iname[80]="";
    strcat(iname, "i");
    strcat(iname,"_");
    strcat(iname, tag);
    for (int i=1; i<4; i++){
        strcat(iname, "_");
        strcat(iname, argv[i]);
    }
    printf("%s\n",iname);
    FILE *ie=fopen(iname, "w");
    fprintf(f, "%0.10e\n", 
		norm*sol[1]
	);
    fclose(ie);



	return 0;


}