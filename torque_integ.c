
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <confuse.h>
// #include <system.h>
// #include <fenv.h>
#include <unistd.h>


#include "torque2_ring.h"
#include "rebound.h"


// static struct argp_option options[] = { 
//     { "etest", 'e', 0, 0, "Compare lines instead of characters."},
//     { "word", 'w', 0, 0, "Compare words instead of characters."},
//     { "nocase", 'i', 0, 0, "Compare case insensitive instead of case sensitive."},
//     { 0 } 
// };




double sig(double a, double a1, double a2, double p){
	double No=(1.0-p)/(pow((a2/a1), 1.0-p)-1.0);
	return No*pow((a/a1), -p);
	// return a2;
}

double* prec_tot(double a1, double a2, double p, double e0, double q, double e_test, double a_test, double omega_test, double b, int N, int flag){
	// Central body
    struct reb_particle star = {0};
    star.m = 1;
    static double sol[4];
    //Number 
    int N0=100;
    int i=0;
    int j=0;
    int k=0;
    struct reb_particle test;
    double deltaMM=2.0*Pi/N;
    double MM=1.1e-4+flag*deltaMM*0.5;
    double deltaA=(a2-a1)/N0;

   	double* torque1;
   	double* force1;
   	double torqueTot[3];
    double forceTot[3];
    for (k=0; k<3; k++){
        forceTot[k]=0;
        torqueTot[k]=0;
    }
   
   	double ex=e_test*cos(omega_test);
   	double ey=e_test*sin(omega_test);
    double jz=pow(a_test*(1.0-e_test*e_test), 0.5);
    
    double a=a1;
    double e=e0;
    double f=0;
    double edotx=0;
    double edoty=0;
    double iedot=0;

    double edotx2=0;
    double edoty2=0;
    double edot2=0;
    double mm=1.0e-30;

    double eddotx=0;
    double eddoty=0;
    double ieddot=0;
    //FILE* forceFile;
    //forceFile=fopen("forces", "w");

    for (j=0; j<N0; j++){
	    for (i=0; i<N; i++){
	    	f=reb_tools_M_to_f(e_test, MM);
	    	test=reb_tools_orbit2d_to_particle(1.0, star, mm, a_test, e_test, omega_test, f);
	    	//test=reb_tools_orbit_to_particle(1.0, star, mm, a_test, e_test, 0.0, omega_test, 0.0, f);
	    	force1=force(a, e, test.x, test.y, test.z, b);
	    	torque1=torque(a, e, test.x, test.y, test.z, b);

	    	MM+=deltaMM;
	    	edotx=force1[1]*jz+(test.vy*torque1[2]-test.vz*torque1[1]);
	    	edoty=-force1[0]*jz-(test.vx*torque1[2]-test.vz*torque1[0]);
	    	iedot+=(ex*edoty-ey*edotx)/(e_test*e_test)*sig(a, a1, a2, p)*deltaA*deltaMM/2.0/Pi;

            edotx2=force1[1]*jz;
            edoty2=-force1[0]*jz;
            edot2+=(ex*edoty2-ey*edotx2)/(e_test*e_test)*sig(a, a1, a2, p)*deltaA*deltaMM/2.0/Pi;

            eddotx=2.0*(force1[1]*torque1[2]-force1[2]*torque1[1]);
            eddoty=-2.0*(force1[0]*torque1[2]-force1[2]*torque1[0]);
            ieddot+=(ex*eddoty-ey*eddotx)/(e_test*e_test)-2.0*(ex*edoty-ey*edotx)/(e_test*e_test*e_test);

            for (k=0; k<3; k++){
                torqueTot[k]+=sig(a, a1, a2, p)*torque1[k]*deltaA*deltaMM/2.0/Pi;
                forceTot[k]+=sig(a, a1, a2, p)*force1[k]*deltaA*deltaMM/2.0/Pi;
            }
	    	
    	}
    //printf("%lf %lf %lf\n", star.x, star.y, star.z);
    a+=deltaA;
    e=pow((a/a1), q)*e0;

	}
	sol[0]=torqueTot[2];
	sol[1]=iedot;
    sol[2]=edot2;
    sol[3]=ieddot;
	return sol;
}


void out(char* pre, char* tag, double out){
    char fname[80]="";
    strcat(fname, pre);
    strcat(fname, "_");
    strcat(fname, tag);
    FILE* f = fopen(fname, "w");
    fprintf(f, "%0.10e\n", out);
    fclose(f);

}

int main(int argc, char* argv[]){
    // feenableexcept(FE_INVALID | FE_OVERFLOW);

    // put ':' in the starting of the 
    // string so that program can  
    //distinguish between '?' and ':'
    // extern char *optarg; 
    // extern int optind;
    int opt; 
    char* tag="a";
    char* e_test="0.7";
    char* a_test="0.99"; 
    char* ang_test="0.0";
    char* b="1e-4";
    char* flag="0";
    char* points="1001";
    char* q_disk="0.0";

    while((opt = getopt(argc, argv, "e:a:o:n:f:q:b:")) != -1)  
    {  
        // printf("%d \n", opt);
        switch(opt)  
        {  
            case 'e':  
                e_test=(optarg);
                break;
            case 'a': 
                a_test=(optarg);
                break;
            case 'o':
                ang_test=(optarg);
                break;
            case 'n':
                points=(optarg);
                break;
            case 'f':
                flag=(optarg);
                if (flag)
                    tag="b";
                break;
            case 'q':
                q_disk=(optarg);
                break;
            case 'b':
                b=(optarg);
                break;
        }  
    }  
    printf("%s\n", e_test);
	double m = 2.5e-7;
	double norm=1000.0*pow(m,2.0);
    double norm2=pow(1000.0*m,2.0);
	double* sol=prec_tot(1.0, 2.0, 1.01, 0.7, atof(q_disk), atof(e_test), atof(a_test), atof(ang_test)*Pi/180.0, atof(b), atoi(points), atoi(flag));

	char tag2[80]="";

    strcat(tag2, tag);
    printf("%s\n", tag2);

    // int i=1;
    strcat(tag2, "_");
    strcat(tag2, e_test);
    strcat(tag2, "_");
    strcat(tag2, a_test);
    strcat(tag2, "_");
    strcat(tag2, ang_test);

    
    printf("%s\n",tag2);
    out("tau", tag2, norm*sol[0]);
    out("i", tag2, norm*sol[1]);
    out("i2", tag2, norm*sol[2]);
    out("idd", tag2, norm2*sol[3]);

	return 0;


}