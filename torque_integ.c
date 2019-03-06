
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "torque2_ring.h"
#include "rebound.h"



double sig(double a, double a1, double a2, double p){
	double No=(1.0-p)/(pow((a2/a1), 1.0-p)-1.0);
	return No*pow((a/a1), -p);
	// return a2;
}

double* force_tot_point(double a1, double a2, double p, double e0, double q, double r1x, double r1y, double r1z, double b){
    // Central body
    //Number of sma bins
    int N0=100;
    int j=0;
    int k=0;
    double deltaA=(a2-a1)/N0;

    double* force1;
    static double forceTot[3];
    for (k=0; k<3; k++){
        forceTot[k]=0;
    }

    double a=a1;
    double e=e0;
    for (j=0; j<N0; j++){
            force1=force(a, e, r1x, r1y, r1z, b);
            for (k=0; k<3; k++){
                forceTot[k]+=sig(a, a1, a2, p)*force1[k]*deltaA;            
            }
    a+=deltaA;
    e=pow((a/a1), q)*e0;

    }

    return forceTot;
}

double* prec_tot(double a1, double a2, double p, double e0, double q, double e_test, double a_test, double omega_test, double b, int N, int flag, 
char* tag){
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
    double torqueOrb[N][3];
    double forceOrb[N][3];
    double ieDotOrb[N];
    double ieDotOrb2[N];
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
    FILE* forceFile;

    char tag2[80]="forces_";
    strcat(tag2, tag);
    forceFile=fopen(tag2, "w");
    for (j=0; j<N0; j++){
	    for (i=0; i<N; i++){
            if (j==0){
                ieDotOrb[i]=0;
                ieDotOrb2[i]=0;
                for (k=0; k<3; k++){
                    torqueOrb[i][k]=0;
                    forceOrb[i][k]=0;
                }
            }
	    	f=reb_tools_M_to_f(e_test, MM);
	    	test=reb_tools_orbit2d_to_particle(1.0, star, mm, a_test, e_test, omega_test, f);
            MM+=deltaMM;
	    	//test=reb_tools_orbit_to_particle(1.0, star, mm, a_test, e_test, 0.0, omega_test, 0.0, f);
	    	force1=force(a, e, test.x, test.y, test.z, b);
	    	torque1=torque(a, e, test.x, test.y, test.z, b);

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
                torqueOrb[i][k]+=sig(a, a1, a2, p)*torque1[k]*deltaA*deltaMM/2.0/Pi;
                forceOrb[i][k]+=sig(a, a1, a2, p)*force1[k]*deltaA*deltaMM/2.0/Pi;
                ieDotOrb[i]+=(ex*edoty-ey*edotx)/(e_test*e_test)*sig(a, a1, a2, p)*deltaA*deltaMM/2.0/Pi;
                ieDotOrb2[i]+=(ex*edoty2-ey*edotx2)/(e_test*e_test)*sig(a, a1, a2, p)*deltaA*deltaMM/2.0/Pi;
            }
	    	
    	}
    //printf("%lf %lf %lf\n", star.x, star.y, star.z);
    a+=deltaA;
    e=pow((a/a1), q)*e0;

	}
    fprintf(forceFile, "x y z vx vy vz tau fx fy ie1 ie2\n");
    MM=1.1e-4+flag*deltaMM*0.5;
    for (i=0; i<N; i++){
        f=reb_tools_M_to_f(e_test, MM);
        test=reb_tools_orbit2d_to_particle(1.0, star, mm, a_test, e_test, omega_test, f);
        MM+=deltaMM;
        fprintf(forceFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", test.x, test.y, test.z, test.vx, test.vy, test.vz, torqueOrb[i][2], 
            forceOrb[i][0], forceOrb[i][1], ieDotOrb[i], ieDotOrb2[i]);
    }
    fclose(forceFile);

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
    char* e_in="0.7";
    char* a_test="0.99"; 
    char* ang_test="0.0";
    char* b="1e-4";
    char* flag="0";
    char* points="1001";
    char* q_disk="0.0";


    static struct option long_options[] =
        {
          /* These options set a flag. */
          {"ein", required_argument, NULL, 'i'}, 
          {"atest", required_argument, NULL, 'a'},
          {"etest", required_argument, NULL, 'e'},
          {"pomega", required_argument, NULL, 'o'},
          {"flag",    required_argument, NULL, 'f'},
          {0, 0, 0, 0}
        };

    while((opt = getopt_long(argc, argv, "e:a:o:n:f:q:b:", long_options, NULL)) != -1)  
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
            case 'i':
                e_in=(optarg);
                break;
            case 'o':
                ang_test=(optarg);
                break;
            case 'n':
                points=(optarg);
                break;
            case 'f':
                flag=(optarg);
                if (atoi(flag)!=0)
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
    char tag2[80]="";
    strcat(tag2, tag);
    printf("%s\n", tag2);

    // int i=1;
    strcat(tag2, "_e");
    strcat(tag2, e_test);
    strcat(tag2, "_a");
    strcat(tag2, a_test);
    strcat(tag2, "_o");
    strcat(tag2, ang_test);
    strcat(tag2, "_q");
    strcat(tag2, q_disk);
    strcat(tag2, "_ein");
    strcat(tag2, e_in);

	double* sol=prec_tot(1.0, 2.0, 1.01, atof(e_in), atof(q_disk), atof(e_test), atof(a_test),
        atof(ang_test)*Pi/180.0, atof(b), atoi(points), atoi(flag), tag2);



    
    printf("%s\n",tag2);
    out("tau", tag2, norm*sol[0]);
    out("i", tag2, norm*sol[1]);
    out("i2", tag2, norm*sol[2]);
    out("idd", tag2, norm2*sol[3]);

	return 0;


}