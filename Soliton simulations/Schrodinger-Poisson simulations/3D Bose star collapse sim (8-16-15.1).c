/*
Author: Eric Cotner, UCLA Physics 2015

This program solves the time-dependent Schrödinger-Newton equation

i\dot{\psi} = -\frac{1}{2m} \nabla^2 \psi + \frac{\lambda}{8m^2} |\psi|^2 \psi + m \phi \psi

which governs the evolution of a self-gravitating Bose-Einstein condensate. Since the condensate in question is formed of axions,
there is a self-interaction term which adds additional nonlinearity to the problem. The algorithm used to solve this equation
is divided into two parts. The first step is to solve for the gravitational potential \phi at each time step using the Poisson
equation

\nabla^2 \phi = 4\pi G m |\psi|^2

This is accomplished by using a successive over-relaxation (SOR) method, which relaxes the solution from an initial guess to
a solution which can be arbitrarily close to the exact answer depending on how long you're willing to wait for an answer. In order
enforce physical boundary conditions (zero potential at infinity) on a finite grid, an asymptotic coordinate transformation is
performed

\xi = \tanh(x/R)

which maps x \in (-\infty,\infty) to the finite domain \xi \in [-1,1]. After transforming all spatial dimensions the domain over
which this equation is solved goes from R^3 -> [-1,1]^3, and can now be discretized and boundary conditions at "infinity" can
be imposed easily.
Once \phi has been solved for at the present timestep, it is then plugged into the Newton-Schrödinger equation. This is solved
using the usual time-dependent grid method, where each new timestep is updated using values from the previous timesteps. Due to
stability issues, the splitting of the time derivative is center-difference, so the first *two* timesteps need to be specified
in order to move forward. The solution for \psi also utilizes an asymptotic coordinate system to enforce physical boundary
conditions (\psi = 0 at infinity). There is also an additional speedup on timesteps past the first due to being able to use the
solution to \phi from the previous timestep as the initial guess for the next timestep.

Once the simulation has run its course, the data for \psi and \phi are written to a text file
which can then be imported to Mathematica and analyzed.

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
//#include <windows.h>  // only need this to use the Beep() function, it's totally optional
#include "./custom.c"

//#define M 30    // number of grid points
#define T 1000  // number of timesteps
#define TRmax 200   // relaxation timesteps
//#define DX 2./M // grid size
//#define DT 1.0e+94  // timestep size (1 GeV^-1 = 6.58e-25 s)
//#define R 3.0e+56   // asymptotic sensitivity parameter (1 GeV^-1 = 2e-19 km)
//* remember stability condition is 2*d*DT/(m*(R*DX)^2) < 1 where d is the number of dimensions *//
#define G 6.708e-39 // Newton's gravitational constant 6.67e-11 in SI units, 6.708e-39 in GeV^-2
//* should be stable for all values, but G*m^2*DT*DX^2*M^2*psi0^2 << 1 is recommended *//
#define w 1.9   // SOR relaxation parameter 1<w<2
#define PI 3.14159265358979

///* GLOBAL VARIABLES */
double a;    // self-coupling constant; negative is attractive, positive is repulsive
//* should be stable for all values, but a*DT*psi0^2/(4m^2) << 1 *//
double m; // axion mass in GeV
double R, DT, DX, N;
int M;


///* CUSTOM FUNCTION DECLARATION *//
// these are the asymptotic coordinates, they're just labeled as the normal coordinates since \xi, \eta and \zeta are not
// supported text characters
double x(int i);
double y(int j);
double z(int k);

int main(){
///* LOCAL VARIABLES */
    int i1=0, i2=0, i3=0, n=0;  // the iterators! i's are spatial, n is time
    int nR=0, n0=0, n1=1, n2=2;   // nR is relaxation for \phi, ni are the timestep indicators for the fields
    double acc_param, acc_goal=1.0e-4;

///* FILE I/0 SETUP */
    FILE *fpInit;
    FILE *fpOutPsi;
    FILE *fpOutPhi;
    fpInit = fopen("./Data files/Initial states/(8-16-15)PythonTestInitial.txt","r");   // initial conditions
    fpOutPsi = fopen("./Python scripts/(8-16-15.1)PythonTestPsi.txt","w+");    // data output file for psi
//    fpOutPhi = fopen("./Data files/Stable state parameter exploration/(8-13-15.1)IntermediateAttractivePhi.txt","w+");    // data output file for phi
    int minimize_data = 0;  // if 1, truncates the number of timesteps to write to file (to reduce file size)
    int timesteps_to_keep = 400;    // number of timesteps to write to file

    // sets up data minimization
    int data_skip;
    if(minimize_data)
        data_skip = T/timesteps_to_keep;  // number of timesteps to skip
    else
        data_skip = 1;

///* SET PARAMETERS FROM FILE */
    //fscanf(fpInit, "R=%lf, DT=%lf, M=%i, m=%lf, a=%lf, N=%e\n\n", &R, &DT, &M, &m, &a, &N);
    fscanf(fpInit, "R=%lf, DT=%lf, M=%i, ", &R, &DT, &M);
    fscanf(fpInit, "m=%lf, a=%lf, N=%lf\n\n", &m, &a, &N);
    printf("R=%e, DT=%e, M=%i, m=%e, a=%e, N=%e\n", R, DT, M, m, a, N);
    // writes the parameters to the output files
    fprintf(fpOutPsi, "{R->");
    fprintf(fpOutPsi, "%.2e", R);
    fprintf(fpOutPsi, ", DT->");
    fprintf(fpOutPsi, "%.2e", DT);
    fprintf(fpOutPsi, ", m->");
    fprintf(fpOutPsi, "%.2e", m);
    fprintf(fpOutPsi, ", a->");
    fprintf(fpOutPsi, "%.2e", a);
    fprintf(fpOutPsi, ", n->");
    fprintf(fpOutPsi, "%.2e", N);
    fprintf(fpOutPsi, ", M->%i, T->%i}\n\n[", M, T);

//    fprintf(fpOutPhi, "{R->");
//    fprintM(fpOutPhi, R);
//    fprintf(fpOutPhi, ", DT->");
//    fprintM(fpOutPhi, DT);
//    fprintf(fpOutPhi, ", m->");
//    fprintM(fpOutPhi, m);
//    fprintf(fpOutPhi, ", a->");
//    fprintM(fpOutPhi, a);
//    fprintf(fpOutPhi, ", n->");
//    fprintM(fpOutPhi, N);
//    fprintf(fpOutPhi, ", M->%i, T->%i}\n\n{", M, T);

    DX=2.0/M;

///* DYNAMIC MEMORY ALLOCATION OF PSI AND PHI ARRAYS */
    double complex ****psi=NULL;
    double complex ****phi=NULL;
    double complex ****phiR=NULL;    // yo dawg i herd u liek pointers
    psi = malloc(M*sizeof(complex double));
    phi = malloc(M*sizeof(complex double));
    phiR = malloc(M*sizeof(complex double));
    for(i1=0; i1<=M; i1++){
        psi[i1] = malloc(M*sizeof(complex double));
        phi[i1] = malloc(M*sizeof(complex double));
        phiR[i1] = malloc(M*sizeof(complex double));
        for(i2=0; i2<=M; i2++){
            psi[i1][i2] = malloc(M*sizeof(complex double));
            phi[i1][i2] = malloc(M*sizeof(complex double));
            phiR[i1][i2] = malloc(M*sizeof(complex double));
            for(i3=0; i3<=M; i3++){
                psi[i1][i2][i3] = malloc(3*sizeof(complex double));
                phi[i1][i2][i3] = malloc(3*sizeof(complex double));
                phiR[i1][i2][i3] = malloc((TRmax+1)*sizeof(complex double));
            }
        }
    }


//    double complex psi[M+1][M+1][M+1][3];   // BEC wavefunction; only need three timesteps, one of which is overwritten
//    double phi[M+1][M+1][M+1][3];   // gravitational potential; only three are needed as above
//    double phiR[M+1][M+1][M+1][TRmax+1];   // relaxation steps for grav. potential


///* INSTABILITY WARNING */
    char inst_flag;
    double inst_param = 6*DT*pow(R*DX,-2)/m;
    if(inst_param >= 1.0){
        while((inst_flag != 'Y') || (inst_flag != 'y')){
            printf("Instability in time evolution detected - proceed? Y/N: ");
            //Beep(523,500);
            gets(&inst_flag);
            if((inst_flag == 'N') || (inst_flag == 'n')){
                printf("Aborting program\n");
                exit(1);
            }
            else if((inst_flag == 'Y') || (inst_flag == 'y'))
                break;
        }
    }
    else if(inst_param >= 0.5 && inst_param < 1.0){
        printf("Careful; instability parameter is %f\n", inst_param);
        //Beep(523,500);

    }


///* SET INITIAL DATA FOR PSI */
    double re, im;

    fseek(fpInit, 1, SEEK_CUR);
    for(i3=0; i3<=M; i3++){
        fseek(fpInit, 1, SEEK_CUR);
        for(i2=0; i2<=M; i2++){
            fseek(fpInit, 1, SEEK_CUR);
            for(i1=0; i1<=M-1; i1++){
                fscanf(fpInit, "{%lf, %lf}, ", &re, &im);
                psi[i1][i2][i3][0] = re + im * I;
                psi[i1][i2][i3][1] = re + im * I;
                //printf("%e + I*%e\n", creal(psi[i1][i2][i3][1]), cimag(psi[i1][i2][i3][1]));
                //printf("i1=%i, i2=%i, i3=%i\n", i1, i2, i3);
            }
            fscanf(fpInit, "{%lf, %lf}}", &re, &im);
            psi[M][i2][i3][0] = re + im * I;
            psi[M][i2][i3][1] = re + im * I;
            fseek(fpInit, 2, SEEK_CUR);
        }
        fseek(fpInit, 1, SEEK_CUR);
    }
    fseek(fpInit, 1, SEEK_CUR);
    fclose(fpInit);

///* BOUNDARY CONDITIONS */
    // psi boundary conditions
    for(n=0; n<=2; n++){
        for(i1=0; i1<=M; i1++){
            for(i2=0; i2<=M; i2++){
                // wavefunction BC at "infinity"
                psi[i1][i2][0][n] = 0. + 0.*I;
                psi[i1][i2][M][n] = 0. + 0.*I;
                psi[i1][0][i2][n] = 0. + 0.*I;
                psi[i1][M][i2][n] = 0. + 0.*I;
                psi[0][i1][i2][n] = 0. + 0.*I;
                psi[M][i1][i2][n] = 0. + 0.*I;
            }
        }
    }
    // phi boundary conditions
    for(nR=0; nR<=TRmax; nR++){
        for(i1=0; i1<=M; i1++){
            for(i2=0; i2<=M; i2++){
                // relaxation gravitational potential BC at "infinity"
                phiR[i1][i2][0][nR] = 0.;
                phiR[i1][i2][M][nR] = 0.;
                phiR[i1][0][i2][nR] = 0.;
                phiR[i1][M][i2][nR] = 0.;
                phiR[0][i1][i2][nR] = 0.;
                phiR[M][i1][i2][nR] = 0.;
            }
        }
    }

///* GRID ALGORITHM AND DATA WRITING *//
    for(n=1, n0=0, n1=1, n2=2; n<=T; n++, n0=(n0+1)%3, n1=(n1+1)%3, n2=(n2+1)%3){
        printf("n=%i out of %i\n", n, T); // prints progress
        ///* SOR GRID ALGORITHM FOR PHI AT EACH TIMESTEP *//
        // set initial conditions for phiR
        if(n==1){   // if this is the first timestep, set initial condition phiR = 0 everywhere
            for(i3=1; i3<=M-1; i3++){
                for(i2=1; i2<=M-1; i2++){
                    for(i1=1; i1<M-1; i1++){
                        phiR[i1][i2][i3][0] = 0.;
                    }
                }
            }
        }
        else{   // if this not the first timestep, use previous iteration of phi as the initial condition for phiR
            for(i3=1; i3<=M-1; i3++){
                for(i2=1; i2<=M-1; i2++){
                    for(i1=1; i1<M-1; i1++){
                        phiR[i1][i2][i3][0] = phi[i1][i2][i3][n0];
                    }
                }
            }
        }

        // update phiR using SOR algorithm
        acc_param = 1.0;
        for(nR=0; nR<TRmax && acc_param >= acc_goal; nR++){
            //printf("n=%i\n",n);
            for(i3=1; i3<=M-1; i3++){
                for(i2=1; i2<=M-1; i2++){
                    for(i1=1; i1<=M-1; i1++){
                        phiR[i1][i2][i3][nR+1] = (1.-w)*phiR[i1][i2][i3][nR] + w/(2.*(pow(1-pow(x(i1),2),2) + pow(1-pow(y(i2),2),2) + pow(1-pow(z(i3),2),2)))
                        * ((1-pow(x(i1),2))*((1-pow(x(i1),2)-DX*x(i1))*phiR[i1+1][i2][i3][nR] + (1-pow(x(i1),2)+DX*x(i1))*phiR[i1-1][i2][i3][nR+1])
                         + (1-pow(y(i2),2))*((1-pow(y(i2),2)-DX*y(i2))*phiR[i1][i2+1][i3][nR] + (1-pow(y(i2),2)+DX*y(i2))*phiR[i1][i2-1][i3][nR+1])
                         + (1-pow(z(i3),2))*((1-pow(z(i3),2)-DX*z(i3))*phiR[i1][i2][i3+1][nR] + (1-pow(z(i3),2)+DX*z(i3))*phiR[i1][i2][i3-1][nR+1])
                         - (pow(R*DX,2)*4*PI*G*m)*pow(cabs(psi[i1][i2][i3][n1]),2));
                    }
                }
            }
            acc_param = cabs((phiR[M/2][M/2][M/2][nR+1] - phiR[M/2][M/2][M/2][nR])/phiR[M/2][M/2][M/2][nR+1]);
            //printf("acc_param = %e\n", acc_param);
            //printf("nR=%i\n", nR);
        }
        //printf("nR=%i\n", nR);

        // set phi = phiR at the last relaxation timestep
        if(n==1){
            for(i3=0; i3<=M; i3++){
                for(i2=0; i2<=M; i2++){
                    for(i1=0; i1<M; i1++){
                        phi[i1][i2][i3][0] = phiR[i1][i2][i3][nR];
                        phi[i1][i2][i3][1] = phiR[i1][i2][i3][nR];
                    }
                }
            }
        }
        else{
            for(i3=0; i3<=M; i3++){
                for(i2=0; i2<=M; i2++){
                    for(i1=0; i1<M; i1++){
                        phi[i1][i2][i3][n1] = phiR[i1][i2][i3][nR];
                        //printf("phi = %e\t", phi[i1][i2][i3][n1]);
                    }
                }
            }
        }
        //printf("nR=%i\n", nR);
        //printf("phi=%e\n", phi[M/2][M/2][M/2][n1]);
        //printf("phiR=%e\n", phiR[M/2][M/2][M/2][nR]);

        ///* GRID ALGORITHM FOR PSI AT EACH TIMESTEP *//
        for(i3=1; i3<=M-1; i3++){
            for(i2=1; i2<=M-1; i2++){          // grid step
                for(i1=1; i1<=M-1; i1++){
                    //printf("n0=%i, n1=%i, n2=%i\n", n0, n1, n2);
                    psi[i1][i2][i3][n2] = psi[i1][i2][i3][n0] + (I*DT*pow(DX*R,-2)/m)
                        *((1-pow(x(i1),2))*((1-pow(x(i1),2)-DX*x(i1))*psi[i1+1][i2][i3][n1] + (1-pow(x(i1),2)+DX*x(i1))*psi[i1-1][i2][i3][n1] - 2*(1-pow(x(i1),2))*psi[i1][i2][i3][n1])
                        + (1-pow(y(i2),2))*((1-pow(y(i2),2)-DX*y(i2))*psi[i1][i2+1][i3][n1] + (1-pow(y(i2),2)+DX*y(i2))*psi[i1][i2-1][i3][n1] - 2*(1-pow(y(i2),2))*psi[i1][i2][i3][n1])
                        + (1-pow(z(i3),2))*((1-pow(z(i3),2)-DX*z(i3))*psi[i1][i2][i3+1][n1] + (1-pow(z(i3),2)+DX*z(i3))*psi[i1][i2][i3-1][n1] - 2*(1-pow(z(i3),2))*psi[i1][i2][i3][n1]))
                        + (I*a*DT*pow(m,-2)/4.)*pow(cabs(psi[i1][i2][i3][n1]),2)*psi[i1][i2][i3][n1]
                        - 2*I*m*DT*phi[i1][i2][i3][n1]*psi[i1][i2][i3][n1];
                    //printf("psi=%e+%e*I\n",creal(psi[i1][i2][i3][n2]), cimag(psi[i1][i2][i3][n2]));
                }
            }
        }

///* INSTABILITY DETECTION */
        if(isnan(psi[M/2][M/2][M/2][n2]) == 1){
            printf("Floating point overflow in psi, possible numerical instability\n");
            //printf("psi=%e+%e*I\n",creal(psi[M/2][M/2][M/2][n2]), cimag(psi[M/2][M/2][M/2][n2]));
            fclose(fpOutPsi);
            fclose(fpOutPhi);
	    getchar();
            exit(2);
        }
        if(isnan(phi[M/2][M/2][M/2][n1]) == 1){
            printf("Floating point overflow in phi, possible numerical instability\n");
            fclose(fpOutPsi);
            fclose(fpOutPhi);
            exit(2);
        }

        ///* WRITE PSI DATA TO FILE IN PYTHON ARRAY FORMAT */
        if((n % data_skip) == 0){  // triggers when n-1 is an integer mulitple of data_skip
            fprintf(fpOutPsi, "[");
            for(i3=0; i3<=M; i3++){
                fprintf(fpOutPsi, "[");
                for(i2=0; i2<=M; i2++){
                    fprintf(fpOutPsi, "[");
                    for(i1=0; i1<=M-1; i1++){
                        //fprintf(fpOutPsi, "%e+%e*I,", creal(psi[i1][i2][i3][n]), cimag(psi[i1][i2][i3][n]));
                        fprintf(fpOutPsi, "%.2e", creal(psi[i1][i2][i3][n1]));
                        if(cimag(psi[i1][i2][i3][n1]) >= 0.)
                            fprintf(fpOutPsi, "+");
                        fprintf(fpOutPsi, "%.2e", cimag(psi[i1][i2][i3][n1]));
                        fprintf(fpOutPsi, "j,");
                    }
                    //fprintf(fpOutPsi, "%e+%e*I", creal(psi[M-1][i2][i3][n]), cimag(psi[i1][i2][i3][n]));
                    fprintf(fpOutPsi, "%.2e", creal(psi[M][i2][i3][n1]));
                    if(cimag(psi[M][i2][i3][n1]) >= 0.)
                        fprintf(fpOutPsi, "+");
                    fprintf(fpOutPsi, "%.2e", cimag(psi[M][i2][i3][n1]));
                    fprintf(fpOutPsi, "j");
                    if(i2==M)
                        fprintf(fpOutPsi, "]");
                    else
                        fprintf(fpOutPsi, "],");
                }
                if(i3==M)
                    fprintf(fpOutPsi, "]");
                else
                    fprintf(fpOutPsi, "],");
            }
            fprintf(fpOutPsi, "],");    // has \n instead of ,
        }

        ///* WRITE PHI DATA TO FILE IN MATHEMATICA ARRAY FORMAT */
//        if((n % data_skip) == 0){  // triggers when n-1 is an integer mulitple of data_skip
//            fprintf(fpOutPhi, "{");
//            for(i3=0; i3<=M; i3++){
//                fprintf(fpOutPhi, "{");
//                for(i2=0; i2<=M; i2++){
//                    fprintf(fpOutPhi, "{");
//                    for(i1=0; i1<=M-1; i1++){
//                        //fprintf(fpOutPhi, "%e,", phi[i1][i2][i3][n1]);
//                        fprintM(fpOutPhi, phi[i1][i2][i3][n1]);
//                        fprintf(fpOutPhi, ",");
//                    }
//                    //fprintM(fpOutPhi, phi[M][i2][i3][n1]); // why isn't this working? output corrupts text file somehow
//                    fprintf(fpOutPhi, "%f", 0.0);
//                    if(i2==M)
//                        fprintf(fpOutPhi, "}");
//                    else
//                        fprintf(fpOutPhi, "},");
//                }
//                if(i3==M)
//                    fprintf(fpOutPhi, "}");
//                else
//                    fprintf(fpOutPhi, "},");
//            }
//            fprintf(fpOutPhi, "},"); // has \n instead of ,
//        }
    }

    fseek(fpOutPsi, -1, SEEK_CUR);
    fprintf(fpOutPsi, "]");
//    fseek(fpOutPhi, -1, SEEK_CUR);
//    fprintf(fpOutPhi, "}");

    fclose(fpOutPsi);
//    fclose(fpOutPhi);

    //Beep(523,500);

    return 0;
}

///* CUSTOM FUNCTION DEFINITION *//
// these are the asymptotic coordinates; x is actually \xi, y is \eta, z is \zeta
double x(int i){
    return DX*i - 1.;
}

double y(int j){
    return DX*j - 1.;
}

double z(int k){
    return DX*k - 1.;
}
