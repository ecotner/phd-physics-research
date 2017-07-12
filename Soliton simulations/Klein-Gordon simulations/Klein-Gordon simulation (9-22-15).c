/*

Author: Eric Cotner, UCLA Physics 2015

This program solves the time-dependent Klein-Gordon equation to determine the stability of certain non-topological solitons.

In this specific simulation, we have two fields, A and B (which have been nondimensionalized). A is a real Higgs field and B is a complex
scalar DM field. Since the B field has a U(1) symmetry, there is a conserved charge associated with it, which prevents solitons from
dissolving into free particles (at least on a classical level).

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "./custom.c"

///* GLOBAL VARIABLE DECLARATION */
double k;   // dimensionless scalar mass in vacuum Higgs vev
double nu;  // dimensionless frequency of the complex scalar field
double Q;   // the charge or particle number of the complex scalar field
double R, DT, DX;
int M;
int T=10000; // max number of timesteps

int main(){

///* LOCAL VARIABLE DECLARATION */
    int i1=0, i2=0, i3=0, n=0, n0=0, n1=1, n2=2;

///* FILE I/0 */
    FILE *fpInitA;
    FILE *fpInitB;
    FILE *fpOutA;
    FILE *fpOutB;

    fpInitA = fopen("./Data files/Initial states/StabilityTestInitialA.txt","r");
    fpInitB = fopen("./Data files/Initial states/StabilityTestInitialB.txt","r");
    fpOutA = fopen("./Data files/Output data/Stability tests/(10-25-15.1)StabilityTestOutputA.txt","w+");
    fpOutB = fopen("./Data files/Output data/Stability tests/(10-25-15.1)StabilityTestOutputB.txt","w+");

    int minimize_data = 1;  // if 1, truncates the number of timesteps to write to file (to reduce file size)
    int timesteps_to_keep = 600;    // target number of timesteps to write to file
    int num_frames = 0; // the actual number of timesteps kept, iterated each time you write to file

    // sets up data minimization
    int data_skip;
    if(minimize_data)
        data_skip = T/timesteps_to_keep;  // number of timesteps to skip
    else
        data_skip = 1;

///* SET PARAMETERS FROM FILE */
    fscanf(fpInitA, "R=%lf, DT=%lf, M=%i, ", &R, &DT, &M);
    fscanf(fpInitA, "k=%lf, nu=%lf, Q=%lf\n\n", &k, &nu, &Q);
    fscanf(fpInitB, "R=%lf, DT=%lf, M=%i, ");
    fscanf(fpInitB, "k=%lf, nu=%lf, Q=%lf\n\n");
    printf("R=%e, DT=%e, M=%i, k=%e, nu=%e, Q=%e\n", R, DT, M, k, nu, Q);

    fprintf(fpOutA, "{R->");
    fprintM(fpOutA, R);
    fprintf(fpOutA, ", DT->");
    fprintM(fpOutA, DT);
    fprintf(fpOutA, ", k->");
    fprintM(fpOutA, k);
    fprintf(fpOutA, ", nu->");
    fprintM(fpOutA, nu);
    fprintf(fpOutA, ", Q->");
    fprintM(fpOutA, Q);
    fprintf(fpOutA, ", M->%i, T->%i}\n\n", M, T);

    fprintf(fpOutB, "{R->");
    fprintM(fpOutB, R);
    fprintf(fpOutB, ", DT->");
    fprintM(fpOutB, DT);
    fprintf(fpOutB, ", k->");
    fprintM(fpOutB, k);
    fprintf(fpOutB, ", nu->");
    fprintM(fpOutB, nu);
    fprintf(fpOutB, ", Q->");
    fprintM(fpOutB, Q);
    fprintf(fpOutB, ", M->%i, T->%i}\n\n", M, T);

    DX = 2*R/M;

///* DYNAMIC MEMORY ALLOCATION OF A AND B ARRAYS */
    double ****A = NULL;
    double complex ****B = NULL;
    A = malloc(M*sizeof(double));
    B = malloc(M*sizeof(complex double));
    for(i1=0;i1<=M;i1++){
        A[i1] = malloc(M*sizeof(double));
        B[i1] = malloc(M*sizeof(complex double));
        for(i2=0; i2<=M; i2++){
            A[i1][i2] = malloc(M*sizeof(double));
            B[i1][i2] = malloc(M*sizeof(complex double));
            for(i3=0; i3<=M; i3++){
                A[i1][i2][i3] = malloc(3*sizeof(double));
                B[i1][i2][i3] = malloc(3*sizeof(complex double));
            }
        }
    }

///* SET UP INITIAL CONDITIONS */
    ///* Initial data for A */
    double re, im;

    fseek(fpInitA, 1, SEEK_CUR);
    for(i3=0; i3<=M; i3++){
        fseek(fpInitA, 1, SEEK_CUR);
        for(i2=0; i2<=M; i2++){
            fseek(fpInitA, 1, SEEK_CUR);
            for(i1=0; i1<=M-1; i1++){
                fscanf(fpInitA, "{%lf, %lf}, ", &re, &im);
                A[i1][i2][i3][0] = re;
                A[i1][i2][i3][1] = re;
//                if(A[i1][i2][i3][0] < 0.9){
//                    printf("A[%i][%i][%i] = %e\n", i1, i2, i3, A[i1][i2][i3][0]);
//                }
            }
            fscanf(fpInitA, "{%lf, %lf}}", &re, &im);
            A[M][i2][i3][0] = re;
            A[M][i2][i3][1] = re;
            fseek(fpInitA, 2, SEEK_CUR);
        }
        fseek(fpInitA, 1, SEEK_CUR);
    }
    fseek(fpInitA, 1, SEEK_CUR);
    fclose(fpInitA);

    ///* Initial data for B */
    fseek(fpInitB, 1, SEEK_CUR);
    for(i3=0; i3<=M; i3++){
        fseek(fpInitB, 1, SEEK_CUR);
        for(i2=0; i2<=M; i2++){
            fseek(fpInitB, 1, SEEK_CUR);
            for(i1=0; i1<=M-1; i1++){
                fscanf(fpInitB, "{%lf, %lf}, ", &re, &im);
                B[i1][i2][i3][0] = re + im * I;
                B[i1][i2][i3][1] = (re + im * I) * (cos(nu*DT) - I*sin(nu*DT)); // if B is to be in an energy eigenstate, it must be oscillating with frequency nu
            }
            fscanf(fpInitB, "{%lf, %lf}}", &re, &im);
            B[M][i2][i3][0] = re + im * I;
            B[M][i2][i3][1] = (re + im * I) * (cos(nu*DT) - I*sin(nu*DT));
            fseek(fpInitB, 2, SEEK_CUR);
        }
        fseek(fpInitB, 1, SEEK_CUR);
    }
    fseek(fpInitB, 1, SEEK_CUR);
    fclose(fpInitB);

///* INSTABILITY WARNING */
    char inst_flag;
    double inst_paramA1 = pow(DT*k*B[M/2][M/2][M/2][0],2) + .5*pow(DT,2) - 2;
    double inst_paramA2 = 12*pow(DT/(R*DX),2) + pow(DT*k*B[M/2][M/2][M/2][0],2) + .5*pow(DT,2) - 2;
    double inst_paramB1 = pow(DT*k,2) - 2;
    double inst_paramB2 = 12*pow(DT/(R*DX),2) + pow(DT*k,2) - 2;
    if((cabs(inst_paramA1)>2.0) || (cabs(inst_paramA2)>2.0) || (cabs(inst_paramB1)>2.0) || (cabs(inst_paramB2)>2.0)){
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

///* SET UP BOUNDARY CONDITIONS */
    ///* A boundary conditions */
    for(n=0; n<=2; n++){
        for(i1=0; i1<=M; i1++){
            for(i2=0; i2<=M; i2++){
                A[i1][i2][0][n] = 1.;
                A[i1][i2][M][n] = 1.;
                A[i1][0][i2][n] = 1.;
                A[i1][M][i2][n] = 1.;
                A[0][i1][i2][n] = 1.;
                A[M][i1][i2][n] = 1.;
            }
        }
    }

    ///* B boundary conditions */
    for(n=0; n<=2; n++){
        for(i1=0; i1<=M; i1++){
            for(i2=0; i2<=M; i2++){
                B[i1][i2][0][n] = 0. + 0.*I;
                B[i1][i2][M][n] = 0. + 0.*I;
                B[i1][0][i2][n] = 0. + 0.*I;
                B[i1][M][i2][n] = 0. + 0.*I;
                B[0][i1][i2][n] = 0. + 0.*I;
                B[M][i1][i2][n] = 0. + 0.*I;
            }
        }
    }

///* RUN GRID ALGORITHM */
    for(n=1, n0=0, n1=1, n2=2; n<=T; n++, n0=(n0+1)%3, n1=(n1+1)%3, n2=(n2+1)%3){
        printf("n=%i out of %i\n", n, T); // prints progress
        for(i1=1; i1<=M-1; i1++){
            for(i2=1; i2<=M-1; i2++){
                for(i3=1; i3<=M-1; i3++){
                    ///* Grid algorithm for A (second-order CTCS scheme) */
                    A[i1][i2][i3][n2] = (2 - 6*pow(DT/DX,2) - pow(DT*k*cabs(B[i1][i2][i3][n1]),2) + .5*pow(DT,2))*A[i1][i2][i3][n1] - A[i1][i2][i3][n0]
                        -.5*pow(DT,2)*pow(A[i1][i2][i3][n1],3) + pow(DT/DX,2)*(A[i1+1][i2][i3][n1] + A[i1-1][i2][i3][n1] + A[i1][i2+1][i3][n1]
                        + A[i1][i2-1][i3][n1] + A[i1][i2][i3+1][n1] + A[i1][i2][i3-1][n1]);
                    ///* Grid algorithm for B (second-order CTCS scheme) */
                    B[i1][i2][i3][n2] = (2 - 6*pow(DT/DX,2) - pow(DT*k*A[i1][i2][i3][n1],2))*B[i1][i2][i3][n1] - B[i1][i2][i3][n0]
                        + pow(DT/DX,2)*(B[i1+1][i2][i3][n1] + B[i1-1][i2][i3][n1] + B[i1][i2+1][i3][n1] + B[i1][i2-1][i3][n1] + B[i1][i2][i3+1][n1]
                        + B[i1][i2][i3-1][n1]);
                }
            }
        }

///* INSTABILITY DETECTION */
        if(isnan(A[M/2][M/2][M/2][n2]) == 1){
            printf("Floating point overflow in A, possible numerical instability\n");
            //printf("psi=%e+%e*I\n",creal(psi[M/2][M/2][M/2][n2]), cimag(psi[M/2][M/2][M/2][n2]));
            fclose(fpOutA);
            fclose(fpOutB);
            exit(2);
        }
        if(isnan(B[M/2][M/2][M/2][n1]) == 1){
            printf("Floating point overflow in B, possible numerical instability\n");
            fclose(fpOutA);
            fclose(fpOutB);
            exit(2);
        }

///* WRITE DATA TO FILE (writes in a Mathematica-readable format that can be read into memory timeslice-by-timeslice */
        ///* Writing B data to file (B is a complex scalar field) */
        if((n % data_skip) == 0){  // triggers when n-1 is an integer multiple of data_skip
            fprintf(fpOutB, "{");
            for(i3=0; i3<=M; i3++){
                fprintf(fpOutB, "{");
                for(i2=0; i2<=M; i2++){
                    fprintf(fpOutB, "{");
                    for(i1=0; i1<=M-1; i1++){
                        fprintM(fpOutB, creal(B[i1][i2][i3][n1]));
                        fprintf(fpOutB, "+");
                        fprintM(fpOutB, cimag(B[i1][i2][i3][n1]));
                        fprintf(fpOutB, "*I,");
                    }
                    fprintM(fpOutB, creal(B[M][i2][i3][n1]));
                    fprintf(fpOutB, "+");
                    fprintM(fpOutB, cimag(B[M][i2][i3][n1]));
                    fprintf(fpOutB, "*I");
                    if(i2==M)
                        fprintf(fpOutB, "}");
                    else
                        fprintf(fpOutB, "},");
                }
                if(i3==M)
                    fprintf(fpOutB, "}");
                else
                    fprintf(fpOutB, "},");
            }
            fprintf(fpOutB, "}\n");    // has \n instead of ,
        }
        ///* Write A data to file (A is a hermitian scalar field) */
        if((n % data_skip) == 0){  // triggers when n-1 is an integer multiple of data_skip
            fprintf(fpOutA, "{");
            for(i3=0; i3<=M; i3++){
                fprintf(fpOutA, "{");
                for(i2=0; i2<=M; i2++){
                    fprintf(fpOutA, "{");
                    for(i1=0; i1<=M-1; i1++){
                        fprintM(fpOutA, A[i1][i2][i3][n1]);
                        fprintf(fpOutA, ",");
                    }
                    fprintM(fpOutA, A[M][i2][i3][n1]);
                    if(i2==M)
                        fprintf(fpOutA, "}");
                    else
                        fprintf(fpOutA, "},");
                }
                if(i3==M)
                    fprintf(fpOutA, "}");
                else
                    fprintf(fpOutA, "},");
            }
            fprintf(fpOutA, "}\n");    // has \n instead of ,
        }

    } // end of n for loop

    fclose(fpOutA);
    fclose(fpOutB);

    return 0;
}
