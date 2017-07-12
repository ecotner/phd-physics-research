#include <stdio.h>
#include <complex.h>
#include <math.h>


/* GLOBAL VARIABLES */
const int M=100; // number of grid points
const int T=3000; // number of timesteps
const double dx=0.01; // grid size
const double dt=1.0; // timestep size
const double a=1.0;


int main(){
/* LOCAL VARIABLES */
    double foo;
    int i=0; int n=0;

    double psi[M][T+1];

/* FILE I/0 SETUP */
    FILE *fpInit;
    FILE *fpOut;
    fpOut = fopen("test3.txt","w+");
    fpInit = fopen("InitData.txt","r");

/* SET INITIAL DATA */          //some kind of error is happening here for large values of T; but this loop never calls T?
    for(i=0; i<=M-1; i++){
        fscanf(fpInit, "%lf", &foo);
        psi[i][0] = foo;
//        printf("%.17f\n", psi[i][0]);
    }
    fclose(fpInit);

/* BOUNDARY CONDITIONS */
/*    for(n=1; n<=T; n++){
        psi[0][n] = 0.0;
        psi[M-1][n] = 0.0;
    }
*/
/* RUN GRID ALGORITHM */
/*    for(n=0; n<=T-1; n++){        // timestep
        printf("n=%i\n", n);
        for(i=1; i<=M-2; i++){  // grid step
            psi[i][n+1] = psi[i][n] + (dt/2)*(psi[i+1][n] - 2*psi[i][n] + psi[i-1][n]); // - dt*a*pow((i*dx - 50),2)*psi[i][n];
        }
    }
*/
/* WRITE DATA TO FILE */
/*    for(n=0; n<=T; n++){
        for(i=0; i<=M-1; i++){
            fprintf(fpOut, "%f\t", psi[i][n]);
        }
        fprintf(fpOut, "\n");
    }
*/
    fclose(fpOut);

    return 0;
}
