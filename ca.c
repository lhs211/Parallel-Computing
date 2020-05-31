/*
Candidate Number: 086696
*/

#include<stdio.h>
#include<math.h>
#include<omp.h>

int main()
{
    // timing variables
    double tstart, tend;

    // start timer if compiled with OpenMP
    #ifdef _OPENMP
        tstart=omp_get_wtime();
    #endif

    // initialize constants
    const double PRECISION = 0.001;
    const int SIZE = 1000;
    const int TIMESTEPS = 1500;
    const double TIME_INTERVAL = 0.05;
    const double GAUSSIAN = 0.03;
    const double VX = 0.01;
    const double VY = 0.01;
    const double X0 = 0.1;
    const double Y0 = 0.1;
    const char filename[] = "u.dat";

    // make 2D array to hold u values 1000x1000
    double data[1000][1000];
    double data2[1000][1000];


    double x_component, y_component;
    double numerator, denominator;
    int i, j;

    // calculate starting values
    #pragma omp parallel for collapse(2) default(none) shared(data) private(i,j,x_component,y_component, numerator, denominator)
    for(i=0; i<SIZE; i++){
        for(j=0; j<SIZE; j++){
            // compute seperate parts of equation
            x_component = pow(((i*PRECISION) - X0),2);
            y_component = pow(((j*PRECISION) - Y0),2);
            numerator = x_component + y_component;
            denominator = 2.0 * pow(GAUSSIAN,2);
            // update array
            data[i][j] = exp((numerator / denominator) * -1.0);
        }
    }

    double x_change;
    double y_change;
    int time;
    // begin timsteps
    for(time=0; time<TIMESTEPS; time++){

        // make a copy of the up to date array
        #pragma omp parallel for collapse(2) default(none) shared(data, data2) private(i,j)
        for(i=0; i<SIZE; i++){
            for(j=0; j<SIZE; j++){
                data2[i][j] = data[i][j];
            }
        }

        // update values in other array
        #pragma omp parallel for collapse(2) default(none) shared(data, data2) private(i,j,x_change,y_change)
        for(i=1; i<SIZE; i++){
            for(j=1; j<SIZE; j++){
                x_change = -1.0 *  VX * ((data2[i][j] - data2[i-1][j]) / PRECISION);
                y_change = -1.0 * VY * ((data2[i][j] - data2[i][j-1]) / PRECISION);
                data[i][j] = data2[i][j] + ((x_change + y_change) * TIME_INTERVAL);
            }
        }

    }

    // write values to file
    FILE *fp = fopen(filename, "w");
    for (i = 0; i < SIZE; i++){
        for(j = 0; j < SIZE; j++){
            fprintf(fp, "%f %f %f\n", i*PRECISION,j*PRECISION,data[i][j]);
        }
    }
    fclose(fp);

    // end timer if compiled with OpenMP
    #ifdef _OPENMP
        tend=omp_get_wtime();
        printf("Time taken to calculate= %f\n", tend-tstart);
    #endif

    return 0;
}