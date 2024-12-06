#include "rands.h"
#include <stdlib.h>
#include <math.h>


int uniform_randInt(int M, int N){
    return floorl(((double)M + ((double)rand() / (RAND_MAX + 1.0)) * (double)(N-M+1)));
}


double uniform_randR(int N){
    return ((double)N * ((double)rand() / (RAND_MAX + 1.0)));
}


void gauss_rand(double *v, int d){
    int i;
    double y1, y2;
    for(i=0;i<d;i++){
        y1 = 1.0 - uniform_randR(1); //y1 uniformly randomly distributed in (0,1]
        y2 = 1.0 - uniform_randR(1); //y2 uniformly randomly distributed in (0,1]
        v[i] = sqrt((-2.0)*log(y1)) * cos(2*M_PI*y2);
        if(v[i] < 0.0)
            v[i] = (-1.0) * v[i];
    }
}

void gauss_rand2(double *v, int d){
    int i;
    double y1, y2;
    for(i=0;i<d;i++){
        y1 = 1.0 - uniform_randR(1); //y1 uniformly randomly distributed in (0,1]
        y2 = 1.0 - uniform_randR(1); //y2 uniformly randomly distributed in (0,1]
        v[i] = sqrt((-2.0)*log(y1)) * cos(2*M_PI*y2);
    }
}

double uniformRandDouble(double N){
    return (N * ((double)rand() / (RAND_MAX + 1.0)));
}
