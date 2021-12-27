#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "equations.h"
#include "genotype.h"

int main(void) {
    double x0 = 15329.0;
    double x[12];
    for(unsigned iter = 1; iter < 12; iter++) {
        x[iter] = 0;
    }

    Phenotype p = (Phenotype) {
        .phi = 0.241158,
        .lambda = 2521.740035,
        .mu = 5.003401,
        .sigma = 984.184144,
        .delta = 12878.810999,
    };

    model_prediction(x0, x, 12, &p);
    for(unsigned iter = 0; iter < 12; iter++) {
        printf("x[%u] = %lf\n", iter, x[iter]);
    }

    printf("%zu\n", sizeof(Genotype));
    printf("%zu\n", sizeof(Phenotype));

    return 0;
}
