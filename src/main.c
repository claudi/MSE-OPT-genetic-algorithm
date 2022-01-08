#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "equations.h"
#include "genetic-algorithm.h"
#include "genotype.h"
#include "randombits.h"

int main(void) {
    randomize();
    Individual best = run_genetic_algorithm(5000);
    Phenotype p = genoype_to_phenotype(best.genotype);

    double x[12] = { 15329.0 };
    const double y[12] = { 15329.0, 14177.0, 13031.0, 9762.0, 11271.0, 8688.0, 7571.0, 6983.0, 4778.0, 2067.0, 1586.0, 793.0 };

    model_prediction(x[0], x, 12, &p);

    for(unsigned iter = 0; iter < 12; iter++) {
        printf("%d\t%lf\t%lf\n", iter, y[iter], x[iter]);
    }

    return 0;
}
