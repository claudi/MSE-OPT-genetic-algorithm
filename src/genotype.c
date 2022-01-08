#include "equations.h"
#include "genotype.h"
#include "randombits.h"
#include <stdint.h>

Phenotype genoype_to_phenotype(const Genotype g) {
    const double phi = g.phi * ((0.35 + 100.0) / ((double) (1UL << PHI_LENGTH) - 1)) - 100.0;
    const double lambda = g.lambda * (30000.0 / ((double) (1UL << LAMBDA_LENGTH) - 1));
    const double mu = g.mu * (20.0 / ((double) (1UL << MU_LENGTH) - 1));
    const double sigma = g.sigma * (1000.0 / ((double) (1UL << SIGMA_LENGTH) - 1));
    const double delta = g.delta * (25000.0 / ((double) (1UL << DELTA_LENGTH) - 1));

    return (Phenotype) {
        .phi = phi,
        .lambda = lambda,
        .mu = mu,
        .sigma = sigma,
        .delta = delta,
    };
}

Genotype get_random_genotype() {
    return (Genotype) {
        .phi = random_U64_length(PHI_LENGTH),
        .lambda = random_U32_length(LAMBDA_LENGTH),
        .mu = random_U32_length(MU_LENGTH),
        .sigma = random_U32_length(SIGMA_LENGTH),
        .delta = random_U16_length(DELTA_LENGTH),
    };
}

/* Generate two children from applying one point crossover on the parameters of
 * two parents
 */
static void one_point_crossover(const Genotype p1, const Genotype p2, Genotype *const c1, Genotype *const c2) {
    unsigned char d = uniform() * (PHI_LENGTH - 1) + 1;
    unsigned char di = PHI_LENGTH - d;
    c1->phi = ((p1.phi >> d) << d) | ((p2.phi << di) >> di);
    c2->phi = ((p2.phi >> d) << d) | ((p1.phi << di) >> di);

    d = uniform() * (LAMBDA_LENGTH - 1) + 1;
    di = PHI_LENGTH - d;
    c1->lambda = ((p1.lambda >> d) << d) | ((p2.lambda << di) >> di);
    c2->lambda = ((p2.lambda >> d) << d) | ((p1.lambda << di) >> di);

    d = uniform() * (MU_LENGTH - 1) + 1;
    di = PHI_LENGTH - d;
    c1->mu = ((p1.mu >> d) << d) | ((p2.mu << di) >> di);
    c2->mu = ((p2.mu >> d) << d) | ((p1.mu << di) >> di);

    d = uniform() * (SIGMA_LENGTH - 1) + 1;
    di = PHI_LENGTH - d;
    c1->sigma = ((p1.sigma >> d) << d) | ((p2.sigma << di) >> di);
    c2->sigma = ((p2.sigma >> d) << d) | ((p1.sigma << di) >> di);

    d = uniform() * (DELTA_LENGTH - 1) + 1;
    di = PHI_LENGTH - d;
    c1->delta = ((p1.delta >> d) << d) | ((p2.delta << di) >> di);
    c2->delta = ((p2.delta >> d) << d) | ((p1.delta << di) >> di);
}

void genotype_crossover(const Genotype p1, const Genotype p2, Genotype *const c1, Genotype *const c2) {
    one_point_crossover(p1, p2, c1, c2);
}

/* Mutate parameters of a Genotype by fliping bits of its members with
 * probability `1 / length`, where `length` is the bitlength of the parameter.
 */
static void bit_flip_mutation(Genotype *const g) {
    // Consider defining a global bitflip probability of 1/(8*sizeof(Genotype))
    const double prob = 0.5;

    for(int iter = 0; iter < PHI_LENGTH; iter++) {
        if(uniform() < ((float) iter + 1) / PHI_LENGTH) {
            g->phi ^= ((uint64_t) 1) << iter;
        }
    }

    for(int iter = 0; iter < LAMBDA_LENGTH; iter++) {
        if(uniform() < ((float) iter + 1) / LAMBDA_LENGTH) {
            g->lambda ^= ((uint32_t) 1) << iter;
        }
    }

    for(int iter = 0; iter < MU_LENGTH; iter++) {
        if(uniform() < ((float) iter + 1) / MU_LENGTH) {
            g->mu ^= ((uint32_t) 1) << iter;
        }
    }

    for(int iter = 0; iter < SIGMA_LENGTH; iter++) {
        if(uniform() < ((float) iter + 1) / SIGMA_LENGTH) {
            g->sigma ^= ((uint32_t) 1) << iter;
        }
    }

    for(int iter = 0; iter < DELTA_LENGTH; iter++) {
        if(uniform() < ((float) iter + 1) / DELTA_LENGTH) {
            g->delta ^= ((uint16_t) 1) << iter;
        }
    }
}

void mutate_genotype(Genotype *const g) {
    bit_flip_mutation(g);
}

double get_genotype_fitness(Genotype const g) {
    const Phenotype p = genoype_to_phenotype(g);

    return get_phenotype_fitness(p);
}
