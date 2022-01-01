#include "equations.h"
#include "genotype.h"
#include "randombits.h"

Phenotype genoype_to_phenotype(const Genotype g) {
    const double phi = g.phi * ((0.35 + 100.0) / ((double) (1UL << phi_length) - 1)) - 100.0;
    const double lambda = g.lambda * (30000.0 / ((double) (1UL << lambda_length) - 1));
    const double mu = g.mu * (20.0 / ((double) (1UL << mu_length) - 1));
    const double sigma = g.sigma * (1000.0 / ((double) (1UL << sigma_length) - 1));
    const double delta = g.delta * (25000.0 / ((double) (1UL << delta_length) - 1));

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
        .phi = random_UL_length(phi_length),
        .lambda = random_U_length(lambda_length),
        .mu = random_U_length(mu_length),
        .sigma = random_U_length(sigma_length),
        .delta = random_US_length(delta_length),
    };
}

void one_point_crossover(const Genotype p1, const Genotype p2, Genotype *const c1, Genotype *const c2) {
    unsigned char d = uniform() * (phi_length - 1) + 1;
    unsigned char di = phi_length - d;
    c1->phi = ((p1.phi >> d) << d) | ((p2.phi << di) >> di);
    c2->phi = ((p2.phi >> d) << d) | ((p1.phi << di) >> di);

    d = uniform() * (lambda_length - 1) + 1;
    di = phi_length - d;
    c1->lambda = ((p1.lambda >> d) << d) | ((p2.lambda << di) >> di);
    c2->lambda = ((p2.lambda >> d) << d) | ((p1.lambda << di) >> di);

    d = uniform() * (mu_length - 1) + 1;
    di = phi_length - d;
    c1->mu = ((p1.mu >> d) << d) | ((p2.mu << di) >> di);
    c2->mu = ((p2.mu >> d) << d) | ((p1.mu << di) >> di);

    d = uniform() * (sigma_length - 1) + 1;
    di = phi_length - d;
    c1->sigma = ((p1.sigma >> d) << d) | ((p2.sigma << di) >> di);
    c2->sigma = ((p2.sigma >> d) << d) | ((p1.sigma << di) >> di);

    d = uniform() * (delta_length - 1) + 1;
    di = phi_length - d;
    c1->delta = ((p1.delta >> d) << d) | ((p2.delta << di) >> di);
    c2->delta = ((p2.delta >> d) << d) | ((p1.delta << di) >> di);
}

void genotype_crossover(const Genotype p1, const Genotype p2, Genotype *const c1, Genotype *const c2) {
    one_point_crossover(p1, p2, c1, c2);
}

void bit_flip_mutation(Genotype *const g) {
    for(int iter = 0; iter < phi_length; iter++) {
        if(uniform() < 1.0f / phi_length) {
            g->phi ^= 1U << iter;
        }
    }

    for(int iter = 0; iter < lambda_length; iter++) {
        if(uniform() < 1.0f / lambda_length) {
            g->lambda ^= 1U << iter;
        }
    }

    for(int iter = 0; iter < mu_length; iter++) {
        if(uniform() < 1.0f / mu_length) {
            g->mu ^= 1U << iter;
        }
    }

    for(int iter = 0; iter < sigma_length; iter++) {
        if(uniform() < 1.0f / sigma_length) {
            g->sigma ^= 1U << iter;
        }
    }

    for(int iter = 0; iter < delta_length; iter++) {
        if(uniform() < 1.0f / delta_length) {
            g->delta ^= 1U << iter;
        }
    }
}

void mutate_genotype(Genotype *const g) {
    bit_flip_mutation(g);
}
