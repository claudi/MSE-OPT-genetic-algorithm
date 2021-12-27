#include "equations.h"
#include "genotype.h"

Phenotype genoype_to_phenotype(const Genotype g) {
    const double phi = g.phi * ((0.35 + 100.0) / ((double) (1UL << 34) - 1)) - 100.0;
    const double lambda = g.lambda * (30000.0 / ((double) (1UL << 25) - 1));
    const double mu = g.mu * (20.0 / ((double) (1UL << 25) - 1));;
    const double sigma = g.sigma * (1000.0 / ((double) (1UL << 17) - 1));
    const double delta = g.delta * (25000.0 / ((double) (1UL << 15) - 1));

    return (Phenotype) {
        .phi = phi,
        .lambda = lambda,
        .mu = mu,
        .sigma = sigma,
        .delta = delta,
    };
}
