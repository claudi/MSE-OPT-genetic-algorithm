#pragma once
#include "genotype.h"

typedef struct {
    Genotype genotype;
    double fitness;
} Individual;

/* Main function to run the genetic algorithm, based in [1].
 */
Individual run_genetic_algorithm(const unsigned n_individuals);
