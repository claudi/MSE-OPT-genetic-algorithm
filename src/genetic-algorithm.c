#include "genetic-algorithm.h"
#include "genotype.h"

static Individual select_individual_with_replacement(const Individual *individuals, const unsigned n_individuals);

/* Mix and match  genotypes of two individuals to form two children genotypes.
 */
static void individual_crossover(const Individual p1, const Individual p2, Individual *const c1, Individual *const c2) {
    genotype_crossover(p1.genotype, p2.genotype, &(c1->genotype), &(c2->genotype));
}

/* Randomly mutate bits of a genotype of an individual.
 */
static void mutate_individual(Individual *const individual) {
    mutate_genotype(&(individual->genotype));
}
