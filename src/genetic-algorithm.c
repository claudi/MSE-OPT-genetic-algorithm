#include "genetic-algorithm.h"
#include "genotype.h"
#include "randombits.h"

/* Get a random individual from the population (with replacement).
 *
 * Disregards their fitness or genotype.
 */
static Individual get_random_individual(const Individual *individuals, const unsigned n_individuals) {
    const unsigned index = uniform() * n_individuals;

    return individuals[index];
}

/* Select the best individual from a tournament of `size` random individuals
 * from the population.
 */
static Individual tournament_selection(const Individual *individuals, const unsigned n_individuals, const unsigned char size) {
    Individual best = get_random_individual(individuals, n_individuals);

    for(unsigned char iter = 1; iter < size; iter++) {
        Individual tmp = get_random_individual(individuals, n_individuals);

        if(tmp.fitness > best.fitness) {
            best = tmp;
        }
    }

    return best;
}

/* Returns a random individual from the population.
 */
static Individual select_individual_with_replacement(const Individual *individuals, const unsigned n_individuals) {
    return get_random_individual(individuals, n_individuals);
}

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
