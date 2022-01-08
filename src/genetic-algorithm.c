#include "genetic-algorithm.h"
#include "equations.h"
#include "genotype.h"
#include "randombits.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Get a random individual from the population (with replacement).
 *
 * Disregards their fitness or genotype.
 */
static Individual select_random_individual(const Individual *individuals, const unsigned n_individuals) {
    const unsigned index = uniform() * n_individuals;

    return individuals[index];
}

/* Select the best individual from a tournament of `size` random individuals
 * from the population.
 */
static Individual tournament_selection(const Individual *individuals, const unsigned n_individuals, const unsigned char size) {
    Individual best = select_random_individual(individuals, n_individuals);

    for(unsigned char iter = 1; iter < size; iter++) {
        Individual tmp = select_random_individual(individuals, n_individuals);

        if(tmp.fitness < best.fitness) {
            best = tmp;
        }
    }

    return best;
}

/* Returns a random individual from the population.
 */
static Individual select_individual_with_replacement(const Individual *individuals, const unsigned n_individuals) {
    return tournament_selection(individuals, n_individuals, 50);
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

/* Generate random individual with valid fitness.
 */
static Individual get_random_individual() {
    Genotype g = get_random_genotype();

    double fitness = get_genotype_fitness(g);
    while(fitness == DBL_MAX) {
        g = get_random_genotype();
        fitness = get_genotype_fitness(g);
    }

    return (Individual) {
        .genotype = g,
        .fitness = fitness,
    };
}

Individual run_genetic_algorithm(const unsigned n_individuals) {
    Individual *individuals = (Individual *) malloc(sizeof(Individual) * n_individuals);;
    Individual *new_individuals = (Individual *) malloc(sizeof(Individual) * n_individuals);;

    for(unsigned iter = 0; iter < n_individuals; iter++) {
        individuals[iter] = get_random_individual();
    }
    Individual best = individuals[0];
    best.fitness = DBL_MAX;

    const unsigned n_generations = 1000;
    unsigned generation = 0;
    while((generation < n_generations) && (best.fitness > 10)) {
        const clock_t start = clock();
        printf("Generation %u\n", generation);

        for(unsigned iter = 0; iter < n_individuals; iter++) {
            if(individuals[iter].fitness < best.fitness) {
                best = individuals[iter];
            }
        }
        printf("\n");

        printf("Best fitness so far: %lf (%lf)\n", best.fitness, sqrt(best.fitness));
        Phenotype p = genoype_to_phenotype(best.genotype);
        printf("\tphi: %f\n\tlambda: %f\n\tmu: %f\n\tsigma: %f\n\tdelta: %f\n",
                p.phi, p.lambda, p.mu, p.sigma, p.delta);

        for(unsigned iter = 0; iter < (n_individuals - 1) / 2; iter++) {
            Individual p1 = select_individual_with_replacement(individuals, n_individuals);
            Individual p2 = select_individual_with_replacement(individuals, n_individuals);
            Individual c1, c2;

            individual_crossover(p1, p2, &c1, &c2);
            mutate_individual(&c1);
            mutate_individual(&c2);

            new_individuals[2 * iter] = c1;
            new_individuals[(2 * iter) + 1] = c2;
        }

        new_individuals[n_individuals - 2] = best;
        new_individuals[n_individuals - 1] = best;

#pragma omp parallel default (none) shared (new_individuals) firstprivate (n_individuals) num_threads (8)
        {
#pragma omp for
            for(unsigned iter = 0; iter < n_individuals; iter++) {
                new_individuals[iter].fitness = get_genotype_fitness(new_individuals[iter].genotype);
            }
        }

        Individual *tmp = individuals;
        individuals = new_individuals;
        new_individuals = tmp;

        generation++;
        const clock_t end = clock();
        printf("time: %lfs\n", ((double) end - start) / ((double) 8 * CLOCKS_PER_SEC));
    }

    free(individuals);
    free(new_individuals);
    return best;
}
