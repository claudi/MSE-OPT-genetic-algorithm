#include "genetic-algorithm.h"

Individual select_individual_with_replacement(const Individual *individuals, const unsigned n_individuals);
Individual individual_crossover(const Individual g1, const Individual g2);
Individual mutate_individual(const Individual g);
