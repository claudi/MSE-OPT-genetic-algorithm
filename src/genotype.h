#pragma once
#include <stdint.h>
#include "equations.h"

#define PHI_LENGTH 34
#define LAMBDA_LENGTH 25
#define MU_LENGTH 25
#define SIGMA_LENGTH 17
#define DELTA_LENGTH 15

/* Structure containing the discretisation of the ODE parameters as unsigned
 * integers of appropriate length.
 *
 * This is done for compatibility with Holland's Convergence Theorem, which
 * works under the setting of genes consisting in unsigned integers expressed
 * in binary.
 *
 * For each parameter we need the theoretical range, an effective search range
 * and a reasonable precision, thus fixing the range and discretisation
 * formula for the genotype.
 *
 * The structure is packed as to occupy the least space possible, and in total
 * we achieve a `sizeof(Phenotype) = 40`, while `sizeof(Genotype) = 16`.
 */
typedef struct {
    /* - Theoretical search range: (-∞, α], where α = 0.3489494085776018 is the
     * neat population growth rate estimated from the first epoch.
     *
     * - Effective search range: [-100, 0.35].
     *
     * - Discretisation step: 1e-8.
     *
     * - Integer search range: [0, 2^34 - 1].
     */
    uint64_t phi : PHI_LENGTH;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 3000].
     *
     * - Discretisation step: 1e-4.
     *
     * - Integer search range: [0, 2^25 - 1].
     */
    uint32_t lambda : LAMBDA_LENGTH;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 20].
     *
     * - Discretisation step: 1e-6.
     *
     * - Integer search range: [0, 2^25 - 1].
     */
    uint32_t mu : MU_LENGTH;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 1000].
     *
     * - Discretisation step: 1e-2.
     *
     * - Integer search range: [0, 2^17 - 1].
     */
    uint32_t sigma : SIGMA_LENGTH;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 25000].
     *
     * - Discretisation step: 1.
     *
     * - Integer search range: [0, 2^15 - 1].
     */
    uint16_t delta : DELTA_LENGTH;
} Genotype;

/* Function to convert from discrertised coefficients, used by the genetic
 * algorithm as unsigned integers, to floating point numbers.
 */
Phenotype genoype_to_phenotype(const Genotype g);

/* Generate random genotype.
 */
Genotype get_random_genotype();

/* Mix and match parts of two genotypes to form two children genotypes.
 */
void genotype_crossover(const Genotype p1, const Genotype p2, Genotype *const c1, Genotype *const c2);

/* Randomly mutate bits of a genotype.
 */
void mutate_genotype(Genotype *const g);

/* Calculate fitness of a genotype through the sum of the squared error between
 * the predictions made from the associated phenotype and the observations.
 */
double get_genotype_fitness(Genotype const g);
