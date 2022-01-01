#include <stdint.h>

/* Generate random `flaot` between `0.0` and `1.0`.
 */
float uniform(void);

/* Utility functions to simplify the use of `ran1` and `idum`.
 */
void randomize(void);

/* Generate random unsigned long integer of at most `width` bits
 */
uint64_t random_U64_length(unsigned char width);

/* Generate random unsigned integer of at most `width` bits
 */
uint32_t random_U32_length(unsigned char width);

/* Generate random unsigned short integer of at most `width` bits
 */
uint16_t random_U16_length(unsigned char width);
