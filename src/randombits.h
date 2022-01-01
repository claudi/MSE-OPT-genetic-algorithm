/* Generate random `flaot` between `0.0` and `1.0`.
 */
float uniform(void);

/* Utility functions to simplify the use of `ran1` and `idum`.
 */
void randomize(void);

/* Generate random unsigned long integer of at most `width` bits
 */
unsigned long random_UL_length(unsigned char width);

/* Generate random unsigned integer of at most `width` bits
 */
unsigned random_U_length(unsigned char width);

/* Generate random unsigned short integer of at most `width` bits
 */
unsigned short random_US_length(unsigned char width);
