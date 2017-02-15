#include "modmul.h"
#include <math.h>
#include <stdint.h>

#define r_bit_length 160 // SHA-1 security standard
#define random_values 100


/* 

Extension a - Pseudo-random number generator

---------------Start ------------------
*/

mpz_t prng_values[random_values];
int prng_index = 0;

// get selection of r_bit_length random numbers from /dev/urandom
void fill_prng_array() {
	FILE *fp;
	uint8_t input[r_bit_length/8];
	for (int i = 0; i < random_values; i++){
		fp = fopen("/dev/urandom", "r"); 	// dev/urandom returns random numbers based on /dev/random
		fread(&input, 1, sizeof(input), fp);
		mpz_init(prng_values[i]);
		mpz_import(prng_values[i], sizeof(input), 1, sizeof(input[0]), 0, 0, input);
	}
	fclose(fp);
}

void clear_prng_array(){
	for (int i = 0; i < random_values; i++){
		mpz_clear(prng_values[i]);
	}
}

int get_prng(){
	prng_index++;
	return prng_index;
}

// Pseudo-random number generator
void prng(mpz_t rand_Num, mpz_t q) {
    mpz_t seed;
    gmp_randstate_t r_state;
	mpz_init(seed);
    //  Create random state - Mersenne Twister algorithm. 
    gmp_randinit_default (r_state);
    // Create seed using randomness from /dev/urandom
	mpz_swap(seed, prng_values[get_prng()]);
    // Random state seeding 
    gmp_randseed(r_state, seed);
	// Generate a uniform random integer in the range 0 to q-1, inclusive.
	// Required parameter for key in ElGamal's Enc Algo.         
    // Generate r using the random state
    // Compute r mod q
 	mpz_urandomm(rand_Num,r_state,q);
	     
    gmp_randclear(r_state);
    mpz_clear(seed);
}

/* 

Extension a - Pseudo-random number generator

---------------END ------------------
*/
