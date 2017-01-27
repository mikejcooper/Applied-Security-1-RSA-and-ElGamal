#include "modmul.h"
#include <time.h>
#include <math.h>
#include "sodium.h"

int prng_bit_length = 128;
mpz_t prng_values[100];
int prng_index = 0;

// get selection of X-bit length random numbers from /dev/urandom
void fill_prng_array() {
	FILE *fp;
	uint8_t input[prng_bit_length/8];
	// /dev/urandom returns random numbers based on /dev/random
	for (int i = 0; i < 100; i++){
		fp = fopen("/dev/urandom", "r");
		fread(&input, 1, sizeof(input), fp);
		mpz_init(prng_values[i]);
		mpz_import(prng_values[i], sizeof(input), 1, sizeof(input[0]), 0, 0, input);
		// gmp_printf( "%zu\n", mpz_sizeinbase(prng_values[i], 2));
  //   	gmp_printf( "%Zd\n", prng_values[i]);
	}
	fclose(fp);
}

void clear_prng_array(){
	for (int i = 0; i < 100; i++){
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
	mpz_set(seed, prng_values[get_prng()]);
    // Random state seeding 
    gmp_randseed(r_state, seed);
	// Generate a uniform random integer in the range 0 to q-1, inclusive - 
	// Required parameter for key in ElGamal's Enc Algo. 
 	mpz_urandomm(rand_Num,r_state,q);
	     
    gmp_randclear(r_state);
    mpz_clear(seed);

    // gmp_printf( "%zu\n", mpz_sizeinbase(rand_Num, 2));
    // gmp_printf( "%Zd\n", rand_Num);
}



