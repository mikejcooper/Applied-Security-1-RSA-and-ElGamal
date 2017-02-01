#include "modmul.h"
#include <math.h>
#define WINDOW_SIZE 4   // Sliding Window Exponentiation

/* 

Extension c - Non-Binary (i.e., "windowed") Exponentiation

---------------START ------------------
*/

int exp_by_squaring_iterative(mpz_t b, mpz_t e, mpz_t N){
	mpz_t tmp, y;
	mpz_inits(tmp,y, NULL);
	if(mpz_cmp_si(e,0) == -1){
		mpz_set_ui(tmp,1);
		mpz_invert(b,b,tmp);
		mpz_mul_si(tmp, e, 2);
		mpz_sub(e, e, tmp);
	}
	if(mpz_cmp_si(e,0) == 0){
    	mpz_set_ui(b,1);
    	return 1;
	}    	
	mpz_set_ui(y,1);
	while (mpz_cmp_si(e,1) > 0){
		mpz_set_ui(tmp,2);
		mpz_mod(tmp, e, tmp);
		if(mpz_cmp_si(tmp,0) == 0){
			mpz_mul(b,b,b);
			mpz_mod(b,b,N);
			mpz_tdiv_q_ui(e,e,2);
		}
		else {
			mpz_mul(y,y,b);
			mpz_mod(y,y,N);
			mpz_mul(b,b,b);
			mpz_mod(b,b,N);
			mpz_sub_ui(e,e,1);
			mpz_tdiv_q_ui(e,e,2);
		}
	}
	mpz_mul(b,b,y);
	mpz_clears(tmp,y, NULL);
	return 1;
}

// r = b^e mod N
void squaring_exp(mpz_t r, mpz_t b, mpz_t e, mpz_t N){	
	exp_by_squaring_iterative(b,e,N);
	mpz_mod(r,b,N);
}

/*
    Returns the bits in a mpz_t from start till end.
*/
unsigned int get_bits(mpz_t exp, int start, int end) {
    mpz_t tmp;
    mpz_init(tmp);
    
    // Right shift: tmp >> end, tmp = exp / 2^end 
    mpz_tdiv_q_2exp(tmp, exp, end);
    
    // Masking - defines which bits to keep
    unsigned int mask = ((1 << (start - end + 1)) - 1);
    unsigned int bits = mpz_get_ui(tmp) & mask;
    mpz_clear(tmp);

    return bits;
}

/*
    Sliding Window Exponentiation - Algorithm Pseudo-code from wikipedia 
    r = b^e mod N
    Window Size: k
*/
void window_exp(mpz_t r, mpz_t b, mpz_t e, mpz_t rho_sq, mp_limb_t omega, mpz_t N) {
    int m = 1 << WINDOW_SIZE; // 2^k
    mpz_t table[m/2], b_sq;
    
    /* -------------- PREPROCESSING: The m-ary Method --------------
        Compute: b^2 = b^1 * b^1, | b^4 = b^2 * b^2, b^8 = b^4 * B^2, b^16 = b^8 * b^2  .... and so on |
                                  | b^4.i = b^2.i * b^2                                                |
                                    Store values in lookup table. 
    */
    mpz_init(table[0]);
    mpz_set(table[0],b);
    mpz_init(b_sq);
    // Compute b_sq = b.b mod N
    mont_multiplication(b_sq, b, b, omega, N);
    // Precompute lookup table
    for(int i = 1; i < m/2; i++) {
        mpz_init(table[i]);
        mont_multiplication(table[i], table[i-1], b_sq, omega, N);
    }
    // -------------- PREPROCESSING FINSHED --------------
    
    mpz_set_ui(r,1);
    mont_number(r, r, rho_sq, omega, N);
    
    // #bits in exponent
    int bits_number = mpz_size(e) * mp_bits_per_limb;
    int i = bits_number - 1;

    while(i > -1) {
        // Traverse through leading zeros
        if (get_bits(e, i, i) == 0) {
            i--;
            mont_multiplication(r, r, r, omega, N);
        } 
        else {
            int s = max(i - WINDOW_SIZE + 1, 0);
            while(get_bits(e, s, s) == 0)
                s++;

            for(int h = 0; h < i - s + 1; h++) {
                mont_multiplication(r, r, r, omega, N);
            }
            int u = get_bits(e, i, s);
            mont_multiplication(r, r, table[(u - 1) / 2], omega, N);
            i = s - 1;  
        } 
    }

    for(int i = 0; i < m/2; i++)
        mpz_clear(table[i]);
    
    mpz_clear(b_sq);
}



/* 

Extension c - Non-Binary (i.e., "windowed") Exponentiation

---------------END ------------------
*/























