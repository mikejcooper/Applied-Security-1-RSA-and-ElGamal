#include "modmul.h"
#include <math.h>


/* 

Extension b - Chinese Remainder Theorem

---------------START ------------------
*/

/* 
CRT Theory:

If 'N = a mod pq' , then 'N = a mod p' and 'N = a mod q'

Therefore find:
	1. x == c^d mod p 
	2. x == c^d mod q
					Then solve the CR problem. 
This works becuase, 
	'c^(p-1) == 1 mod p' (likewise for q)
	Therefore at worst we have to compute c^(p-1) || c^(q-1) istead of c^d

*/
void crt(mpz_t N, mpz_t d, mpz_t p, mpz_t q, mpz_t dP, mpz_t dQ, mpz_t iP, mpz_t iQ, mpz_t c, mpz_t m){

	mpz_t a , b, s, f, h, i;
	mpz_inits(a, b, s, f, h, i, NULL);
	/* 
	Precomputed:
		dP = d mod p-1
		dQ = d mmod q-1
	Compute: 
		a == c^dP mod p
		b == c^dQ mod q
	*/
	mpz_powm_sec(a, c, dP, p);
	mpz_powm_sec(b, c, dQ, q);
	/*
	Precomputed: 
		qi = q^-1 mod p
	Compute:
		h == qi.(a - b) mod p
		where
			s = a - b
			f = iq.c
	*/
	mpz_sub(s, a, b);
	mpz_mul(f, iQ, s);
	mpz_mod(h, f, p);
	/*
	Compute:
		m = b + h.q
	*/
	mpz_mul(i, h, q);
	mpz_add(m, b, i);

	mpz_clears(a, b, s, f, h, i, NULL);
}
/* 

Extension b - Chinese Remainder Theorem

---------------END ------------------
*/


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
void modular_exponentiation(mpz_t r, mpz_t b, mpz_t e, mpz_t N){	
	exp_by_squaring_iterative(b,e,N);
	mpz_mod(r,b,N);
}


/* 

Extension c - Non-Binary (i.e., "windowed") Exponentiation

---------------END ------------------
*/

/* 

Extension d - Montgomery multiplication

---------------Start ------------------
*/




/* 

Extension d - Montgomery multiplication

---------------END ------------------
*/






















