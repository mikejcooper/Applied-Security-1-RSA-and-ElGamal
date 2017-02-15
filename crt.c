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
    mpz_init(a); mpz_init(b); mpz_init(s); mpz_init(f); mpz_init(h); mpz_init(i); 

	/* 
	Precomputed:
		dP = d mod p-1
		dQ = d mmod q-1
	Compute: 
		a == c^dP mod p
		b == c^dQ mod q
	*/
	mpz_powm(a, c, dP, p);
	mpz_powm(b, c, dQ, q);
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

	mpz_clear(a); mpz_clear(b); mpz_clear(s); mpz_clear(f); mpz_clear(h); mpz_clear(i); 

}

/* 

Extension b - Chinese Remainder Theorem

---------------END ------------------
*/






















