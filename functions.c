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
-------EXTRA---------
We want to compute m = c^d mod n. 
If we know the pair (c^d mod p, c^d mod q) 
then the CRT tells us there is a unique value of c^d mod n in the range [0, n-1].

To recover x from its CRT representation (x1, x2) we use Garner's formula
*/

void crt_special(mpz_t x2, mpz_t x1, mpz_t q, mpz_t p, mpz_t m){

	mpz_t iQ, xs, xsiQ, x, h, hq;
	mpz_inits(iQ, xs, xsiQ, x, h, hq, NULL);
	
	/* 
	Compute: 
		x = x2 + hq
		where 
			hq = h.q
			h = xsiQ mod p     (h = ((x1 - x2))((1/q) mod p)) mod p)
			xsiQ = xs*iQ
			iQ = 1/q mod p
			xs = x1 - x2

	*/
	mpz_sub(xs, x1, x2);
	mpz_invert(iQ, q, p);
	mpz_mul(xsiQ, xs, iQ);
	mpz_mod(h, xsiQ, p);	
	mpz_mul(hq, h, q);
	mpz_add(m, x2, hq);

	mpz_clears(iQ, xs, xsiQ, x, h, hq, NULL);
}

void crt_special_case(mpz_t N, mpz_t d, mpz_t p, mpz_t q, mpz_t dP, mpz_t dQ, mpz_t iP, mpz_t iQ, mpz_t c, mpz_t m){
	mpz_t x1, x2;
	mpz_inits(x1, x2, NULL);
	/* 
	Compute: 
		x1 = c^d mod p
		x2 = c^d mod q
	*/
	mpz_powm_sec(x1, c, d, p);
	mpz_powm_sec(x2, c, d, q);
	/* 
	Special case assumes (x1,x2) are known
	We could computer x1, x2 using Euler's Theorem
	*/
	crt_special(x2, x1, q, p, m);

	mpz_clears(x1, x2, NULL);
}

/* 

Extension b - Chinese Remainder Theorem

---------------END ------------------
*/


/* 

Extension c - Non-Binary (i.e., "windowed") Exponentiation

---------------START ------------------
*/


void binaryToHex(const char *inStr, char *outStr) {
    // outStr must be at least strlen(inStr)/4 + 1 bytes.
    static char hex[] = "0123456789ABCDEF";
    int len = strlen(inStr) / 4;
    int i = strlen(inStr) % 4;
    char current = 0;
    if(i) { // handle not multiple of 4
        while(i--) {
            current = (current << 1) + (*inStr - '0');
            inStr++;
        }
        *outStr = hex[(int)current];
        ++outStr;
    }
    while(len--) {
        current = 0;
        for(i = 0; i < 4; ++i) {
            current = (current << 1) + (*inStr - '0');
            inStr++;
        }
        *outStr = hex[(int)current];
        ++outStr;
    }
    *outStr = 0; // null byte
}

void modular_exponentiation_window2(int m, mpz_t r, mpz_t b, mpz_t e, mpz_t N){
	mpz_t tmp, running_total;
	mpz_inits(tmp,running_total, NULL);

	int bit_length = mpz_sizeinbase(e, 2);
	int window_size = (m >= bit_length) ? bit_length : m;
	int bits_per_m = ceil(bit_length/window_size + 0.499);
	printf("%d\n",bits_per_m );
	
	// (MSB) 2^x .... 2^0 (LSB)    (Big Endian)
	char* binary = mpz_get_str(NULL, 2, e);
	// printf("%s\n", binary);

	// for each m groups 
	for (int i = 0; i < bit_length; i+=bits_per_m){
		char group[3000] = "";
		int block_bits = bits_per_m;
		//if odd bit legth, deal with trailing bits
		if (i == bit_length*(window_size-1)){
			block_bits = bits_per_m*window_size - bit_length;
		}
		// for bit in group 
		for(int j = 0; j < block_bits; j++){
			group[j] = binary[i + j];
		}

		// printf("%s\n", group);
		char group2[3000] = "";
		binaryToHex(group, group2);
		// printf("%s\n", group2);
		// char* buffer = NULL;
		// size_t bufferSize = 0;
		// FILE* myStream = open_memstream(&buffer, &bufferSize);

		// CONVERT BINARY TO GMP NUMBER

		mpz_mul(running_total, running_total, tmp);
		// running_total *= Wi;
		mpz_powm_sec(running_total, b, tmp, N);
		// printf("%lld\n", running_total);
	}

	mpz_set(r,tmp);
	mpz_clears(tmp,running_total, NULL);
}


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

// x^n
// 	   if n < 0 then
//       x := 1 / x;
//       n := -n;
//     if n = 0 then return 1
//     y := 1;
//     while n > 1 do
//       if n is even then 
//         x := x * x;
//         n := n / 2;
//       else
//         y := x * y;
//         x := x * x;
//         n := (n – 1) / 2;
//     return x * y


void modular_exponentiation_window(int m, mpz_t r, mpz_t b, mpz_t e, mpz_t N){
	mpz_t tmp, running_total;
	mpz_inits(tmp,running_total, NULL);

	int bit_length = mpz_sizeinbase(e, 2);
	int window_size = (m >= bit_length) ? bit_length : m;
	int bits_per_m = ceil(bit_length/window_size + 0.499);
	printf("%d\n",bits_per_m );
	
	// (MSB) 2^x .... 2^0 (LSB)    (Big Endian)
	char* binary = mpz_get_str(NULL, 2, e);
	// printf("%s\n", binary);

	// for each m groups 
	for (int i = 0; i < bit_length; i+=bits_per_m){
		char group[3000] = "";
		int block_bits = bits_per_m;
		//if odd bit legth, deal with trailing bits
		if (i == bit_length*(window_size-1)){
			block_bits = bits_per_m*window_size - bit_length;
		}
		// for bit in group 
		for(int j = 0; j < block_bits; j++){
			group[j] = binary[i + j];
		}

		// printf("%s\n", group);
		char group2[3000] = "";
		binaryToHex(group, group2);
		// printf("%s\n", group2);
		// char* buffer = NULL;
		// size_t bufferSize = 0;
		// FILE* myStream = open_memstream(&buffer, &bufferSize);

		// CONVERT BINARY TO GMP NUMBER

		mpz_mul(running_total, running_total, tmp);
		// running_total *= Wi;
		mpz_powm_sec(running_total, b, tmp, N);
		// printf("%lld\n", running_total);
	}

	mpz_set(r,tmp);
	mpz_clears(tmp,running_total, NULL);
}


// r = b^e mod N
void modular_exponentiation(mpz_t r, mpz_t b, mpz_t e, mpz_t N){
	// int winodow_size = 2;
	// modular_exponentiation_window(winodow_size,r,b,e,N);
	
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






















