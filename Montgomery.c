#include "modmul.h"
#include <math.h>

/* 

Extension d - Montgomery multiplication

---------------Start ------------------
*/


/*
    Returns the larger of two integers
*/
int max(int x, int y) {
    return x > y ? x : y;
}

/*
    Montgomery multiplication
    r = a.b mod N
*/
void mont_multiplication(mpz_t r, mpz_t a, mpz_t b, mp_limb_t w, mpz_t N) {
    mp_limb_t u, b_i, a_0, r_0;
    mpz_t t; // Using temp t therefore r and x can be passed as same variable (adkin to GMP functions)
    mpz_init(t);
    mpz_set_ui(t, 0);  

    for (mp_size_t i = 0; i < mpz_size(N); i++) {
        a_0 = mpz_getlimbn(a, 0); // 0-th limb a
        r_0 = mpz_getlimbn(t, 0); // 0-th limb t
        b_i = mpz_getlimbn(b, i); // i-th limb b
        u = (r_0 + b_i * a_0) * w; // u = (r_0 + b_i.a_0).w (mod b)

        mpz_addmul_ui(t, a, b_i); // t = t + b_i.a
        mpz_addmul_ui(t, N, u);   // t = t + u.N
        mpz_tdiv_q_2exp(t, t, mp_bits_per_limb);  // t = t / b = t + b_i.a + u.N) / b
    }
    
    // if (t > N) then t = t - N
    if(mpz_cmp(t,N) >= 0)
        mpz_sub(t,t,N);
    
    mpz_swap(r, t);
    mpz_clear(t);
}

/*
    Computing Omega w: 
        ω = −N^−1 (mod ρ)
        ρ = b^k, for the smallest kst. b^k > N
*/
void mont_omega(mp_limb_t* w, mpz_t N) {
    *w = 1;  // w = 1 mod b
    mp_limb_t b = mpz_getlimbn(N, 0);  // 0-th limb N
    
    for (mp_size_t i = 1; i <= mp_bits_per_limb; i++)
        *w *= *w * b; // w = w * w * N mod ρ

    // w = -w mod ρ
    *w = -*w;
}

/*
    Computing r^2
*/
void mont_r_sq(mpz_t r_sq, mpz_t N) {
    mpz_set_ui(r_sq, 1); // r_sq = 1 (mod N)
    // Range: 1 to  2 * size(N) * w
    for (mp_size_t i = 1; i < 2 * mpz_size(N) * mp_bits_per_limb + 1; i++) {
        mpz_add(r_sq, r_sq, r_sq);
        // Alternative for mpz_mod        
        if (mpz_cmp(r_sq,N) >= 0) // if r_sq > N
            mpz_sub(r_sq, r_sq, N);
    }
}

/*
  Montgomery reduction
  b = mp_bits_per_limb
     Montgomery reduction, also known as REDC, is an algorithm that simultaneously computes 
     the product by R′ and reduces modulo N more quickly than the naive method. The speed is
     because all computations are done using only reduction and divisions 
     with respect to R, not N:
*/
void mont_REDC(mpz_t rop, mpz_t t, mp_limb_t omega, mpz_t N) {
    mpz_t r, bN;     // Using temp t therefore r and x can be passed as same variable (adkin to GMP functions)
    mp_limb_t u, r_i;
    mpz_init(r); mpz_init(bN);
    mpz_set(r, t);
    
    for (mp_size_t i = 0; i < mpz_size(N); i++) {
        r_i = mpz_getlimbn(r,i);
        u = r_i * omega;   // u = r_i * omega (mod b)
        // r = r + u.N.(b^i)
        mpz_mul_2exp(bN, N, mp_bits_per_limb * i);
        mpz_addmul_ui(r,bN,u);
    }

    // r = r / b^size(N)
    mpz_tdiv_q_2exp(r, r, mp_bits_per_limb * mpz_size(N) );
    
    // Alternative for mpz_mod        
    if(mpz_cmp(r,N) >= 0) // if r > N
        mpz_sub(r,r,N);
    
    mpz_swap(rop, r); 
    mpz_clear(r); mpz_clear(bN);
}

/*
    Convert to Montgomery number r, where r < N
*/
void mont_number(mpz_t r, mpz_t num, mpz_t r_sq, mp_limb_t omega, mpz_t N) {
    // r = num * rho mod N
    mont_multiplication(r, num, r_sq, omega, N);
}


/* 

Extension d - Montgomery multiplication

---------------END ------------------
*/