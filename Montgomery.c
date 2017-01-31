#include "modmul.h"
#include <math.h>

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
void mont_multiplication(mpz_t r, mpz_t a, mpz_t b, mp_limb_t omega, mpz_t N) {
    mp_limb_t u, b_i, a_0, t_0;
    mpz_t t;
    mpz_init(t);
    mpz_set_ui(t, 0);
    
    // l_N - mpz_size(N)
    for (mp_size_t i = 0; i < mpz_size(N); i++) {
        
        // u = (t_0 + b_i.a_0).omega (mod b)
        b_i = mpz_getlimbn(b, i); // i-th limb b
        a_0 = mpz_getlimbn(a, 0); // 0-th limb a
        t_0 = mpz_getlimbn(t, 0); // 0-th limb t
        u = (t_0 + b_i * a_0) * omega;
        
        // t = (t + b_i.a + u.N) / b
        mpz_addmul_ui(t, a, b_i); // t = t + b_i.a
        mpz_addmul_ui(t, N, u);   // t = t + u.N
        mpz_tdiv_q_2exp(t, t, mp_bits_per_limb);  // t = t / b
    }
    
    // if (t > N) then t = t - N
    if(mpz_cmp(t,N) >= 0)
        mpz_sub(t,t,N);
    
    mpz_swap(r, t);
    mpz_clear(t);
}

/*****************************
    mont preprocessing
******************************/

/*
    Computing omega
*/
void mont_omega(mp_limb_t* omega, mpz_t N) {
    // omega <- 1 (mod b)
    *omega = 1;
    
    // b is the 0th limb of N, 
    mp_limb_t b = mpz_getlimbn(N, 0);
    
    for (mp_size_t i = 1; i <= mp_bits_per_limb; i++)
        // omega <- omega * omega * N (mod b)
        *omega *= *omega * b;
    
    // omega = -omega (mod b)
    *omega = -*omega;
}

/*
    Computing rho^2
*/
void mont_rho_sq(mpz_t rho_sq, mpz_t N) {
    // rho_sq <- 1 (mod N)
    mpz_set_ui(rho_sq, 1);
    
    // upto 2 * l_N * w
    for (mp_size_t i = 1; i < 2 * mpz_size(N) * mp_bits_per_limb + 1; i++) {
        // rho^2 <- rho^2 + rho^2
        mpz_add(rho_sq, rho_sq, rho_sq);
        
        // modular reduction instead of mpz_mod
        // if rho^2 > N, rho^2 <- rho^2 - N
        if (mpz_cmp(rho_sq,N) >= 0)
            mpz_sub(rho_sq, rho_sq, N);
    }
}

/*
    Convert a number into a mont number
        num should be < N
*/
void mont_number(mpz_t rop, mpz_t num, mpz_t rho_sq, mp_limb_t omega, mpz_t N) {
    // r <- mont_num = num * rho (mod N)
    mont_multiplication(rop, num, rho_sq, omega, N);
}

/*
    mont reduction
*/
void mont_reduction(mpz_t rop, mpz_t t, mp_limb_t omega, mpz_t N) {
    mpz_t r;
    mpz_init(r);
    
    // r <- t
    mpz_set(r, t);
    
    mpz_t b_times_N;
    mpz_init(b_times_N);
    mp_limb_t u, r_i;
    
    // l_N - mpz_size
    for (mp_size_t i = 0; i < mpz_size(N); i++) {
        r_i = mpz_getlimbn(r,i);
 
        // u <- r_i*omega (mod b)
        u = r_i * omega;
        
        // r <- r + (u*N*(b^i))
        mpz_mul_2exp(b_times_N, N, mp_bits_per_limb * i);
        mpz_addmul_ui(r,b_times_N,u);
    }
    
    // clear temp var
    mpz_clear(b_times_N);
    
    // r <- r / b ^ (l_N)
    mpz_tdiv_q_2exp(r, r, mp_bits_per_limb * mpz_size(N) );
    
    // if r > N, r <- r - N
    if(mpz_cmp(r,N) >= 0)
        mpz_sub(r,r,N);
    
    mpz_swap(rop, r); // mpz_swap is O(1), while mpz_set is O(n) where n is the number of limbs
    
    // clear temp var
    mpz_clear(r);
}

/*
    get window
    i - start
    k - number of bits
*/
unsigned int get_window(mpz_t exp, int i, int k) {
    mpz_t rop;
    mpz_init(rop);
    
    // rop <- exp / 2^i -- shift exp right by i
    mpz_tdiv_q_2exp(rop, exp, i);
    
    // masking off the incorrect bits
    unsigned int window = mpz_get_ui(rop) & ((1 << k) - 1);
    
    mpz_clear(rop);
    
    return window;
}

/*
    return the bits in a number from start till end,
    start and end inclusive!
*/
unsigned int return_bits(mpz_t exp, int start, int end) {
    mpz_t rop;
    mpz_init(rop);
    
    // rop <- exp / 2^end -- shift exp right until the end
    mpz_tdiv_q_2exp(rop, exp, end);
    
    // masking off the incorrect bits
    unsigned int bits = mpz_get_ui(rop) & ((1 << (start - end + 1)) - 1);
    
    mpz_clear(rop);
    
    return bits;
}

/*
    sliding window exponentiation
    rop - t
    base - x
    exp - y
    mod - n
    window - m = 2^k
*/
void window_exp(mpz_t rop, mpz_t base, mpz_t exp, mpz_t rho_sq, mp_limb_t omega, mpz_t mod) {
    int k = 8, m = 1 << k;
    mpz_t table[m/2], base_squared;
    
    ////////////////////////////////////
    // PREPROCESSING
    // base case
    mpz_init(table[0]);
    mpz_set(table[0],base);
    mpz_init(base_squared);
    
    // compute base^2 
    mont_multiplication(base_squared, base, base, omega, mod);
    
    // precompute look-up table
    for(int i = 1; i < m/2; i++) {
        mpz_init(table[i]);
        mont_multiplication(table[i], table[i-1], base_squared, omega, mod);
    }
    ////////////////////////////////////
    
    mpz_set_ui(rop,1);
    mont_number(rop, rop, rho_sq, omega, mod);
    
    // the number of bits
    int bits_number = mpz_size(exp) * mp_bits_per_limb;
    
    // where to start from
    int i = bits_number - 1;
    
    while(i >= 0) {
        
        // traverse through leading zeros
        if(return_bits(exp, i, i) == 0) {
            i--;
            mont_multiplication(rop, rop, rop, omega, mod);
            continue;
        }
        
        // end of the "initial" window
        int l = max(i - k + 1, 0);
        
        // line 8
        while(return_bits(exp, l, l) == 0)
            l++;
        
        // u is the one and true window
        int u = return_bits(exp, i, l);
        
        // line 11 (squaring t) - place adjustment
        for(int j = 0; j < i - l + 1; j++) {
            mont_multiplication(rop, rop, rop, omega, mod);
        }
        
        // table lookup
        // t <- t * T[(u-1)/2]
        mont_multiplication(rop, rop, table[(u - 1) / 2], omega, mod);
        
        // update i
        i = l - 1;   
    }
    
    // clear precomputed table, and base_squared
    for(int i = 0; i < m/2; i++)
        mpz_clear(table[i]);
    
    mpz_clear(base_squared);
}
