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
*/
void montgomery_multiplication(mpz_t rop, mpz_t x, mpz_t y, mp_limb_t omega, mpz_t N) {
    
    // work with a temp var instead of rop, so that the same variable can be passed as x and rop (similarly to native GMP functions)
    mpz_t r;
    mpz_init(r);
    
    // r <- 0
    mpz_set_ui(r, 0);
    mp_limb_t u, y_i, x_0, r_0;
    
    // l_N - mpz_size(N)
    for (mp_size_t i = 0; i < mpz_size(N); i++) {
        
        // u <- (r_0 + y_i*x_0)*omega (mod b)
        y_i = mpz_getlimbn(y, i); // i-th limb of y
        x_0 = mpz_getlimbn(x, 0); // 0-th limb of x
        r_0 = mpz_getlimbn(r, 0); // 0-th limb of r
        u = (r_0 + y_i * x_0) * omega;
        
        // r <- (r + y_i*x + u*N)/b
        mpz_addmul_ui(r, x, y_i); // r <- r + y_i*x
        mpz_addmul_ui(r, N, u);   // r <- r + u*N
        mpz_tdiv_q_2exp(r, r, mp_bits_per_limb);  // r <- r/b
    }
    
    // if r > N, r <- r - N
    if(mpz_cmp(r,N) >= 0)
        mpz_sub(r,r,N);
    
    mpz_swap(rop, r); // mpz_swap is O(1), while mpz_set is O(n) where n is the number of limbs
    
    // clear temp var
    mpz_clear(r);
}

/*****************************
    Montgomery preprocessing
******************************/

/*
    Computing omega
*/
void montgomery_omega(mp_limb_t* omega, mpz_t N) {
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
void montgomery_rho_sq(mpz_t rho_sq, mpz_t N) {
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
    Convert a number into a montgomery number
        num should be < N
*/
void montgomery_number(mpz_t rop, mpz_t num, mpz_t rho_sq, mp_limb_t omega, mpz_t N) {
    // r <- mont_num = num * rho (mod N)
    montgomery_multiplication(rop, num, rho_sq, omega, N);
}

/*
    Montgomery reduction
*/
void montgomery_reduction(mpz_t rop, mpz_t t, mp_limb_t omega, mpz_t N) {
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
    montgomery_multiplication(base_squared, base, base, omega, mod);
    
    // precompute look-up table
    for(int i = 1; i < m/2; i++) {
        mpz_init(table[i]);
        montgomery_multiplication(table[i], table[i-1], base_squared, omega, mod);
    }
    ////////////////////////////////////
    
    mpz_set_ui(rop,1);
    montgomery_number(rop, rop, rho_sq, omega, mod);
    
    // the number of bits
    int bits_number = mpz_size(exp) * mp_bits_per_limb;
    
    // where to start from
    int i = bits_number - 1;
    
    while(i >= 0) {
        
        // traverse through leading zeros
        if(return_bits(exp, i, i) == 0) {
            i--;
            montgomery_multiplication(rop, rop, rop, omega, mod);
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
            montgomery_multiplication(rop, rop, rop, omega, mod);
        }
        
        // table lookup
        // t <- t * T[(u-1)/2]
        montgomery_multiplication(rop, rop, table[(u - 1) / 2], omega, mod);
        
        // update i
        i = l - 1;   
    }
    
    // clear precomputed table, and base_squared
    for(int i = 0; i < m/2; i++)
        mpz_clear(table[i]);
    
    mpz_clear(base_squared);
}
