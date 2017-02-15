#include "modmul.h"

/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c, then
- write the ciphertext c to stdout.
*/

void stage1() {
  // Stdin parameters
  mpz_t N, e, m, c;
  mpz_init(N); mpz_init(e); mpz_init(m); mpz_init(c); 
  // Montgomery parameters
  mpz_t r_sq;
  mpz_init(r_sq);
  mp_limb_t omega;
  
  while (!feof(stdin)) {
    gmp_scanf("%Zx\n%Zx\n%Zx", N, e, m);

    // Iteration 1
    // mpz_powm_sec(c,m,e,N);

    // Iteration 2
    // squaring_exp(c,m,e,N);
    
    // Iteration 3
    // Initialise Montgomery parameters
    mont_omega(&omega, N);
    mont_r_sq(r_sq, N);
    
    // Convert m to Montgomery number
    mont_number(m, m, r_sq, omega, N);
    
    // Exponentiation: c = m^e
    window_exp(c, m, e, r_sq, omega, N);
    
    // Perform Montgomery reduction: c = c mod N
    mont_REDC(c, c, omega, N);
    

    gmp_printf( "%ZX\n", c);
  }
  mpz_clear(r_sq); mpz_clear(N); mpz_clear(e); mpz_clear(m); mpz_clear(c); 
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m, then
- write the plaintext m to stdout.
*/

void stage2() {
  // Stdin parameters
  mpz_t N, d, p, q, dP, dQ, iP, iQ, c, m;
  mpz_init(N); mpz_init(d); mpz_init(p); mpz_init(q); mpz_init(dP); 
  mpz_init(dQ); mpz_init(iP); mpz_init(iQ); mpz_init(c); mpz_init(m); 

  // Montgomery parameters
  mpz_t r_sq, mP, mQ, cP, cQ, t;
  mpz_init(r_sq); mpz_init(mP); mpz_init(mQ); mpz_init(cP); mpz_init(cQ); mpz_init(t); 
  mp_limb_t omega;
  
  while (!feof(stdin)){
    gmp_scanf("%Zx%Zx%Zx%Zx%Zx%Zx%Zx%Zx%Zx", N, d, p, q, dP, dQ, iP, iQ, c);

    // Iteration 1
    // mpz_powm_sec(m,c,d,N);
    
    // Iteration 2
    // crt(N, d, p, q, dP, dQ, iP, iQ, c, m);

    // Iteration 3
    /* 
    CRT - Garner Algorithm
    Compute: 
            mP = c^(d mod (p-1) (mod p)    ==     mP = c^(d mod phi(p)) (mod p)
            mQ = c^(d mod (q-1) (mod q)    ==     mQ = c^(d mod phi(q)) (mod q)
    */

    // Initialise Montgomery parameters p. Compute mP = c^(d mod (p-1) (mod p)
    mont_omega(&omega, p);
    mont_r_sq(r_sq, p);
    mpz_mod(cP, c, p);
    mont_number(cP, cP, r_sq, omega, p);
    window_exp(mP, cP, dP, r_sq, omega, p);
    mont_REDC(mP, mP, omega, p);
    // Initialise Montgomery parameters q. Compute mQ = c^(d mod (q-1) (mod q)
    mont_omega(&omega, q);
    mont_r_sq(r_sq, q);
    mpz_mod(cQ, c, q);
    mont_number(cQ, cQ, r_sq, omega, q);
    window_exp(mQ, cQ, dQ, r_sq, omega, q);
    mont_REDC(mQ, mQ, omega, q); 
    
    // p.t == mQ - mP (mod q)
    mpz_sub(t, mQ, mP);
    // If mQ - mP is neg, convert to positive modulo
    mpz_mod(t, t, q);
    
    // Convert t and iP to Montgomery numbers
    mont_number(t, t, r_sq, omega, q);
    mont_number(iP, iP, r_sq, omega, q);
    
    // t = p^(-1).(mQ - mP) (mod q)
    mont_multiplication(t, t, iP, omega, q);
    
    // Perform Montgomery reduction: t = t mod q
    mont_REDC(t, t, omega, q);
    
    // m = mP + p*t
    mpz_mul(m, p, t);
    mpz_add(m, m, mP);

    gmp_printf( "%ZX\n", m);

  }

  mpz_clear(N); mpz_clear(d); mpz_clear(p); mpz_clear(q); mpz_clear(dP); 
  mpz_clear(dQ); mpz_clear(iP); mpz_clear(iQ); mpz_clear(c); mpz_clear(m); 
  mpz_clear(r_sq); mpz_clear(mP); mpz_clear(mQ); mpz_clear(cP); mpz_clear(cQ); mpz_clear(t); 
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2), then
- write the ciphertext c to stdout.
*/

void stage3(int is_prng) {

  mpz_t p, q, g, h, m, c1, c2, key;
  mpz_init(p); mpz_init(q); mpz_init(g); mpz_init(h); mpz_init(m); mpz_init(c1); mpz_init(c2); mpz_init(key); 
  
  // Montgomery parameters
  mpz_t r_sq;
  mpz_init(r_sq);
  mp_limb_t omega;

  while (!feof(stdin)){
    gmp_scanf("%Zx%Zx%Zx%Zx%Zx", p, q, g, h, m);
      
    if (is_prng == 1){
      prng(key,q);
    }
    else {
      mpz_set_ui (key, 1);
    }
    
    // Iteration 1
    // mpz_powm_sec(c1,g,key,p);
    // mpz_powm_sec(c2,h,key,p);
    // mpz_mul (c2, m, c2);
    // mpz_mod(c2,c2,p);

    // Iteration 2
    // squaring_exp(c1,g,key,p);
    // squaring_exp(c2,h,key,p);
    // mpz_mul (c2, m, c2);
    // mpz_mod(c2,c2,p);



    // Enc: generate random r of size
    //      c1 = g^r    c2 = m*h^r
    
    // Initialise Montgomery parameters p.
    mont_omega(&omega, p);
    mont_r_sq(r_sq, p);
    
    // Convert g and h to Montgomery numbers
    mont_number(g, g, r_sq, omega, p);
    mont_number(h, h, r_sq, omega, p);
    
    // Compute c1 = g^key mod p
    window_exp(c1, g, key, r_sq, omega, p);
    
    // Compute c2 = h^key mod p
    window_exp(c2, h, key, r_sq, omega, p);
    
    // Convert c2 to a Montgomery number
    mont_number(c2, c2, r_sq, omega, p);
    // Compute c2 = m.c2
    mont_multiplication(c2, c2, m, omega, p);
    
    // Perform Montgomery reduction: ci = ci mod p
    mont_REDC(c2, c2, omega, p);
    mont_REDC(c1, c1, omega, p);

    gmp_printf("%ZX\n%ZX\n", c1, c2);   
  }

  mpz_clear(r_sq); mpz_clear(p); mpz_clear(q); mpz_clear(g); mpz_clear(h); mpz_clear(m); mpz_clear(c1); mpz_clear(c2); mpz_clear(key); 


}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m, then
- write the plaintext m to stdout.
*/

void stage4() {

  mpz_t p, q, g, x, c1, c2, m;
  mpz_init(p); mpz_init(q); mpz_init(g); mpz_init(x); mpz_init(c1); mpz_init(c2); mpz_init(m); 

  // Montgomery parameters
  mpz_t r_sq;
  mpz_init(r_sq);
  mp_limb_t omega;

  while (!feof(stdin)){
    gmp_scanf("%Zx%Zx%Zx%Zx%Zx%Zx", p, q, g, x, c1, c2);
    
    // Iteration 1
    // mpz_powm_sec(c1,c1,x,p);
    // mpz_invert (c1, c1, p);
    // mpz_mul (c2, c1, c2);
    // mpz_mod(m, c2, p);
    
    // Iteration 2
    // squaring_exp(c1,c1,x,p);
    // mpz_invert (c1, c1, p);
    // mpz_mul (c2, c1, c2);
    // mpz_mod(m, c2, p);


    // Initialise Montgomery parameters for p
    mont_omega(&omega, p);
    mont_r_sq(r_sq, p);
    
    // Convert c1 and c2 to Montgomery numbers
    mont_number(c1, c1, r_sq, omega, p);
    mont_number(c2, c2, r_sq, omega, p);
    
    // -x = q - x
    mpz_sub(x, q, x);
    
    // compute c1^(-x)
    window_exp(c1, c1, x, r_sq, omega, p);
    
    // c2 = c1.c2 
    mont_multiplication(c2, c2, c1, omega, p);
    
    // Perform Montgomery reduction
    mont_REDC(m, c2, omega, p);

    gmp_printf( "%ZX\n", m);
  }

  mpz_clear(r_sq); mpz_clear(p); mpz_clear(q); mpz_clear(g); mpz_clear(x); mpz_clear(c1); mpz_clear(c2); mpz_clear(m); 
}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  fill_prng_array();

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "NOPRNGstage3" )) {
    stage3(0);
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3(1);
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }
  else if( !strcmp( argv[ 1 ], "test" ) ) {
    gmp_printf( "Testing on %s ....\n", argv[2]);
    test(strcat(argv[2],".output"), "playground.output");
  }
  else {
    abort();
  }

  clear_prng_array();

  return 0;
}
