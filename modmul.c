#include "modmul.h"


/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c, then
- write the ciphertext c to stdout.
*/

void stage1() {

  mpz_t N, e, m, c;
  mpz_inits( N, e, m, c, NULL);
  
  while (!feof(stdin)){
    mpz_inp_str(N,stdin,16);
    mpz_inp_str(e,stdin,16);
    mpz_inp_str(m,stdin,16);

    // Requires N to be odd.
    mpz_powm_sec(c,m,e,N);
    // mpz_powm_ui??

    gmp_printf( "%ZX\n", c);
  }

  mpz_clears( N, e, m, c, NULL);
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m, then
- write the plaintext m to stdout.
*/

void stage2() {

  mpz_t N, d, p, q, dp, dq, ip, iq, c, m;
  mpz_inits( N, d, p, q, dp, dq, ip, iq, c, m, NULL);
  
  while (!feof(stdin)){
    mpz_inp_str(N,stdin,16);
    mpz_inp_str(d,stdin,16);
    mpz_inp_str(p,stdin,16);
    mpz_inp_str(q,stdin,16);
    mpz_inp_str(dp,stdin,16);
    mpz_inp_str(dq,stdin,16);
    mpz_inp_str(ip,stdin,16);
    mpz_inp_str(iq,stdin,16);
    mpz_inp_str(c,stdin,16);

    // Requires N to be odd.
    mpz_powm_sec(m,c,d,N);
    // mpz_powm_ui??

    gmp_printf( "%ZX\n", m);
  }

  mpz_clears( N, d, p, q, dp, dq, ip, iq, c, m, NULL);
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2), then
- write the ciphertext c to stdout.
*/

void stage3(int is_prng) {

  mpz_t p, q, g, h, m, c1, c2, key;
  mpz_inits(  p, q, g, h, m, c1, c2, key, NULL);
  

  while (!feof(stdin)){
    mpz_inp_str(p,stdin,16);
    mpz_inp_str(q,stdin,16);
    mpz_inp_str(g,stdin,16);
    mpz_inp_str(h,stdin,16);
    mpz_inp_str(m,stdin,16);
    
    if (is_prng == 1){
      prng(key,q);
    }
    else {
      mpz_set_ui (key, 1);
    }

    mpz_powm_sec(c1,g,key,p);

    mpz_mul (m, m, h);
    mpz_powm_sec(c2,m,key,p);

    gmp_printf( "%ZX\n", c1);
    gmp_printf( "%ZX\n", c2);

    // gmp_printf( "%zu\n", mpz_sizeinbase(q, 2));
    // gmp_printf( "%Zd\n", q);
    
  }

  mpz_clears(  p, q, g, h, m, c1, c2, key, NULL);

}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m, then
- write the plaintext m to stdout.
*/

void stage4() {

  mpz_t p, q, g, x, c1, c2, m;
  mpz_inits(  p, q, g, x, c1, c2, m, NULL);

  while (!feof(stdin)){
    mpz_inp_str(p,stdin,16);
    mpz_inp_str(q,stdin,16);
    mpz_inp_str(g,stdin,16);
    mpz_inp_str(x,stdin,16);
    mpz_inp_str(c1,stdin,16);
    mpz_inp_str(c2,stdin,16);

    mpz_powm_sec(c1,c1,x,p);
    mpz_invert (c1, c1, p);
    mpz_mul (c2, c1, c2);
    mpz_mod(m, c2, p);

    gmp_printf( "%ZX\n", m);
  }

  mpz_clears(  p, q, g, x, c1, c2, m, NULL);

}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  fill_prng_array();
  // if( 2 != argc ) {
  //   abort();
  // }

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3a" )) {
    stage3(1);
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3(0);
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
