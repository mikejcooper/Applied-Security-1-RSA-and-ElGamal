#ifndef __MODMUL_H
#define __MODMUL_H

#include  <stdio.h>
#include <stdlib.h>

#include <string.h>
#include    <gmp.h>


void test(char* file1, char* file2);
void stage1();

void prng(mpz_t rand_Num, mpz_t q);

void fill_prng_array();
void clear_prng_array();


void crt(mpz_t N, mpz_t d, mpz_t p, mpz_t q, mpz_t dP, mpz_t dQ, mpz_t ip, mpz_t iq, mpz_t c, mpz_t m);
void crt_special_case(mpz_t N, mpz_t d, mpz_t p, mpz_t q, mpz_t dP, mpz_t dQ, mpz_t ip, mpz_t iq, mpz_t c, mpz_t m);

// void sub_gmp(mpz_t a, double b);
// void div_gmp(mpz_t a, double b);

void modular_exponentiation(mpz_t r, mpz_t a, mpz_t b, mpz_t N);

void Montgomery_multiplication(mpz_t r, mpz_t a, mpz_t b, mpz_t N);


#endif

