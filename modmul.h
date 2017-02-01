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

void modular_exponentiation(mpz_t r, mpz_t a, mpz_t b, mpz_t N);
void window_exp(mpz_t r, mpz_t b, mpz_t e, mpz_t rho_sq, mp_limb_t omega, mpz_t N);

void mont_multiplication(mpz_t r, mpz_t a, mpz_t b, mp_limb_t omega, mpz_t N);
void mont_omega(mp_limb_t* omega, mpz_t N);
void mont_r_sq(mpz_t rho_sq, mpz_t N);
void mont_number(mpz_t r, mpz_t num, mpz_t rho_sq, mp_limb_t omega, mpz_t N);
void mont_REDC(mpz_t r, mpz_t t, mp_limb_t omega, mpz_t N);
int max(int x, int y);

#endif

