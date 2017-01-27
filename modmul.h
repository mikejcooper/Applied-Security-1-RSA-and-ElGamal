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


#endif

