#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>

#define printfmpz(x,y) printf("%s :", (x)); fmpz_print((y)); printf("\n");

// Déclaration de fonction
void fmpz_nextprime(fmpz_t rop, fmpz_t op);
void division_polynomial(fq_poly_t *tab, fq_t a, fq_t b, fq_poly_t ecc, ulong k, fq_ctx_t fq);
void schoof(fmpz_t card, fq_t a, fq_t b, fmpz_t q, fq_ctx_t fq);