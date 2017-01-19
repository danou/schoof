#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>

//Structure de liste chainée
struct Element
{
    fq_poly_t fk;
    fmpz_t k;
    struct Element * next;
};
 
typedef struct Element  Element;
typedef Element*  List;


// Déclaration de fonction
void division_polynomial(fq_poly_t fk, fq_poly_t ecc, fq_poly_t frob, fq_poly_t f1, fq_poly_t f2, fq_poly_t f3, fq_poly_t f4, fmpz_t k, fq_ctx_t fq);
void schoof(fmpz_t card, fq_t a, fq_t b, fmpz_t q, fq_ctx_t fq);