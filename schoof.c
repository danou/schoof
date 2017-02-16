#include "schoof.h"

int main(int agrc, char** argv)
{
    // Déclaration
    fmpz_t q, a, b, card;
    fq_ctx_t fq;
    fq_t A, B, delta, tmp_fq;

    // Allocation
    fmpz_init(q); fmpz_init(a); fmpz_init(b); fmpz_init(card);

    // Initialisation
    if(fmpz_set_str(q, argv[1], 10)  == -1)
    {
        fprintf(stderr, "q is not an integer.\n"); exit(EXIT_FAILURE);
    }

    if(fmpz_is_prime(q) != 1)
    {
        fprintf(stderr, "q is not prime.\n"); exit(EXIT_FAILURE);
    }

    if(fmpz_set_str(a, argv[2], 10) == -1)
    {
        fprintf(stderr, "a is not an integer.\n"); exit(EXIT_FAILURE);
    }

    if(fmpz_set_str(b, argv[3], 10) == -1)
    {
        fprintf(stderr, "b is not an integer.\n"); exit(EXIT_FAILURE);
    }

    // Transformation de a et b dans Fq
    fq_ctx_init(fq, q, 1, "X");

    fq_init(A, fq); fq_init(B, fq);

    fq_set_fmpz(A, a, fq); fq_set_fmpz(B, b, fq);

    // Libération mémoire
    fmpz_clear(a); fmpz_clear(b);

    // Test du discriminant pour savoir s'il s'agit bien d'une courbe elliptique
    fq_init(delta, fq); fq_init(tmp_fq, fq);
    fq_sqr(tmp_fq, B, fq);
    fq_mul_ui(tmp_fq, tmp_fq, 27, fq);
    fq_pow_ui(delta, A, 3, fq);
    fq_mul_ui(delta, delta, 4, fq);
    fq_add(delta, delta, tmp_fq, fq);
    fq_mul_ui(delta, delta, 16, fq);
    fq_neg(delta, delta, fq);

    if(fq_is_zero(delta, fq))
    {
        fprintf(stderr, "The input parameters don't give us a elliptic curve.\n"); exit(EXIT_FAILURE);
    }

    // Libération mémoire
    fq_clear(delta, fq); fq_clear(tmp_fq, fq);

    // Fonction schoof
    schoof(card, A, B, q, fq);

    // Affichage du résultat 
    printf("Le cardinal est "); fmpz_print(card); printf("\n");

    // Libération mémoire
    fq_clear(A, fq); fq_clear(B, fq);
    fq_ctx_clear(fq);
    fmpz_clear(q); fmpz_clear(card);

    exit(EXIT_SUCCESS);
}