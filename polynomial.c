#include "schoof.h"

int main(int agrc, char** argv)
{
    // Déclaration
    ulong k, i;
    fmpz_t q, a, b;
    fq_ctx_t fq;
    fq_t A, B, delta, tmp_fq, one;
    fq_poly_t ecc;
    fq_poly_t* tab = NULL;

    // Allocation
    fmpz_init(q); fmpz_init(a); fmpz_init(b);

    // Initialisation

   /*f(atoi(argv[0]) != 3)
    {
        fprintf(stderr, "The program need 3 inptus\n"); exit(EXIT_FAILURE);
    }*/

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

    if((k = atoi(argv[4])) == -1)
    {
        fprintf(stderr, "k is not an integer.\n"); exit(EXIT_FAILURE);
    }

    // Transformation de a et b dans Fq
    fq_ctx_init(fq, q, 1, "X");

    fq_init(A, fq); fq_init(B, fq); fq_init(one, fq);

    fq_set_fmpz(A, a, fq); fq_set_fmpz(B, b, fq); fq_one(one, fq);

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

    // Initialisation de l'équation de la courbe elliptique
    fq_poly_init(ecc, fq);
    fq_poly_set_fq(ecc, B, fq); 
    fq_poly_set_coeff(ecc, 3, one, fq);
    fq_poly_set_coeff(ecc, 1, A, fq); 

    // Construction du tableau de polynôme de division
    tab = malloc((k + 1) * sizeof(fq_poly_t));
    if(tab == NULL)
    {
        fprintf(stderr, "You need more memory.\n"); exit(EXIT_FAILURE);
    }
    for(i = 0; i <= k ; i++) fq_poly_init(tab[i], fq);

    // Remplissage du tableau
    division_polynomial(tab, A, B, ecc, k, fq);

    // Affichage du résultat 
    for(i = 0; i <= k ; i++) { printf("f_%lu : \n", i) ; fq_poly_print_pretty(tab[i], "x", fq) ; printf("\n\n"); }

    // Libération de la mémoire du tableau :
    for(i = 0; i <= k + 1; i++) fq_poly_clear(tab[i], fq);
    free(tab);

    // Libération mémoire
    fq_poly_clear(ecc, fq);
    fq_clear(A, fq); fq_clear(B, fq); fq_clear(one, fq);
    fq_ctx_clear(fq);
    fmpz_clear(q);

    exit(EXIT_SUCCESS);
}