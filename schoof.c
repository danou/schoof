#include "schoof.h"

/* Algorithme de Schoof :
 *
 * ENTREE :
 *  Entier q premier tel que Fq un corps fini à q éléments
 *  Entiers a,b tels que E: y² = x³ + ax + b une courbe elliptique sur Fq
 *
 * SORTIE :
 *  Entier card tel que #E(Fq) = card
 * 
*/ 

void schoof(fmpz_t card, fq_t a, fq_t b, fmpz_t q, fq_ctx_t fq)
{

    // Déclaration :
    fmpz_t M, l, trace, sqrt, tmp, k;
    fq_t tmp_fq, tmp1_fq, one; // constantes temporaires;
    fq_poly_t ecc; // ecc = X³ + aX + b
    fq_poly_t gcd_poly; // pgcd
    fq_poly_t frob; // frob = X^n - x
    fq_poly_t tmp_poly;
    fq_poly_t f0, f1, f2, f3, f4, f5, fk, fl;
    
    
    
    
    
    
    FILE* data = NULL;
    
    data = fopen("test.txt", "wr");
    
    
    
    
    
    
    
    

    // Initialisation :
    fmpz_init(trace);
    fmpz_init(sqrt);
    fmpz_init(tmp);
    fmpz_init(k);
    fmpz_init_set_ui(M, 2); 
    fmpz_init_set_ui(l, 3);
    fmpz_sqrt(sqrt, q);
    fmpz_mul_ui(sqrt, sqrt, 4);

    fq_init(tmp_fq, fq);
    fq_init(tmp1_fq, fq);
    fq_init(one, fq);
    fq_one(tmp_fq, fq);
    fq_one(one, fq);

    fq_poly_init(ecc, fq);
    fq_poly_init(gcd_poly, fq);
    fq_poly_init(frob, fq);
    fq_poly_init(tmp_poly, fq);
    fq_poly_init(f0, fq);
    fq_poly_init(f1, fq);
    fq_poly_init(f2, fq);
    fq_poly_init(f3, fq);
    fq_poly_init(f4, fq);
    fq_poly_init(f5, fq);
    fq_poly_init(fk, fq);
    fq_poly_init(fl, fq);

    // Initialisation de l'équation de la courbe elliptique
    fq_poly_set_fq(ecc, b, fq); 
    fq_poly_set_coeff(ecc, 3, one, fq);
    fq_poly_set_coeff(ecc, 1, a, fq); 

    // Initialisation du polynome X^q-X
    fq_poly_set_coeff(frob, fmpz_get_si(q), one, fq); // je transforme q en slong
    fq_neg(tmp_fq, tmp_fq, fq); // tmp_fq = -1
    fq_poly_set_coeff(frob, 1, tmp_fq, fq);

    // Initialisation de f0
    fq_poly_zero(f0, fq); // f0 = 0
    fq_poly_one(f1, fq);
    fq_set_ui(tmp_fq, 2, fq); 
    fq_poly_set_fq(f2, tmp_fq, fq); // tmp_fq = 2 et f2 = 2

    // Initialisation de f3
    fq_set_ui(tmp_fq, 3, fq); // tmp_fq = 3
    fq_poly_set_coeff(f3, 4, tmp_fq, fq);
    fq_mul_ui(tmp_fq, a, 6, fq); // tmp_fq = 6a
    fq_poly_set_coeff(f3, 2, tmp_fq, fq);
    fq_mul_ui(tmp_fq, b, 12, fq); // tmp_fq = 12b
    fq_poly_set_coeff(f3, 1, tmp_fq, fq);
    fq_mul(tmp_fq, a, a, fq); // tmp_fq = a²
    fq_poly_set_coeff(f3, 0, tmp_fq, fq);

    // Initialisation de f4
    fq_mul(tmp_fq, tmp_fq, a, fq); // tmp_fq = a³
    fq_mul(tmp1_fq, b, b, fq); // tmp1_fq = b²
    fq_mul_ui(tmp1_fq, tmp1_fq, 8, fq); // tmp1_fq = 8b²
    fq_add(tmp_fq, tmp_fq, tmp1_fq, fq); // tmp_fq = a³+8b²
    fq_mul_ui(tmp_fq, tmp_fq, 4, fq); // tmp_fq = 4(a³+8b²)
    fq_neg(tmp_fq, tmp_fq,fq); // tmp_fq = -4(a³+8b²)
    fq_poly_set_coeff(f4, 0, tmp_fq, fq);
    fq_mul(tmp_fq, a, b, fq); // tmp_fq = // tmp_fq = ab
    fq_mul_ui(tmp_fq, tmp_fq, 16,fq); // tmp_fq = 16ab
    fq_neg(tmp_fq, tmp_fq, fq); // tmp_fq = -16ab
    fq_poly_set_coeff(f4, 1, tmp_fq, fq);
    fq_mul(tmp_fq, a, a, fq); // tmp_fq = a²
    fq_mul_ui(tmp_fq, tmp_fq, 20, fq); // tmp_fq = 20a²
    fq_neg(tmp_fq, tmp_fq, fq); // tmp_fq = -20a²
    fq_poly_set_coeff(f4, 2, tmp_fq, fq); 
    fq_mul_ui(tmp_fq, b, 80, fq); // tmp_fq = 80b
    fq_poly_set_coeff(f4, 3, tmp_fq, fq);
    fq_mul_ui(tmp_fq, a, 20, fq); // tmp_fq = 20a
    fq_poly_set_coeff(f4, 4, tmp_fq, fq);
    fq_set_ui(tmp_fq, 4, fq); // tmp_fq = 4
    fq_poly_set_coeff(f4, 6, tmp_fq, fq); 

    // Initialisation de f5
    fq_poly_pow(f5, f2, 3, fq);
    fq_poly_mul(f5, f5, f4, fq);
    fq_poly_sqr(tmp_poly, ecc, fq);
    fq_poly_mul(f5, f5, tmp_poly, fq); // c'est la réduction et élimination des y dans le polynôme de division
    fq_poly_pow(tmp_poly, f3, 3, fq);
    fq_poly_mul(tmp_poly, tmp_poly, f1, fq);
    fq_poly_sub(f5, f5, tmp_poly, fq);
    
    fq_poly_print_pretty(f5, "X", fq); 
    
    
    
    
    
    //fq_poly_fprint(data, f0, fq);
    //fq_poly_fprint(data, f1, fq);
    //fq_poly_fprint(data, f2, fq);
    //fq_poly_fprint(data, f3, fq);
    //fq_poly_fprint(data, f4, fq);
    
    fq_poly_fprint_pretty(data, f3, "X", fq);

    
    
    

    fclose(data);
    exit(-1);
    
    
    
    
    
    
    
    
    
    
    
    
    

    // Cas t modulo 2 :
    fq_poly_gcd(gcd_poly, frob, ecc, fq);

    if(!fq_poly_is_one(gcd_poly, fq)) fmpz_set_ui(trace, 1); // si pgcd=1 alors t=0[M] sinon t=1[M] avec ici M = 2

    // Cas général :
    while(fmpz_cmp(M, sqrt) < 0)
    {
        //printf("Une boucle while\n");

        fmpz_mod(k, q, l);

        // Test de (phi_l)²P = +-k>P

        // Précalculs des polynôme pour le pgcd
        if(fmpz_is_even(k)) 
        {
            division_polynomial(fk, ecc, frob, f1, f2, f3, f4, k, fq);
            fq_poly_pow(gcd_poly, fk, 2, fq);
            fq_poly_mul(gcd_poly, gcd_poly, frob, fq);
            fq_poly_mul(gcd_poly, gcd_poly, ecc, fq);
            fmpz_sub_ui(k, k, 1);
            division_polynomial(tmp_poly, ecc, frob, f1, f2, f3, f4, k, fq);
            fmpz_add_ui(k, k, 2);
            division_polynomial(fk, ecc, frob, f1, f2, f3, f4, k, fq);
            fq_poly_mul(tmp_poly, tmp_poly, fk, fq);
            fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
        }
        else
        {
            division_polynomial(fk, ecc, frob, f1, f2, f3, f4, k, fq);
            fq_poly_pow(gcd_poly, fk, 2, fq);
            fq_poly_mul(gcd_poly, gcd_poly, frob, fq);
            fmpz_sub_ui(k, k, 1);
            division_polynomial(tmp_poly, ecc, frob, f1, f2, f3, f4, k, fq);
            fmpz_add_ui(k, k, 2);
            division_polynomial(fk, ecc, frob, f1, f2, f3, f4, k, fq);
            fq_poly_mul(tmp_poly, tmp_poly, fk, fq);
            fq_poly_mul(tmp_poly, tmp_poly, ecc, fq);
            fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
        }

        division_polynomial(fl, ecc, frob, f1, f2, f3, f4, l, fq);
        fq_poly_print_pretty(fl, "X", fq);
        printf("\n");
        fmpz_print(l);
        printf("\n");

        // Test du pgcd
        fq_poly_gcd(gcd_poly, gcd_poly, ecc, fq);
        if(!fq_poly_is_one(gcd_poly, fq))
        {
            // Test de q carré modulo l
            if(fmpz_jacobi(q, l) == -1)
            {
                fmpz_zero(tmp); // tmp = 0
                fmpz_CRT(trace, trace, M, tmp, l, 0); // on fait le théorème chinois
            }
            else
            {
            }
        }
        else
        {

        }

        // Incrémentation de la boucle
        fmpz_mul(M, M, l);
        fmpz_add_ui(l, l, 2);
        while(fmpz_is_prime(l) != 1) fmpz_add_ui(l, l, 1); // on fait la boucle jusqu'à obtenir un nombre premier, il devrait il y avoir des améliorations possible
        fmpz_print(l);
        printf("\n");
    }

    // Retourne le résultat dans card
    fmpz_set_ui(card, 1);
    fmpz_add(card, card, q);
    fmpz_sub(card, card, trace);

    // Libération de la mémoire :
    fq_poly_clear(ecc, fq); fq_poly_clear(gcd_poly, fq); fq_poly_clear(frob, fq); fq_poly_clear(tmp_poly, fq);
    fq_poly_clear(f0, fq); fq_poly_clear(f1, fq); fq_poly_clear(f2, fq); fq_poly_clear(f3, fq); fq_poly_clear(f4, fq);
    fq_poly_clear(fk, fq); fq_poly_clear(fl, fq);
    fq_clear(tmp_fq, fq); fq_clear(tmp1_fq, fq); fq_clear(one, fq);
    fmpz_clear(M); fmpz_clear(l); fmpz_clear(sqrt); fmpz_clear(trace); fmpz_clear(tmp), fmpz_clear(k);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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