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
    ulong lmax, i, k, l, tho;
    mpz_t lmax_mpz, M_mpz, sqrt_mpz;
    fmpz_t M, l_fmpz, trace, sqrt, tmp, k_fmpz, tho_fmpz;
    fq_t tmp_fq, tmp1_fq, one; // constantes temporaires;
    fq_poly_t ecc; // ecc = X³ + aX + b
    fq_poly_t gcd_poly; // pgcd
    fq_poly_t frob; // frob = X^q - X
    fq_poly_t frob2; // frob2 = X^q² - X
    fq_poly_t tmp_poly;
    fq_poly_t fk, fl;
    fq_poly_t lambda;
    fq_poly_t* tab = NULL;

    // Initialisation :
    mpz_init(sqrt_mpz);
    mpz_init_set_ui(M_mpz, 2);
    mpz_init_set_ui(lmax_mpz, 3);

    // Initialisation :
    fmpz_init(trace);
    fmpz_init(sqrt);
    fmpz_init(tmp);
    fmpz_init(k_fmpz);
    fmpz_init(tho_fmpz);
    fmpz_init_set_ui(M, 2); 
    fmpz_init_set_ui(l_fmpz, 3);
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
    fq_poly_init(frob2, fq);
    fq_poly_init(tmp_poly, fq);
    fq_poly_init(fk, fq);
    fq_poly_init(fl, fq);
    fq_poly_init(lambda, fq);

    // Initialisation de l'équation de la courbe elliptique
    fq_poly_set_fq(ecc, b, fq); 
    fq_poly_set_coeff(ecc, 3, one, fq);
    fq_poly_set_coeff(ecc, 1, a, fq); 

    // Initialisation du polynôme X^q-X
    fq_poly_set_coeff(frob, fmpz_get_si(q), one, fq); // je transforme q en slong
    fq_neg(tmp_fq, tmp_fq, fq); // tmp_fq = -1
    fq_poly_set_coeff(frob, 1, tmp_fq, fq);

    // Initialisation du polynôme X^q²-X
    fq_poly_set_coeff(frob, fmpz_get_si(q) ^ 2, one, fq); // je transforme q en slong puissance 2
    fq_neg(tmp_fq, tmp_fq, fq); // tmp_fq = -1
    fq_poly_set_coeff(frob, 1, tmp_fq, fq);

    //Initialisation de lmax pour initialiser le tableau
    gmp_printf("M %Zd\n", M_mpz);
    gmp_printf("lmax %Zd\n", lmax_mpz);
    while(mpz_cmp(M_mpz, sqrt_mpz) < 0)
    {
        mpz_mul(M_mpz, M_mpz, lmax_mpz);
        gmp_printf("M %Zd\n", M_mpz);
        mpz_nextprime(lmax_mpz, lmax_mpz);
        gmp_printf("lmax %Zd\n", lmax_mpz);
    }
    lmax = 7;//mpz_get_ui(lmax_mpz);

    printf("lmax %lu\n", lmax);

    // Construction du tableau de polynôme de division
    tab = malloc((lmax +1 ) * sizeof(fq_poly_t));
    if(tab == NULL)
    {
        fprintf(stderr, "You need more memory.\n"); exit(EXIT_FAILURE);
    }
    for(i = 0; i <= lmax; i++) fq_poly_init(tab[i], fq);
    printf("lmax %lu", lmax);

    // Remplissage du tableau
    printf("polynôme de division\n");
    division_polynomial(tab, a, b, ecc, lmax, fq);

    //for(i = 0; i <= lmax ; i++) { printf("%lu : ", i) ; fq_poly_print_pretty(tab[i], "X", fq) ; printf("\n"); }


    // Cas l = 2 :
    fq_poly_gcd(gcd_poly, frob, ecc, fq);

    if(!fq_poly_is_one(gcd_poly, fq)) fmpz_set_ui(trace, 1); // si pgcd = 1 alors t = 0 [M], sinon t = 1 [M] avec ici M = 2

    // Cas général :
    while(fmpz_cmp(M, sqrt) < 0)
    {
        fmpz_mod(k_fmpz, q, l_fmpz); // k = q [l]

        k = fmpz_get_ui(k_fmpz);
        l = fmpz_get_ui(l_fmpz);

        /*
         *
         * Cas 1 : Test de (phi_l)²P = +-kP
         *
         */

        // Pré-calculs des polynôme pour le pgcd
        if(k & 0x1)
        {
            // Cas k impair
            fq_poly_sqr(gcd_poly, tab[k], fq);
            fq_poly_mul(gcd_poly, gcd_poly, frob2, fq);
            fq_poly_mul(tmp_poly, tab[k - 1], tab[k + 1], fq);
            fq_poly_mul(tmp_poly, tmp_poly, ecc, fq);
            fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
        }
        else
        {
            // Cas k pair
            fq_poly_sqr(gcd_poly, tab[k], fq);
            fq_poly_mul(gcd_poly, gcd_poly, frob2, fq);
            fq_poly_mul(gcd_poly, gcd_poly, ecc, fq);
            fq_poly_mul(tmp_poly, tab[k - 1], tab[k + 1], fq);
            fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
        }

        // Test du pgcd
        fq_poly_gcd(gcd_poly, gcd_poly, tab[l], fq);
        if(!fq_poly_is_one(gcd_poly, fq))
        {
            // Test de q carré modulo l
            if(fmpz_jacobi(q, l_fmpz) == -1)
            {
                fmpz_zero(tmp); // tmp = 0
                fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0); // On fait le théorème chinois
            }
            else
            {
                // Calcul de q = w² modulo l
                fmpz_sqrtmod(k_fmpz, q, l_fmpz); // On a mis la variable k à la place w pour libérer de la place
                k = fmpz_get_ui(k_fmpz);

                // Pré-calculs des polynôme pour le pgcd
                if(k & 0x1) // On test si w est pair
                {
                    // Cas impair
                    fq_poly_sqr(gcd_poly, tab[k], fq);
                    fq_poly_mul(gcd_poly, gcd_poly, frob, fq);
                    fq_poly_mul(tmp_poly, tab[k - 1], tab[k + 1], fq);
                    fq_poly_mul(tmp_poly, tmp_poly, ecc, fq);
                    fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                }
                else
                {
                    // Cas pair
                    fq_poly_sqr(gcd_poly, tab[k], fq);
                    fq_poly_mul(gcd_poly, gcd_poly, frob, fq);
                    fq_poly_mul(gcd_poly, gcd_poly, ecc, fq);
                    fq_poly_mul(tmp_poly, tab[k - 1], tab[k + 1], fq);
                    fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                }

                // Test du pgcd
                fq_poly_gcd(gcd_poly, gcd_poly, tab[l], fq);
                if(fq_poly_is_one(gcd_poly, fq))
                {
                    // Cas dans lequel w n'est pas une valeur propre de phi_l
                    fmpz_zero(tmp); // tmp = 0
                    fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0); // On fait le théorème chinois
                }
                else
                {
                    // Cas dans lequel phi_l(P) = +-wP

                    // Pré-calculs des polynôme pour le pgcd


                    fmpz_mul_ui(tmp, k_fmpz, 2); // tmp = 2w


                    if(k & 0x1) // On test si w est pair
                    {
                        // Cas impair
                        fq_poly_pow(gcd_poly, ecc, (fmpz_get_si(q) + 3) / 2 , fq); // gcd_poly = (X³+aX+b)^((q+3)/2)                        
                    }
                    else
                    {
                        // Cas pair
                        fq_poly_pow(gcd_poly, ecc, (fmpz_get_si(q) - 1) / 2 , fq); // gcd_poly = (X³+aX+b)^((q-1)/2)
                    }

                    fq_set_ui(tmp_fq, 4, fq); // tmp_fq = 4
                    fq_poly_scalar_mul_fq(gcd_poly, gcd_poly, tmp_fq, fq);
                    fq_poly_pow(tmp_poly, tab[k], 3, fq);
                    fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_sqr(tmp_poly, tab[k + 2], fq);
                    fq_poly_mul(tmp_poly, tmp_poly, tab[k - 1], fq);
                    fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_sqr(tmp_poly, tab[k - 2], fq);
                    fq_poly_mul(tmp_poly, tmp_poly, tab[k + 1], fq);
                    fq_poly_add(gcd_poly, gcd_poly, tmp_poly, fq);

                    // Test du pgcd
                    fq_poly_gcd(gcd_poly, gcd_poly, tab[l], fq);
                    if(fq_poly_is_one(gcd_poly, fq)) fmpz_neg(tmp, tmp); // si pgcd = 1 alors tmp = -2w [l], sinon tmp = 2w [l]

                    // On fait le théorème chinois
                    fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0);
                }
            }
        }

        /*
         *
         * Cas 2 - Test de (phi_l)²P != +-kP
         *
         */

        else
        {
            // On teste tous les tho possible tel que 0 < tho < l

            for(tho = 1; tho < l; tho++)
            {
                
                
                
                // if mod égale 0 alors break;
            }

        }

        // Incrémentation de la boucle
        fmpz_mul(M, M, l_fmpz);
        fmpz_add_ui(l_fmpz, l_fmpz, 2);
        while(fmpz_is_prime(l_fmpz) != 1) fmpz_add_ui(l_fmpz, l_fmpz, 1); // on fait la boucle jusqu'à obtenir un nombre premier, il devrait il y avoir des améliorations possible
        fmpz_print(l_fmpz);
        printf("\n");
    }

    // Retourne le résultat dans card
    fmpz_set_ui(card, 1);
    fmpz_add(card, card, q);
    fmpz_sub(card, card, trace);

    // Libération de la mémoire :
    fq_poly_clear(ecc, fq); fq_poly_clear(gcd_poly, fq); fq_poly_clear(frob, fq); fq_poly_clear(tmp_poly, fq); fq_poly_clear(lambda, fq);
    fq_poly_clear(fk, fq); fq_poly_clear(fl, fq);
    fq_clear(tmp_fq, fq); fq_clear(tmp1_fq, fq); fq_clear(one, fq);
    fmpz_clear(M); fmpz_clear(l_fmpz); fmpz_clear(sqrt); fmpz_clear(trace); fmpz_clear(tmp); fmpz_clear(k_fmpz); fmpz_clear(tho_fmpz);
    mpz_clear(sqrt_mpz); mpz_clear(lmax_mpz); mpz_clear(M_mpz);
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