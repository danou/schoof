#include "schoof.h"

/* Fonction prochain nombre premier
 *
 * ENTREE :
 *  Entier op tel que op > 2
 * SORTIE :
 *  Entier rop
 *
*/


void fmpz_nextprime(fmpz_t rop, fmpz_t op)
{
    fmpz_add_ui(rop, op, 2);
    while(!fmpz_is_prime(rop))
    {
        fmpz_add_ui(rop, rop, 2);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Fonction polynômes de division :
 *
 * ENTREE :
 *  Tableau de Fq-polynôme tab et k la taille du tableau
 *  Entiers a,b tels que E: y² = x³ + ax + b une courbe elliptique sur Fq
 *  Fq-polynôme eec représentant la courbe elliptique
 *  Corps fini fq à q éléments
 *
 * SORTIE :
 *  Tableau de Fq-polynôme tab rempli de k polynôme de division
 *
*/

void division_polynomial(fq_poly_t *tab, fq_t a, fq_t b, fq_poly_t ecc, ulong k, fq_ctx_t fq)
{
    // Déclaration
    ulong i, n;
    fq_t tmp, tmp1;
    fq_poly_t tmp_poly;

    printf("Fonction polynôme de division\n");

    // Initialisation
    fq_init(tmp, fq);
    fq_init(tmp1, fq);

    // On traite à part les cas de f0 à f4
    // Initialisation de f0
    fq_poly_zero(tab[0], fq); // f0 = 0
    fq_poly_one(tab[1], fq); // f1 = 1
    fq_set_ui(tmp, 2, fq); // tmp = 2 
    fq_poly_set_fq(tab[2], tmp, fq); // f2 = 2

    // Initialisation de f3
    fq_set_ui(tmp, 3, fq); // tmp = 3
    fq_poly_set_coeff(tab[3], 4, tmp, fq);
    fq_mul_ui(tmp, a, 6, fq); // tmp = 6a
    fq_poly_set_coeff(tab[3], 2, tmp, fq);
    fq_mul_ui(tmp, b, 12, fq); // tmp = 12b
    fq_poly_set_coeff(tab[3], 1, tmp, fq);
    fq_sqr(tmp, a, fq); // tmp = a²
    fq_neg(tmp, tmp, fq); // tmp = -a²
    fq_poly_set_coeff(tab[3], 0, tmp, fq);

    // Initialisation de f4
    fq_mul(tmp, tmp, a, fq); // tmp = a³
    fq_sqr(tmp1, b, fq); // tmp1 = b²
    fq_mul_ui(tmp1, tmp1, 8, fq); // tmp1 = 8b²
    fq_add(tmp, tmp, tmp1, fq); // tmp = a³+8b²
    fq_mul_ui(tmp, tmp, 4, fq); // tmp = 4(a³+8b²)
    fq_neg(tmp, tmp,fq); // tmp = -4(a³+8b²)
    fq_poly_set_coeff(tab[4], 0, tmp, fq);
    fq_mul(tmp, a, b, fq); // tmp = // tmp = ab
    fq_mul_ui(tmp, tmp, 16,fq); // tmp = 16ab
    fq_neg(tmp, tmp, fq); // tmp = -16ab
    fq_poly_set_coeff(tab[4], 1, tmp, fq);
    fq_mul(tmp, a, a, fq); // tmp = a²
    fq_mul_ui(tmp, tmp, 20, fq); // tmp = 20a²
    fq_neg(tmp, tmp, fq); // tmp = -20a²
    fq_poly_set_coeff(tab[4], 2, tmp, fq); 
    fq_mul_ui(tmp, b, 80, fq); // tmp = 80b
    fq_poly_set_coeff(tab[4], 3, tmp, fq);
    fq_mul_ui(tmp, a, 20, fq); // tmp = 20a
    fq_poly_set_coeff(tab[4], 4, tmp, fq);
    fq_set_ui(tmp, 4, fq); // tmp = 4
    fq_poly_set_coeff(tab[4], 6, tmp, fq);

    // Libération mémoire
    fq_clear(tmp, fq); fq_clear(tmp1, fq);

    // Initialisation
    fq_poly_init(tmp_poly, fq);

    // Remplissage du tableau
    for(i = 5; i <= k; i++)
    {
        // Itération sur le i
        if(i & 0x1)
        {
            // Cas i impair
            n = i - 1;
            n >>= 1; // i = 2n + 1
            if(n & 0x1)
            {
                // Cas n impair
                fq_poly_pow(tmp_poly, tab[n + 1], 3, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tab[n - 1], fq);
                fq_poly_sqr(tab[i], ecc, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tab[i], fq); // c'est la réduction et élimination des y dans le polynôme de division
                fq_poly_pow(tab[i], tab[n], 3, fq);
                fq_poly_mul(tab[i], tab[i], tab[n + 2], fq);
                fq_poly_sub(tab[i], tab[i], tmp_poly, fq);
            }
            else
            {
                // Cas n pair
                fq_poly_pow(tab[i], tab[n], 3, fq);
                fq_poly_mul(tab[i], tab[i], tab[n + 2], fq);
                fq_poly_sqr(tmp_poly, ecc, fq);
                fq_poly_mul(tab[i], tab[i], tmp_poly, fq); // c'est la réduction et élimination des y dans le polynôme de division
                fq_poly_pow(tmp_poly, tab[n + 1], 3, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tab[n - 1], fq);
                fq_poly_sub(tab[i], tab[i], tmp_poly, fq);
            }
        }
        else
        {
            // Cas i pair
            n = i >> 1; // i = 2n

            fq_poly_sqr(tab[i], tab[n - 1], fq);
            fq_poly_mul(tab[i], tab [i], tab[n + 2], fq);
            fq_poly_sqr(tmp_poly, tab[n + 1], fq);
            fq_poly_mul(tmp_poly, tmp_poly, tab[n - 2], fq);
            fq_poly_sub(tab[i], tab[i], tmp_poly, fq);
            fq_poly_mul(tab[i], tab[i], tab[n], fq);
        }
    }

    // Libération de la mémoire
    fq_poly_clear(tmp_poly, fq);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
    fmpz_t M, l_fmpz, lmax_fmpz, trace, sqrt, tmp, k_fmpz, tho_fmpz;
    fq_t tmp_fq, tmp1_fq, one; // constantes temporaires;
    fq_poly_t ecc; // ecc = X³ + aX + b
    fq_poly_t gcd_poly; // pgcd
    fq_poly_t frob; // frob = X^q - X
    fq_poly_t frob2; // frob2 = X^q² - X
    fq_poly_t frob3; // frob3 = X^q² + X^q + X
    fq_poly_t tmp_poly;
    fq_poly_t tmp1_poly;
    fq_poly_t alpha;
    fq_poly_t beta;
    fq_poly_t* tab = NULL;

    // Initialisation :
    fmpz_init(trace);
    fmpz_init(sqrt);
    fmpz_init(tmp);
    fmpz_init(k_fmpz);
    fmpz_init(tho_fmpz);
    fmpz_init(lmax_fmpz);
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
    fq_poly_init(frob3, fq);
    fq_poly_init(tmp_poly, fq);
    fq_poly_init(tmp1_poly,fq);
    fq_poly_init(alpha, fq);
    fq_poly_init(beta, fq);

    // Initialisation de l'équation de la courbe elliptique
    fq_poly_set_fq(ecc, b, fq); 
    fq_poly_set_coeff(ecc, 3, one, fq);
    fq_poly_set_coeff(ecc, 1, a, fq); 

    // Initialisation du polynôme X^q-X
    fq_poly_set_coeff(frob, fmpz_get_si(q), one, fq); // je transforme q en slong
    fq_neg(tmp_fq, one, fq); // tmp_fq = -1
    fq_poly_set_coeff(frob, 1, tmp_fq, fq);

    // Initialisation du polynôme X^q²-X
    fq_poly_set_coeff(frob2, fmpz_get_si(q) ^ 2, one, fq); // je transforme q en slong puissance 2
    fq_poly_set_coeff(frob2, 1, tmp_fq, fq);

    // Initialisation du polynôme X^q² + X^q + X
    fq_poly_set_coeff(frob3, fmpz_get_si(q) ^ 2, one, fq); // je transforme q en slong puissance 2
    fq_poly_set_coeff(frob3, fmpz_get_si(q), one, fq);
    fq_poly_set_coeff(frob3, 1, one, fq);

    //Initialisation de lmax pour initialiser le tableau
    while(fmpz_cmp(M, sqrt) < 0)
    {
        fmpz_mul(M, M, l_fmpz);
        fmpz_set(lmax_fmpz, l_fmpz);
        fmpz_nextprime(l_fmpz, l_fmpz);
        if(fmpz_equal(l_fmpz, q)) fmpz_nextprime(l_fmpz, l_fmpz);
    }
    lmax = fmpz_get_ui(lmax_fmpz);
    printf("lmax = %lu\n", lmax);

    // Construction du tableau de polynôme de division
    tab = malloc((lmax + 3) * sizeof(fq_poly_t));
    if(tab == NULL)
    {
        fprintf(stderr, "You need more memory.\n"); exit(EXIT_FAILURE);
    }
    for(i = 0; i <= lmax + 2; i++) fq_poly_init(tab[i], fq);

    // Remplissage du tableau
    division_polynomial(tab, a, b, ecc, lmax + 2, fq);

    for(i = 0; i <= lmax + 2 ; i++) { printf("%lu : ", i) ; fq_poly_print_pretty(tab[i], "X", fq) ; printf("\n"); }


    // Cas l = 2 :
    fq_poly_gcd(gcd_poly, frob, ecc, fq);
    if(fq_poly_is_one(gcd_poly, fq)) fmpz_set_ui(trace, 1); // si pgcd = 1 alors t = 1 [M], sinon t = 0 [M] avec ici M = 2
    printf("t = "); fmpz_print(trace) ; printf(" [2]\n");


    // Cas général :
    fmpz_set_ui(M, 2);
    fmpz_set_ui(l_fmpz, 3);
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

        //fq_poly_print_pretty(gcd_poly, "X", fq);printf("\n");

        if(!fq_poly_is_one(gcd_poly, fq))
        {
            // Test de q carré modulo l
            if(fmpz_jacobi(q, l_fmpz) == -1)
            {
                fmpz_zero(tmp); // tmp = 0
                fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0); // On fait le théorème chinois
                printf("t = "); fmpz_print(tmp); printf(" [%lu]\n", l);
                printf("trace = "); fmpz_print(trace); printf("\n");
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
                    printf("t = "); fmpz_print(tmp); printf(" [%lu]\n", l);
                    printf("trace = "); fmpz_print(trace); printf("\n");
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
                    //printf("k = %lu\n",k);

                    // Effet de bords
                    if(k == 0x1) 
                    {
                        fq_poly_one(tmp_poly, fq);
                        fq_poly_neg(tmp_poly, tmp_poly, fq);
                        fq_poly_sqr(tmp_poly, tmp_poly, fq); // k = 1, donc k - 2 = -1 et tab[-1] = -1
                    }
                    else
                    {
                        fq_poly_sqr(tmp_poly, tab[k - 2], fq);
                    }
                    fq_poly_mul(tmp_poly, tmp_poly, tab[k + 1], fq);
                    fq_poly_add(gcd_poly, gcd_poly, tmp_poly, fq);

                    // Test du pgcd
                    fq_poly_gcd(gcd_poly, gcd_poly, tab[l], fq);
                    if(fq_poly_is_one(gcd_poly, fq)) fmpz_neg(tmp, tmp); // si pgcd = 1 alors tmp = -2w [l], sinon tmp = 2w [l]

                    // On fait le théorème chinois
                    fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0);
                    printf("t = "); fmpz_print(tmp); printf(" [%lu]\n", l);
                    printf("trace = "); fmpz_print(trace); printf("\n");
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
            printf("Cas de base\n");

            // Initialisation de alpha
            fmpz_mul(tmp, q, q);
            fmpz_add_ui(tmp, tmp, 1);
            fmpz_cdiv_q_ui(tmp, tmp, 2); //tmp = (q² + 1)/2

            fq_set_ui(tmp_fq, 4, fq); // tmp_fq = 4

            fq_poly_pow(alpha, ecc, fmpz_get_si(tmp), fq);
            fq_poly_scalar_mul_fq(alpha, alpha, tmp_fq, fq);
            fq_poly_pow(tmp_poly, tab[k], 3, fq);
            fq_poly_mul(alpha, alpha, tmp_poly, fq);
            fq_poly_sqr(tmp_poly, tab[k + 1], fq);
            fq_poly_mul(tmp_poly, tmp_poly, tab[k - 1], fq);
            fq_poly_sub(alpha, tmp_poly, alpha, fq);
            fq_poly_sqr(tmp_poly, tab[k - 1], fq);
            fq_poly_mul(tmp_poly, tmp_poly, tab[k + 2], fq);
            fq_poly_sub(alpha, tmp_poly, alpha, fq);

            // Initialisation de beta
            fq_poly_neg(beta, frob2, fq);
            fq_poly_sqr(tmp_poly, tab[k], fq);
            fq_poly_mul(beta, beta, tmp_poly, fq);
            fq_poly_mul(tmp_poly, tab[k - 1], tab[k+1], fq);
            fq_poly_sub(beta, beta, tmp_poly, fq);
            fq_poly_scalar_mul_fq(beta, beta, tmp_fq, fq);
            fq_poly_mul(beta, beta, tab[k], fq);

            // On teste tous les tho possible tel que 0 < tho < l
            for(tho = 1; tho < l/2 ; tho++)
            {
                fq_poly_mul(gcd_poly, tab[k - 1], tab[k + 1], fq);
                fq_poly_mul(gcd_poly, frob3, tab[k], fq);
                fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                fq_poly_sqr(tmp_poly, beta, fq);
                fq_poly_mul(tmp_poly, tmp_poly, ecc, fq);
                fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);

                fq_poly_sqr(tmp_poly, tab[k], fq);
                fq_poly_sqr(tmp1_poly, alpha, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                fq_poly_add(gcd_poly, gcd_poly, tmp_poly, fq);
                fq_poly_pow(tmp_poly, tab[tho], 2 * fmpz_get_si(q), fq);
                fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);

                fq_poly_sqr(tmp_poly, tab[k], fq);
                fq_poly_sqr(tmp1_poly,beta, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                fq_poly_sqr(tmp1_poly, ecc, fq);
                fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                fq_poly_pow(tmp1_poly, tab[tho - 1], fmpz_get_si(q), fq);
                fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                fq_poly_pow(tmp1_poly, tab[tho + 1], fmpz_get_si(q), fq);
                fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);


                fq_poly_rem(tmp_poly, gcd_poly, tab[l], fq);

                //fq_poly_print_pretty(tmp_poly, "X", fq);printf("\n");

                if(fq_poly_divides(tmp_poly, gcd_poly, tab[l], fq))
                {
                    fq_set_ui(tmp_fq, 2, fq);
                    fq_poly_set_coeff(gcd_poly, fmpz_get_si(q) ^ 2, tmp_fq, fq);
                    fq_poly_set_coeff(gcd_poly, 1, one, fq);
                    fq_poly_sqr(tmp_poly, tab[k], fq);
                    fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_mul(tmp_poly, tab[k - 1], tab[k + 1], fq);
                    fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_mul(gcd_poly, gcd_poly, alpha, fq);
                    fq_poly_sqr(tmp_poly, tab[k], fq);
                    fq_poly_pow(tmp1_poly, ecc, fmpz_get_si(tmp), fq);
                    fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                    fq_poly_sub(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_pow(tmp_poly, tab[tho], 3 * fmpz_get_si(q), fq);
                    fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_poly_pow(tmp_poly, tab[tho], fmpz_get_si(q) / 2, fq);
                    fq_poly_mul(gcd_poly, gcd_poly, tmp_poly, fq);
                    fq_set_ui(tmp_fq, 4, fq);
                    fq_poly_scalar_mul_fq(gcd_poly, gcd_poly, tmp_fq, fq);

                    fq_poly_sqr(tmp_poly, tab[tho - 1], fq);
                    fq_poly_mul(tmp_poly, tmp_poly, tab[tho + 2], fq);
                    fq_poly_sqr(tmp1_poly, tab[tho + 1], fq);

                    // Effet de bords
                    if(tho == 0x1)
                    {
                        fq_poly_add(tmp_poly, tmp_poly, tmp1_poly, fq); // tho = 1 , donc tho - 2 = -1 et tab[-1] = -1
                    }
                    else
                    {
                        fq_poly_mul(tmp1_poly, tmp1_poly, tab[tho - 2], fq);
                        fq_poly_sub(tmp_poly, tmp_poly, tmp1_poly, fq);
                    }
                    fq_poly_pow(tmp_poly, tmp_poly, fmpz_get_si(q), fq);
                    fq_poly_sqr(tmp1_poly, tab[k], fq);
                    fq_poly_mul(tmp_poly, tmp_poly, tmp1_poly, fq);
                    fq_poly_mul(tmp_poly, tmp_poly, beta, fq);

                    fmpz_set_ui(tmp, tho);

                    fq_poly_rem(tmp_poly, gcd_poly, tab[l], fq);

                    //fq_poly_print_pretty(tmp_poly, "X", fq);printf("\n");


                    if(fq_poly_divides(tmp_poly, gcd_poly, tab[l], fq)) fmpz_neg(tmp, tmp);

                    fmpz_CRT(trace, trace, M, tmp, l_fmpz, 0);
                    printf("t = "); fmpz_print(tmp); printf(" [%lu]\n", l);
                    printf("trace = "); fmpz_print(trace); printf("\n");
                    break;
                }
            }

        }

        // Incrémentation de la boucle
        fmpz_mul(M, M, l_fmpz);
        // on fait la boucle jusqu'à obtenir un nombre premier, il devrait il y avoir des améliorations possible
        fmpz_nextprime(l_fmpz, l_fmpz);
        if(fmpz_equal(l_fmpz, q)) fmpz_nextprime(l_fmpz, l_fmpz);
        printf("\n");
    }

    printf("trace = "); fmpz_print(trace); printf("\n");


    // Retourne le résultat dans card
    fmpz_set_ui(card, 1);
    fmpz_add(card, card, q);
    fmpz_sub(card, card, trace);

    // Libération de la mémoire du tableau :
    for(i = 0; i <= lmax + 2; i++) fq_poly_clear(tab[i], fq);
    free(tab);

    // Libération de la mémoire des variables :
    fq_poly_clear(ecc, fq); fq_poly_clear(gcd_poly, fq); fq_poly_clear(frob, fq); fq_poly_clear(frob2, fq); fq_poly_clear(frob3, fq); 
    fq_poly_clear(tmp_poly, fq); fq_poly_clear(tmp1_poly, fq); fq_poly_clear(alpha, fq); fq_poly_clear(beta, fq);
    fq_clear(tmp_fq, fq); fq_clear(tmp1_fq, fq); fq_clear(one, fq);
    fmpz_clear(M); fmpz_clear(l_fmpz); fmpz_clear(sqrt); fmpz_clear(trace); fmpz_clear(tmp); fmpz_clear(k_fmpz); fmpz_clear(tho_fmpz);
}
