#include "schoof.h"

/* Fonction du calcul du :
 *
 * ENTREE :
 *
 *
 *
 *
 * SORTIE :
 *  Polynôme
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

    // On traite à part les cas de f0 à f5
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
