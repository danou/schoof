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

/*void division_polynomial2(fq_poly_t fk, fq_poly_t ecc, fq_poly_t frob, fq_poly_t* tab, fq_ctx_t fq)
{
    fmpz_t n, tmp;
    fq_poly_t fn_2, fn_1, fn, fn1, fn2, tmp_poly;

    if(tab[k] != NULL)
    {
        fq_poly_set(fk, tab[k], fq);
    }
    else
    {
        fmpz_init(n);
        fmpz_init(tmp);

        fq_poly_init(fn_1, fq);
        fq_poly_init(fn, fq);
        fq_poly_init(fn1, fq);
        fq_poly_init(fn2, fq);
        fq_poly_init(tmp_poly, fq);

        fmpz_divexact_ui(n, k, 2);

        fmpz_sub_ui(tmp, n, 1);
        division_polynomial2(fn_1, ecc, frob, tab, tmp, fq);
        division_polynomial2(fn, ecc, frob, tab, n, fq);
        fmpz_add_ui(tmp, n, 1);
        division_polynomial2(fn1, ecc, frob, tab, tmp, fq);
        fmpz_add_ui(tmp, n, 2);
        division_polynomial2(fn2, ecc, frob, tab, tmp, fq);

        if(fmpz_is_even(k))
        {
            fq_poly_init(fn_2, fq);

            fmpz_sub_ui(tmp, n, 2);
            division_polynomial2(fn_2, ecc, frob, tab, tmp, fq);

            fq_poly_sqr(fk, fn_1, fq);
            fq_poly_mul(fk, fk, fn2, fq);
            fq_poly_sqr(tmp_poly, fn1, fq);
            fq_poly_mul(tmp_poly, tmp_poly, fn_2, fq);
            fq_poly_sub(fk, fk, tmp_poly, fq);
            fq_poly_mul(fk, fk, fn, fq);

            fq_poly_clear(fn_2, fq);
        }
        else
        {
            if(fmpz_is_even(n))
            {
                fq_poly_pow(fk, fn, 3, fq);
                fq_poly_mul(fk, fk, fn2, fq);
                fq_poly_sqr(tmp_poly, ecc, fq);
                fq_poly_mul(fk, fk, tmp_poly, fq); // c'est la réduction et élimination des y dans le polynôme de division
                fq_poly_pow(tmp_poly, fn1, 3, fq);
                fq_poly_mul(tmp_poly, tmp_poly, fn_1, fq);
                fq_poly_sub(fk, fk, tmp_poly, fq);
            }
            else
            {
                fq_poly_pow(tmp_poly, fn1, 3, fq);
                fq_poly_mul(tmp_poly, tmp_poly, fn_1, fq);
                fq_poly_sqr(fk, ecc, fq);
                fq_poly_mul(tmp_poly, tmp_poly, fk, fq); // c'est la réduction et élimination des y dans le polynôme de division
                fq_poly_pow(fk, fn, 3, fq);
                fq_poly_mul(fk, fk, fn2, fq);
                fq_poly_sub(fk, fk, tmp_poly, fq);
            }
        }
        
        fq_poly_set(tmp_poly, fk, fq);
        tab[k] = tmp_poly;
        
        fq_poly_clear(fn_1, fq); fq_poly_clear(fn, fq); fq_poly_clear(fn1, fq); fq_poly_clear(fn2, fq); //fq_poly_clear(tmp_poly, fq);
        fmpz_clear(n); fmpz_clear(tmp);

    }
}*/

void division_polynomial(fq_poly_t fk, fq_poly_t ecc, fq_poly_t frob, fq_poly_t f1, fq_poly_t f2, fq_poly_t f3, fq_poly_t f4, fmpz_t k, fq_ctx_t fq)
{
    fmpz_t n, tmp;
    fq_poly_t fn_2, fn_1, fn, fn1, fn2, tmp_poly;

    switch(fmpz_get_ui(k))
    {
        case 1 :
            fq_poly_set(fk, f1, fq); break;
            
        case 2 :
            fq_poly_set(fk, f2, fq); break;
            
        case 3 :
            fq_poly_set(fk, f3, fq); break;
            
        case 4 :
            fq_poly_set(fk, f4, fq); break;
            
        default :
            printf("cas default :"); fmpz_print(k); printf("\n");
            fmpz_init(n);
            fmpz_init(tmp);

            fq_poly_init(fn_1, fq);
            fq_poly_init(fn, fq);
            fq_poly_init(fn1, fq);
            fq_poly_init(fn2, fq);
            fq_poly_init(tmp_poly, fq);

            fmpz_divexact_ui(n, k, 2);
            fmpz_print(n); printf("\n");

            fmpz_sub_ui(tmp, n, 1); 
            division_polynomial(fn_1, ecc, frob, f1, f2, f3, f4, tmp, fq); 
            division_polynomial(fn, ecc, frob, f1, f2, f3, f4, n, fq); break;
            fmpz_add_ui(tmp, n, 1);
            division_polynomial(fn1, ecc, frob, f1, f2, f3, f4, tmp, fq);
            fmpz_add_ui(tmp, n, 2);
            division_polynomial(fn2, ecc, frob, f1, f2, f3, f4, tmp, fq);

            if(fmpz_is_even(k))
            {
                fq_poly_init(fn_2, fq);

                fmpz_sub_ui(tmp, n, 2); 
                division_polynomial(fn_2, ecc, frob, f1, f2, f3, f4, tmp, fq);

                fq_poly_sqr(fk, fn_1, fq);
                fq_poly_mul(fk, fk, fn2, fq);
                fq_poly_sqr(tmp_poly, fn1, fq);
                fq_poly_mul(tmp_poly, tmp_poly, fn_2, fq);
                fq_poly_sub(fk, fk, tmp_poly, fq);
                fq_poly_mul(fk, fk, fn, fq);

                fq_poly_clear(fn_2, fq);
            }
            else
            {
                if(fmpz_is_even(n))
                {
                    fq_poly_pow(fk, fn, 3, fq);
                    fq_poly_mul(fk, fk, fn2, fq);
                    fq_poly_sqr(tmp_poly, ecc, fq);
                    fq_poly_mul(fk, fk, tmp_poly, fq); // c'est la réduction et élimination des y dans le polynôme de division
                    fq_poly_pow(tmp_poly, fn1, 3, fq);
                    fq_poly_mul(tmp_poly, tmp_poly, fn_1, fq);
                    fq_poly_sub(fk, fk, tmp_poly, fq);
                }
                else
                {
                    fq_poly_pow(tmp_poly, fn1, 3, fq);
                    fq_poly_mul(tmp_poly, tmp_poly, fn_1, fq);
                    fq_poly_sqr(fk, ecc, fq);
                    fq_poly_mul(tmp_poly, tmp_poly, fk, fq); // c'est la réduction et élimination des y dans le polynôme de division
                    fq_poly_pow(fk, fn, 3, fq);
                    fq_poly_mul(fk, fk, fn2, fq);
                    fq_poly_sub(fk, fk, tmp_poly, fq);
                }
            }

            fq_poly_clear(fn_1, fq); fq_poly_clear(fn, fq); fq_poly_clear(fn1, fq); fq_poly_clear(fn2, fq); fq_poly_clear(tmp_poly, fq);
            fmpz_clear(n); fmpz_clear(tmp);

            break;
    }
}