Schoof's algorithm
========

It's the sources of my project.

## main.c
## function.c
## schoof.c
### sage.txt

Warning !!!!

This new version :
- the division polynomials are correct : I forgot to divide by two where k is even,
- the improvement of schoof are correct : I wrote in fq_poly_set_coeff(..., fmpz_get_si(q)^2,fq) instead of fq_poly_set_coeff(..., fmpz_get_si(q)*fmpz_get_si(q),fq),
- In the general case, I forgot to report y where f_k, f_tho, alpha, beta are even or odd. It's work for some cases.

I give a example in sage.txt where the algorithm works, but I need more time to fix the errors.
I continue to look for the errors and correct them. 


