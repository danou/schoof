#include "schoof.h"

int push(List* list, fmpz_t k, fq_poly_t fk)
{
    fmpz_t k0;
    Element *elem = (Element*)malloc(sizeof(Element));
    if(!elem) exit(EXIT_FAILURE);
    fmpz_init_set(k0, k);
    elem->k = k0;
    elem->fk = fk;
    elem->next = *list;
    *list = elem;
}

