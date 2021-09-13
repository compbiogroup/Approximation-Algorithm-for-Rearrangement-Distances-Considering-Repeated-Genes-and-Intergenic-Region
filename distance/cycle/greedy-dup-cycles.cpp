#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include "greedy-dup-cycles.hpp"

/* obs: os ciclos sao sempre no que chamamos de canonico que nada mais é que o
   indice de menor valor seguido pelo seu vizinho também de menor valor */

short *copy_cycle(short *cycle) {
    short i, *copyc, tam;
    tam = cycle[0];
    copyc = (short *)malloc((tam+2)*sizeof(short));
    for (i = 0; i <= tam; i++)
        copyc[i] = cycle[i];
    return copyc;
}

short **copy_adj_list(short **adjlist, short tam) {
    short i, j, **copyadj;

    copyadj = (short **)malloc(tam * sizeof(short*));

    for (i = 0; i < tam; i++) {
        copyadj[i] = (short *)malloc((adjlist[i][0]+1) * sizeof(short));
    }
    for (i = 0; i < tam; i++) {
        for (j = 0; j <= adjlist[i][0]; j++) {
            copyadj[i][j] = adjlist[i][j];
        }
    }
    return copyadj;
}

void free_adj_list(short **adjlist, short tam) {
    short i;

    for (i = 0; i < tam; i++)
        free(adjlist[i]);

    free(adjlist);
}

void free_tuple(adj_gray_black *c, int size) {
    short i;
    for (i = 0; i < size; i++) {
        free(c[i].currc);
        free_adj_list(c[i].adj_listc, c[i].n);
        free_adj_list(c[i].adj_listb, c[i].n);
    }
    return;
}

adj_gray_black copy_tuple(short *cycle, short **adj_listc, short **adj_listb, short tam_adj) {
    adj_gray_black *cl;

    cl = (adj_gray_black *)malloc(sizeof(adj_gray_black));
    cl->currc = copy_cycle(cycle);
    cl->n = tam_adj;
    cl->adj_listc = copy_adj_list(adj_listc, tam_adj);
    cl->adj_listb = copy_adj_list(adj_listb, tam_adj);

    return *cl;
}

short find_element(short **adj_list, short el1, short el2) {
    short pos = 0, i;
    for(i = 1; i <= adj_list[el1][0]; i++) 
        if (el2 == adj_list[el1][i])
            pos = i;
    return pos;
}

/* inicia a busca em largura, considera o ciclo que esta sendo construido (que comeca com o primeiro vertice enviado pelas
    funcoes na main), as duas listas de adjacencias cinza e preta, se vai ser random ou nao e se comeca pelas cinzas ou pretas. via de regra
     eu sempre comeco pelas pretas. ele retorna uma estrutura la da segunda linha, com o ciclo que ele encontrou e o que sobrou de cinzas e pretas */
adj_gray_black *bfs_cycle(short *cycle, short **adj_listc, short **adj_listb, int tam, bool isRandom, bool is_gray) {
    adj_gray_black *list = NULL, *list2 = NULL, *result;
    short i, *list_order, first_element, last_element, next_element, aux_find, aux_min;
    int max_list2, size_list = 0, size_list2 = 0, aux_i, aux_j, aux_k, removeaux1, removeaux2, removeaux3, removeaux4, headtailcorresp1, headtailcorresp2;

    list = (adj_gray_black *)malloc(100*sizeof(adj_gray_black));
    list_order = (short *)malloc(1*sizeof(short));
    list_order[0] = 0;


    list[size_list].currc = copy_cycle(cycle);
    list[size_list].adj_listc = copy_adj_list(adj_listc, tam);
    list[size_list].adj_listb = copy_adj_list(adj_listb, tam);
    list[size_list].n = tam;
    size_list = size_list + 1;

    while(size_list > 0) {
        list2 = (adj_gray_black *)malloc(1000*sizeof(adj_gray_black));
        size_list2 = 0;
        max_list2 = 1000;

        for (aux_k = 0; aux_k < size_list && max_list2 < 20000; aux_k++) {
            i = list_order[aux_k];
            first_element = list[i].currc[1];
            last_element = list[i].currc[list[i].currc[0]];

            if (is_gray) {
                /* procura o elemento na posicao zero na lista de adjacentes do ultimo elemento
                   se for diferente de zero, entao ele achou e fechamos o ciclo */
                aux_find = find_element(list[i].adj_listc, last_element, first_element);

                /* a cycle is formed! we can stop and return this cycle */
                if (aux_find > 0) {

                    /* remove first_element and last_element list of gray edges */
                    if (list[i].adj_listc[last_element][0] > 1) {
                        for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[last_element][0]; removeaux1++) {
                            removeaux2 = list[i].adj_listc[last_element][removeaux1];
                            removeaux3 = find_element(list[i].adj_listc, removeaux2, last_element);
                            removeaux4 = list[i].adj_listc[removeaux2][0];
                            if (removeaux3 > 0) {
                                list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                            }
                        }
                        if (last_element % 2 == 1)
                            headtailcorresp1 = last_element + 1;
                        else
                            headtailcorresp1 = last_element - 1;
                        if (first_element % 2 == 1)
                            headtailcorresp2 = first_element + 1;
                        else
                            headtailcorresp2 = first_element - 1;

                        if (list[i].adj_listc[headtailcorresp1][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[headtailcorresp1][0]; removeaux1++) {
                                removeaux2 = list[i].adj_listc[headtailcorresp1][removeaux1];
                                if (removeaux2 != headtailcorresp2) {
                                    removeaux3 = find_element(list[i].adj_listc, removeaux2, headtailcorresp1);
                                    removeaux4 = list[i].adj_listc[removeaux2][0];
                                    if (removeaux3 > 0) {
                                        list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                        list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                                    }
                                }
                            }
                            list[i].adj_listc[headtailcorresp1][0] = 1;
                            list[i].adj_listc[headtailcorresp1][1] = headtailcorresp2;
                        }

                        if (list[i].adj_listc[headtailcorresp2][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[headtailcorresp2][0]; removeaux1++) {
                                removeaux2 = list[i].adj_listc[headtailcorresp2][removeaux1];
                                if (removeaux2 != headtailcorresp1) {
                                    removeaux3 = find_element(list[i].adj_listc, removeaux2, headtailcorresp2);
                                    removeaux4 = list[i].adj_listc[removeaux2][0];
                                    if (removeaux3 > 0) {
                                        list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                        list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                                    }
                                }
                            }
                            list[i].adj_listc[headtailcorresp2][0] = 1;
                            list[i].adj_listc[headtailcorresp2][1] = headtailcorresp1;
                        }

                    }
                    list[i].adj_listc[last_element][0] = 0;

                    if (list[i].adj_listc[first_element][0] > 1) {
                        for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[first_element][0]; removeaux1++) {
                            removeaux2 = list[i].adj_listc[first_element][removeaux1];
                            removeaux3 = find_element(list[i].adj_listc, removeaux2, first_element);
                            removeaux4 = list[i].adj_listc[removeaux2][0];
                            if (removeaux3 > 0) {
                                list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                            }
                        }
                        if (last_element % 2 == 1)
                            headtailcorresp2 = last_element + 1;
                        else
                            headtailcorresp2 = last_element - 1;
                        if (first_element % 2 == 1)
                            headtailcorresp1 = first_element + 1;
                        else
                            headtailcorresp1 = first_element - 1;

                        if (list[i].adj_listc[headtailcorresp1][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[headtailcorresp1][0]; removeaux1++) {
                                removeaux2 = list[i].adj_listc[headtailcorresp1][removeaux1];
                                if (removeaux2 != headtailcorresp2) {
                                    removeaux3 = find_element(list[i].adj_listc, removeaux2, headtailcorresp1);
                                    removeaux4 = list[i].adj_listc[removeaux2][0];
                                    if (removeaux3 > 0) {
                                        list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                        list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                                    }
                                }
                            }
                            list[i].adj_listc[headtailcorresp1][0] = 1;
                            list[i].adj_listc[headtailcorresp1][1] = headtailcorresp2;
                        }

                        if (list[i].adj_listc[headtailcorresp2][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list[i].adj_listc[headtailcorresp2][0]; removeaux1++) {
                                removeaux2 = list[i].adj_listc[headtailcorresp2][removeaux1];
                                if (removeaux2 != headtailcorresp1) {
                                    removeaux3 = find_element(list[i].adj_listc, removeaux2, headtailcorresp2);
                                    removeaux4 = list[i].adj_listc[removeaux2][0];
                                    if (removeaux3 > 0) {
                                        list[i].adj_listc[removeaux2][removeaux3] = list[i].adj_listc[removeaux2][removeaux4];
                                        list[i].adj_listc[removeaux2][0] = list[i].adj_listc[removeaux2][0] - 1;
                                    }
                                }
                            }
                            list[i].adj_listc[headtailcorresp2][0] = 1;
                            list[i].adj_listc[headtailcorresp2][1] = headtailcorresp1;
                        }



                    }
                    list[i].adj_listc[first_element][0] = 0;

                    /* copy and transform this cycle in the canonical representation (lowest element followed by its lowest neighbour) */
                    aux_min = 1;
                    for (aux_i = 1; aux_i <= list[i].currc[0]; aux_i++)
                        if (list[i].currc[aux_i] < list[i].currc[aux_min])
                            aux_min = aux_i;
                    result = (adj_gray_black *)malloc(sizeof(adj_gray_black));
                    result->currc = copy_cycle(list[i].currc);
                    /* if not canonical, start the cycle with element at position aux_min */
                    if (aux_min > 1) {
                        for (aux_i = aux_min; aux_i <= list[i].currc[0]; aux_i++)
                            result->currc[1+aux_i-aux_min] = list[i].currc[aux_i];
                        /* copy the elements before aux_min to the end */
                        for (aux_i = 1; aux_i < aux_min; aux_i++)
                            result->currc[list[i].currc[0]-aux_min+aux_i+1] = list[i].currc[aux_i];
                    }

                    /* reverse the list if second element is greater than the last one */
                    if (result->currc[2] > result->currc[list[i].currc[0]])
                        for (aux_i = 2; aux_i < 2+(list[i].currc[0]-1)/2; aux_i++) {
                            aux_min = result->currc[aux_i];
                            result->currc[aux_i] = result->currc[list[i].currc[0]-aux_i+2];
                            result->currc[list[i].currc[0]-aux_i+2] = aux_min;
                        }
                    /* copy the final adj_list */
                    result->adj_listc = copy_adj_list(list[i].adj_listc, list[i].n);
                    result->adj_listb = copy_adj_list(list[i].adj_listb, list[i].n);

                    /* free the current list */
                    free_tuple(list, size_list);
                    free_tuple(list2, size_list2);
                    free(list);
                    free(list2);
                    free(list_order);
                    return result;
                } else {
                    for (aux_i = 1; aux_i <= list[i].adj_listc[last_element][0]; aux_i++) {
                        list2[size_list2].currc = copy_cycle(list[i].currc);
                        list2[size_list2].adj_listc = copy_adj_list(list[i].adj_listc, list[i].n);
                        list2[size_list2].adj_listb = copy_adj_list(list[i].adj_listb, list[i].n);
                        list2[size_list2].n = list[i].n;

                        /* add next_element to the cycle list */
                        next_element = list[i].adj_listc[last_element][aux_i];
                        list2[size_list2].currc[0] = list2[size_list2].currc[0] + 1;
                        list2[size_list2].currc[list2[size_list2].currc[0]] = next_element;

                        /* remove next_element and last_element gray edges list */
                        if (list[i].adj_listc[last_element][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[last_element][0]; removeaux1++) {
                                removeaux2 = list2[size_list2].adj_listc[last_element][removeaux1];
                                removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, last_element);
                                removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                if (removeaux3 > 0) {
                                    list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                    list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                }
                            }

                            if (last_element % 2 == 1)
                                headtailcorresp1 = last_element + 1;
                            else
                                headtailcorresp1 = last_element - 1;
                            if (next_element % 2 == 1)
                                headtailcorresp2 = next_element + 1;
                            else
                                headtailcorresp2 = next_element - 1;
                            if (list2[size_list2].adj_listc[headtailcorresp1][0] > 1) {
                                for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[headtailcorresp1][0]; removeaux1++) {
                                    removeaux2 = list2[size_list2].adj_listc[headtailcorresp1][removeaux1];
                                    if (removeaux2 != headtailcorresp2) {
                                        removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, headtailcorresp1);
                                        removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                        if (removeaux3 > 0) {
                                            list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                            list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                        }
                                    }
                                }
                                list2[size_list2].adj_listc[headtailcorresp1][0] = 1;
                                list2[size_list2].adj_listc[headtailcorresp1][1] = headtailcorresp2;
                            }

                            if (list2[size_list2].adj_listc[headtailcorresp2][0] > 1) {
                                for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[headtailcorresp2][0]; removeaux1++) {
                                    removeaux2 = list2[size_list2].adj_listc[headtailcorresp2][removeaux1];
                                    if (removeaux2 != headtailcorresp1) {
                                        removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, headtailcorresp2);
                                        removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                        if (removeaux3 > 0) {
                                            list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                            list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                        }
                                    }
                                }
                                list2[size_list2].adj_listc[headtailcorresp2][0] = 1;
                                list2[size_list2].adj_listc[headtailcorresp2][1] = headtailcorresp1;
                            }

                        }
                        list2[size_list2].adj_listc[last_element][0] = 0;
                        if (list[i].adj_listc[next_element][0] > 1) {
                            for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[next_element][0]; removeaux1++) {
                                removeaux2 = list2[size_list2].adj_listc[next_element][removeaux1];
                                removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, next_element);
                                removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                if (removeaux3 > 0) {
                                    list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                    list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                }
                            }

                            if (last_element % 2 == 1)
                                headtailcorresp2 = last_element + 1;
                            else
                                headtailcorresp2 = last_element - 1;
                            if (next_element % 2 == 1)
                                headtailcorresp1 = next_element + 1;
                            else
                                headtailcorresp1 = next_element - 1;

                            if (list2[size_list2].adj_listc[headtailcorresp1][0] > 1) {
                                for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[headtailcorresp1][0]; removeaux1++) {
                                    removeaux2 = list2[size_list2].adj_listc[headtailcorresp1][removeaux1];
                                    if (removeaux2 != headtailcorresp2) {
                                        removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, headtailcorresp1);
                                        removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                        if (removeaux3 > 0) {
                                            list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                            list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                        }
                                    }
                                }
                            }

                            if (list2[size_list2].adj_listc[headtailcorresp2][0] > 1) {
                                for (removeaux1 = 1; removeaux1 <= list2[size_list2].adj_listc[headtailcorresp2][0]; removeaux1++) {
                                    removeaux2 = list2[size_list2].adj_listc[headtailcorresp2][removeaux1];
                                    if (removeaux2 != headtailcorresp1) {
                                        removeaux3 = find_element(list2[size_list2].adj_listc, removeaux2, headtailcorresp2);
                                        removeaux4 = list2[size_list2].adj_listc[removeaux2][0];
                                        if (removeaux3 > 0) {
                                            list2[size_list2].adj_listc[removeaux2][removeaux3] = list2[size_list2].adj_listc[removeaux2][removeaux4];
                                            list2[size_list2].adj_listc[removeaux2][0] = list2[size_list2].adj_listc[removeaux2][0] - 1;
                                        }
                                    }
                                }
                            }

                        }
                        list2[size_list2].adj_listc[next_element][0] = 0;


                        /* increase the size of list2 */
                        size_list2 = size_list2 + 1;

                        /* check if max_list2 is out of space; if so, double its size */
                        if (size_list2 == max_list2) {
                            max_list2 = 2*max_list2;
                            list2 = (adj_gray_black *)realloc(list2,max_list2*sizeof(adj_gray_black));
                        }
                    }
                }
            } else { /*is_black*/

                aux_find = find_element(list[i].adj_listb, last_element, first_element);

                /* a cycle is formed! we can stop and return this cycle */
                if (aux_find > 0) {

                    /* remove first_element and last_element list of gray edges */
                    list[i].adj_listb[last_element][0] = 0;
                    list[i].adj_listb[first_element][0] = 0;

                    /* copy and transform this cycle in the canonical representation (lowest element followed by its lowest neighbour) */
                    aux_min = 1;
                    for (aux_i = 1; aux_i <= list[i].currc[0]; aux_i++)
                        if (list[i].currc[aux_i] < list[i].currc[aux_min])
                            aux_min = aux_i;
                    result = (adj_gray_black *)malloc(sizeof(adj_gray_black));
                    result->currc = copy_cycle(list[i].currc);
                    /* if not canonical, start the cycle with element at position aux_min */
                    if (aux_min > 1) {
                        for (aux_i = aux_min; aux_i <= list[i].currc[0]; aux_i++)
                            result->currc[1+aux_i-aux_min] = list[i].currc[aux_i];
                        /* copy the elements before aux_min to the end */
                        for (aux_i = 1; aux_i < aux_min; aux_i++)
                            result->currc[list[i].currc[0]-aux_min+aux_i+1] = list[i].currc[aux_i];
                    }

                    /* reverse the list if second element is greater than the last one */
                    if (result->currc[2] > result->currc[list[i].currc[0]])
                        for (aux_i = 2; aux_i < 2+(list[i].currc[0]-1)/2; aux_i++) {
                            aux_min = result->currc[aux_i];
                            result->currc[aux_i] = result->currc[list[i].currc[0]-aux_i+2];
                            result->currc[list[i].currc[0]-aux_i+2] = aux_min;
                        }
                    /* copy the final adj_list */
                    result->adj_listc = copy_adj_list(list[i].adj_listc, list[i].n);
                    result->adj_listb = copy_adj_list(list[i].adj_listb, list[i].n);

                    /* free the current list */
                    free_tuple(list, size_list);
                    free_tuple(list2, size_list2);
                    free(list);
                    free(list2);
                    free(list_order);
                    return result;
                } else {
                    for (aux_i = 1; aux_i <= list[i].adj_listb[last_element][0]; aux_i++) {
                        list2[size_list2].currc = copy_cycle(list[i].currc);
                        list2[size_list2].adj_listc = copy_adj_list(list[i].adj_listc, list[i].n);
                        list2[size_list2].adj_listb = copy_adj_list(list[i].adj_listb, list[i].n);
                        list2[size_list2].n = list[i].n;

                        /* add next_element to the cycle list */
                        next_element = list[i].adj_listb[last_element][aux_i];
                        list2[size_list2].currc[0] = list2[size_list2].currc[0] + 1;
                        list2[size_list2].currc[list2[size_list2].currc[0]] = next_element;

                        /* remove next_element and last_element gray edges list */
                        list2[size_list2].adj_listb[last_element][0] = 0;
                        list2[size_list2].adj_listb[next_element][0] = 0;

                        /* increase the size of list2 */
                        size_list2 = size_list2 + 1;

                        /* check if max_list2 is out of space; if so, double its size */
                        if (size_list2 == max_list2) {
                            max_list2 = 2*max_list2;
                            list2 = (adj_gray_black *)realloc(list2,max_list2*sizeof(adj_gray_black));
                        }
                    }
                }
            }
        }
        is_gray = !(is_gray);
        free_tuple(list, size_list);
        free(list);
        list = list2;
        size_list = size_list2;

        free(list_order);
        list_order = (short *)malloc(size_list*sizeof(short));
        for (aux_j = 0; aux_j < size_list; aux_j++)
            list_order[aux_j] = aux_j;

        /*shuffle the list if isRandom is True */
        if (isRandom && size_list > 1) {
            for (aux_i = 0; aux_i < size_list; aux_i++) {
                aux_j = rand() % size_list;
                aux_k = list_order[aux_j];
                list_order[aux_j] = list_order[aux_i];
                list_order[aux_i] = aux_k;
            }
        }
    }
    printf("Ops! Something wrong happened...\n");
    exit(0);
}

/* essa funcao remove ciclos, ela funcionava em python e eu so traduzi mas nunca usei
    a ideia era basicamente remover os ciclos do Lin pra garantir a aproximacao, mas 
      eu traduzi a versao pra ciclos gerais, logo nao existem arestas pretas e cinzas
     no seu caso acho que enviar as duas listas de adjacencias e ir alterando deve ser
     suficiente. pra saber se comeca com aresta preta ou cinza basta verificar se 
     pos[1] - pos[0] = 1 ou -1, caso sim comece pelas pretas e caso contrario pelas cinzas */
bool remove_cycle(short **cycle, short ***adjlist, short tam_c, short tam_adj) {
    bool isValid = true;
    short a, b, i, j, pos, *localc, **localadj;

    localc = copy_cycle(*cycle);
    localadj = copy_adj_list(*adjlist, tam_adj);

    for (i = 0; isValid && i < tam_c; i++) {
        a = localc[i];
        b = localc[(i+1)%tam_c];
        pos = find_element(localadj, a, b);
        if (pos > 0) {
            for (j = pos; j < localadj[a][0]; j++) {
                localadj[a][j] = localadj[a][j+1];
            }
            /* quickest one */
            /* localadj[a][pos] = localadj[a][localadj[a][0]]; */
            localadj[a][0] = localadj[a][0] - 1;
        } else
            isValid = false;
        pos = find_element(localadj, b, a);
        if (pos > 0) {
            for (j = pos; j < localadj[b][0]; j++) {
                localadj[b][j] = localadj[b][j+1];
            }
            /* quickest one */
            /* localadj[a][pos] = localadj[a][localadj[a][0]]; */
            localadj[b][0] = localadj[b][0] - 1;
        } else
            isValid = false;
    }
    if (isValid) {
        free(*cycle);
        for (i = 0; i < tam_adj; i++)
            free(*adjlist[i]);
        free(*adjlist);
        cycle = &localc;
        adjlist = &localadj;
    } else {
        free(localc);
        for (i = 0; i < tam_adj; i++)
            free(localadj[i]);
        free(localadj);
    }
    return isValid;
}

short **copy_decomp(short **curr_decomp, short size_curr_decomp) {
    short i, **decomp_copy;
    decomp_copy = (short **)malloc((size_curr_decomp)*sizeof(short*));
    for (i = 0; i < size_curr_decomp; i++) {
        decomp_copy[i] = copy_cycle(curr_decomp[i]);
    }
    return decomp_copy;
}

void print_decomp(short **best_decomp, short size_best_decomp, char* out_file, double tempo) {
    short i,j, tamc;
    FILE *outfile;

    outfile = fopen (out_file, "a+");
    fprintf(outfile, "%d [", size_best_decomp);
    for (i = 0; i < size_best_decomp; i++) {
        fprintf(outfile, "[%d,", best_decomp[i][1]);
        tamc = best_decomp[i][0];
        for (j = 2; j < tamc; j++)
            fprintf(outfile, "%d,", best_decomp[i][j]);
        fprintf(outfile, "%d]", best_decomp[i][tamc]);
        if ((i+1) < size_best_decomp)
            fprintf(outfile, ",");
    }
    if (tempo < 0)
        fprintf(outfile, "]\n");
    else
        fprintf(outfile, "] %.2f\n", tempo);
    fclose(outfile);
}

void print_all_decomp(short ***all_decomps, short *sizes_all_decomps, short number_of_decomps, char* out_all_file) {
    short i;

    for (i = 0; i < number_of_decomps; i++)
        print_decomp(all_decomps[i], sizes_all_decomps[i], out_all_file, -1);
}

void free_decomp(short **best_decomp, short size_best_decomp) {
    short i;
    for (i = 0; i < size_best_decomp; i++)
        free(best_decomp[i]);
    free(best_decomp);

}

void free_all_decomp(short ***all_decomps, short *sizes_all_decomps, short number_of_decomps) {
    short i;
    for (i = 0; i < number_of_decomps; i++)
        free_decomp(all_decomps[i], sizes_all_decomps[i]);
    free(all_decomps);

}

/* int main(int argc, char* argv[]) { */
/*     short ot_sizes, *orig_perm, *target_perm, *permutation, **op_pos, **tp_pos, op_max = -1, tp_max = -1; */
/*     short vertices, edges, i, j, k, adj1, adj2, max_cycles, last_best, *available, size_available; */
/*     short allcont, *bestc, **bestadjc, **bestadjb, **curr_adj_allc, **curr_adj_allb; */
/*     short **curr_adjc, **curr_adjb, **original_adjc, **original_adjb, size_curr_decomp, **curr_decomp, size_best_decomp = 0, **best_decomp = NULL, ***all_decomp, *sizes_all_decomp, *curr_cycle; */
/*     short size_mandatory = 0, **mandatory; */
/*     char out_first_file[100] = "\0", out_file[100] = "\0", out_all_file[100] = "\0", path[50] = "/home/andre/decomp/", input_seed[15] = "\0"; */
/*     int is_random, test_all, number_of_iter, return_all_decompositions; */
/*     FILE *infile; */
/*     double time_spent; */
/*     clock_t begin, tend; */
/*     adj_gray_black *result; */
/*     char* token; */ 

/*     begin = clock(); */

/*     ot_sizes = atoi(argv[1]); */

/*     orig_perm = (short *)malloc((ot_sizes+2) * sizeof(short)); */
/*     target_perm = (short *)malloc((ot_sizes+2) * sizeof(short)); */

/*     token  = strtok(argv[2], ","); */

/*     i = 1; */
/*     orig_perm[0] = 0; */
/*     while (token != NULL) { */ 
/*         orig_perm[i] = atoi(token); */
/*         if (orig_perm[i] > op_max) */
/*             op_max = orig_perm[i]; */
/*         token = strtok(NULL, ","); */ 
/*         i = i + 1; */
/*     } */ 

/*     orig_perm[i] = op_max + 1; */

/*     token = strtok(argv[3], ","); */

/*     target_perm[0] = 0; */
/*     i = 1; */
/*     while (token != NULL) { */ 
/*         target_perm[i] = atoi(token); */
/*         if (target_perm[i] > tp_max) */
/*             tp_max = target_perm[i]; */
/*         token = strtok(NULL, ","); */  
/*         i = i + 1; */
/*     } */ 

/*     target_perm[i] = tp_max + 1; */

/*     op_pos = (short **)malloc((op_max+2) * sizeof(short*)); */
/*     tp_pos = (short **)malloc((tp_max+2) * sizeof(short*)); */


/*     for (i = 0; i < op_max+2; i++) { */
/*         op_pos[i] = (short *)malloc((op_max+2) * sizeof(short)); */
/*         op_pos[i][0] = 0; */
/*     } */
/*     for (i = 0; i < tp_max+2; i++) { */
/*         tp_pos[i] = (short *)malloc((tp_max+2) * sizeof(short)); */
/*         tp_pos[i][0] = 0; */
/*     } */

/*     for (i = 0; i < ot_sizes+2; i++) { */
/*         j = orig_perm[i]; */
/*         k = op_pos[j][0] + 1; */
/*         op_pos[j][k] = i; */
/*         op_pos[j][0] = k; */
/*     } */

/*     for (i = 0; i < ot_sizes+2; i++) { */
/*         j = target_perm[i]; */
/*         k = tp_pos[j][0] + 1; */
/*         tp_pos[j][k] = i; */
/*         tp_pos[j][0] = k; */
/*     } */

/*     vertices = 2*ot_sizes+4; */

/*     original_adjc = (short **)malloc(vertices * sizeof(short*)); */
/*     for (i = 0; i < vertices; i++) { */
/*         original_adjc[i] = (short *)malloc((ot_sizes+1) * sizeof(short)); */
/*         original_adjc[i][0] = 0; */
/*     } */

/*     original_adjb = (short **)malloc(vertices * sizeof(short*)); */
/*     for (i = 0; i < vertices; i++) { */
/*         original_adjb[i] = (short *)malloc(3*sizeof(short)); */
/*         original_adjb[i][0] = 1; */
/*         if (i % 2 == 0) */
/*             original_adjb[i][1] = i+1; */
/*         else */
/*             original_adjb[i][1] = i-1; */
/*     } */

/*     for (i = 0; i < (ot_sizes+2); i++) { */
/*         j = orig_perm[i]; */
/*         for (k = 1; k <= tp_pos[j][0]; k++) { */
/*             adj1 = original_adjc[i][0] + 1; */
/*             original_adjc[i][adj1] = ot_sizes+2+tp_pos[j][k]; */
/*             original_adjc[i][0] = adj1; */
/*         } */
/*     } */

/*     for (i = 0; i < (ot_sizes+2); i++) { */
/*         j = target_perm[i]; */
/*         for (k = 1; k <= op_pos[j][0]; k++) { */
/*             adj1 = original_adjc[i+ot_sizes+2][0] + 1; */
/*             original_adjc[i+ot_sizes+2][adj1] = op_pos[j][k]; */
/*             original_adjc[i+ot_sizes+2][0] = adj1; */
/*         } */
/*     } */



/*     strcat(out_file,argv[4]); */
/*     strcat(out_all_file,argv[4]); */
/*     strcat(out_all_file,"-all"); */
/*     strcat(out_first_file,argv[4]); */
/*     strcat(out_first_file,"-first"); */
/*     is_random = atoi(argv[5]); */
/*     test_all = atoi(argv[6]); */
/*     number_of_iter = atoi(argv[7]); */
/*     return_all_decompositions = atoi(argv[8]); */

/*     curr_adjc = copy_adj_list(original_adjc, vertices); */
/*     curr_adjb = copy_adj_list(original_adjb, vertices); */
/*     curr_decomp = (short **)malloc((vertices/2)*sizeof(short*)); */
/*     size_curr_decomp = 0; */
/*     srand(95584000+vertices*1000+ot_sizes*10+atoi(argv[9])); */

    /* Let us start with a trivial solution, using a linear variation (a.k.a FIST) */
    /* for (i = 0; i < vertices; i++) { */
    /*     while (curr_adjb[i][0] > 0) { */
    /*         curr_cycle = (short *)malloc(3*sizeof(short)); */
    /*         curr_cycle[0] = 1; */
    /*         curr_cycle[1] = i; */
    /*         result = bfs_cycle(curr_cycle, curr_adjc, curr_adjb, vertices, false, false); */
    /*         curr_decomp[size_curr_decomp] = copy_cycle(result->currc); */
    /*         size_curr_decomp = size_curr_decomp + 1; */
    /*         free(curr_cycle); */
    /*         free_adj_list(curr_adjc, vertices); */
    /*         free_adj_list(curr_adjb, vertices); */
    /*         curr_adjc = copy_adj_list(result->adj_listc, vertices); */
    /*         curr_adjb = copy_adj_list(result->adj_listb, vertices); */

    /*         free(result->currc); */
    /*         free_adj_list(result->adj_listc, vertices); */
    /*         free_adj_list(result->adj_listb, vertices); */
    /*         free(result); */
    /*     } */
    /* } */

    /* free_adj_list(curr_adjc, vertices); */
    /* free_adj_list(curr_adjb, vertices); */

    /* best_decomp = copy_decomp(curr_decomp, size_curr_decomp); */
    /* free_decomp(curr_decomp, size_curr_decomp); */
    /* size_best_decomp = size_curr_decomp; */

    /* If the FIRST variation was the goal, we can return its decomposition */
    /* if (is_random == 0 && test_all == 0) { */

    /*     tend = clock(); */
    /*     time_spent = (double)(tend - begin) / CLOCKS_PER_SEC; */

    /*     print_decomp(best_decomp, size_best_decomp, out_file, time_spent); */
    /*     free_decomp(best_decomp, size_best_decomp); */
    /*     free_adj_list(original_adjb, vertices); */
    /*     free_adj_list(original_adjc, vertices); */

    /*     free(orig_perm); */
    /*     free(target_perm); */
    /*     for (i = 0; i < op_max+2; i++) { */
    /*         free(op_pos[i]); */
    /*         free(tp_pos[i]); */
    /*     } */
    /*     free(op_pos); */
    /*     free(tp_pos); */

    /*     return(0); */
    /* } */

    /* tend = clock(); */
    /* time_spent = (double)(tend - begin) / CLOCKS_PER_SEC; */
    /* print_decomp(best_decomp, size_best_decomp, out_first_file, time_spent); */

    /* if (return_all_decompositions == 1) { */
    /*     sizes_all_decomp = (short *)malloc((number_of_iter+1)*sizeof(short)); */
    /*     all_decomp = (short ***)malloc((number_of_iter+1)*sizeof(short**)); */
    /*     sizes_all_decomp[0] = size_best_decomp; */
    /*     all_decomp[0] = copy_decomp(best_decomp, size_best_decomp); */
    /* } */

    /* max_cycles = vertices/2; */

    /* for (i = 1; i <= number_of_iter; i++) { */

    /*     curr_adjc = copy_adj_list(original_adjc, vertices); */
    /*     curr_adjb = copy_adj_list(original_adjb, vertices); */
    /*     curr_decomp = (short **)malloc((max_cycles)*sizeof(short*)); */
    /*     available = (short *)malloc(vertices * sizeof(short)); */
    /*     size_available = 0; */
    /*     size_curr_decomp = 0; */

    /*     for (j = 0; j < ot_sizes+2; j++) { */
    /*         if (curr_adjb[j][0] > 0) { */
    /*             available[size_available] = j; */
    /*             size_available = size_available + 1; */
    /*         } */
    /*     } */
    /*     last_best = 3; */

    /*     while(size_available > 0) { */

    /*         if (test_all == 1) { */
    /*             bestc = NULL; */
    /*             bestadjc = NULL; */
    /*             bestadjb = NULL; */
    /*             allcont = 1; */

    /*             for (j = 0; allcont == 1 && j < size_available; j++) { */
    /*                 curr_adj_allc = copy_adj_list(curr_adjc, vertices); */
    /*                 curr_adj_allb = copy_adj_list(curr_adjb, vertices); */
    /*                 curr_cycle = (short *)malloc(3*sizeof(short)); */
    /*                 curr_cycle[0] = 1; */
    /*                 curr_cycle[1] = available[j]; */

    /*                 result = bfs_cycle(curr_cycle, curr_adj_allc, curr_adj_allb, vertices, is_random, false); */

    /*                 if (bestc == NULL || bestc[0] < result->currc[0]) { */
    /*                     if (bestc != NULL) */
    /*                         free(bestc); */
    /*                     if (bestadjc != NULL) */
    /*                         free_adj_list(bestadjc, vertices); */
    /*                     if (bestadjb != NULL) */
    /*                         free_adj_list(bestadjb, vertices); */
    /*                     bestc = copy_cycle(result->currc); */
    /*                     bestadjc = copy_adj_list(result->adj_listc, vertices); */
    /*                     bestadjb = copy_adj_list(result->adj_listb, vertices); */
    /*                 } */
    /*                 free_adj_list(curr_adj_allc, vertices); */
    /*                 free_adj_list(curr_adj_allb, vertices); */
    /*                 free(curr_cycle); */
    /*                 if (result->currc[0] == last_best) */
    /*                     allcont = 0; */
    /*                 free(result->currc); */
    /*                 free_adj_list(result->adj_listc, vertices); */
    /*                 free_adj_list(result->adj_listb, vertices); */
    /*                 free(result); */
    /*             } */
    /*         } */
    /*         if (is_random == 1 && test_all == 0) { */
    /*             j = rand() % size_available; */
    /*             bestadjc = copy_adj_list(curr_adjc, vertices); */
    /*             bestadjb = copy_adj_list(curr_adjb, vertices); */
    /*             bestc = (short *)malloc(3*sizeof(short)); */
    /*             bestc[0] = 1; */
    /*             bestc[1] = available[j]; */
    /*             result = bfs_cycle(bestc, bestadjc, bestadjb, vertices, is_random, false); */
    /*             free(bestc); */
    /*             bestc = copy_cycle(result->currc); */
    /*             free_adj_list(bestadjc, vertices); */
    /*             bestadjc = copy_adj_list(result->adj_listc, vertices); */
    /*             free_adj_list(bestadjb, vertices); */
    /*             bestadjb = copy_adj_list(result->adj_listb, vertices); */
    /*             free(result->currc); */
    /*             free_adj_list(result->adj_listc, vertices); */
    /*             free_adj_list(result->adj_listb, vertices); */
    /*             free(result); */
    /*         } */
    /*         free_adj_list(curr_adjc, vertices); */
    /*         free_adj_list(curr_adjb, vertices); */
    /*         curr_adjc = copy_adj_list(bestadjc, vertices); */
    /*         curr_adjb = copy_adj_list(bestadjb, vertices); */
    /*         free_adj_list(bestadjc, vertices); */
    /*         free_adj_list(bestadjb, vertices); */
    /*         curr_decomp[size_curr_decomp] = copy_cycle(bestc); */
    /*         size_curr_decomp = size_curr_decomp + 1; */
    /*         last_best = bestc[0]; */
    /*         free(bestc); */


    /*         size_available = 0; */
    /*         for (j = 0; j < (ot_sizes+2); j++) { */
    /*             if (curr_adjb[j][0] > 0) { */
    /*                 available[size_available] = j; */
    /*                 size_available = size_available + 1; */
    /*             } */
    /*         } */

    /*     } */

    /*     if (size_best_decomp < size_curr_decomp) { */
    /*         free_decomp(best_decomp, size_best_decomp); */
    /*         best_decomp = copy_decomp(curr_decomp, size_curr_decomp); */
    /*         size_best_decomp = size_curr_decomp; */
    /*     } */

    /*     if (return_all_decompositions == 1) { */
    /*         /1* cc = mandatory_valid + final_cycles *1/ */
    /*         sizes_all_decomp[i] = size_curr_decomp; */
    /*         all_decomp[i] = copy_decomp(curr_decomp, size_curr_decomp); */


    /*     } */
    /*     free_decomp(curr_decomp, size_curr_decomp); */
    /*     free(available); */
    /*     free_adj_list(curr_adjc, vertices); */
    /*     free_adj_list(curr_adjb, vertices); */
    /* } */

    /* free_adj_list(original_adjb, vertices); */
    /* free_adj_list(original_adjc, vertices); */

    /* cc = best_final_cycles + mandatory_valid
    tamcc = len(cc) */
    /*return tamcc, str(cc).replace(" ",""), all_decomps, tend-tstart*/

    /* if (return_all_decompositions == 1) { */
    /*     print_all_decomp(all_decomp, sizes_all_decomp, number_of_iter+1, out_all_file); */
    /*     free_all_decomp(all_decomp, sizes_all_decomp, number_of_iter+1); */
    /*     free(sizes_all_decomp); */
    /* } */

    /* tend = clock(); */
    /* time_spent = (double)(tend - begin) / CLOCKS_PER_SEC; */

    /* print_decomp(best_decomp, size_best_decomp, out_file, time_spent); */
    /* free_decomp(best_decomp, size_best_decomp); */

    /* free(orig_perm); */
    /* free(target_perm); */
    /* for (i = 0; i < op_max+2; i++) { */
    /*     free(op_pos[i]); */
    /*     free(tp_pos[i]); */
    /* } */
    /* free(op_pos); */
    /* free(tp_pos); */
    
    /* return(0); */
/* } */

