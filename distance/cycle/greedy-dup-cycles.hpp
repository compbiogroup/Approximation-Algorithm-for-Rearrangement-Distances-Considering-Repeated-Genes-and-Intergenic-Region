#pragma once

/* estrutura que salva um ciclo atual, o numero de vertices que existem (sempre 2n+4),
    e as listas de adjacencias entre vertices considerando arestas cinzas e arestas pretas */
typedef struct adj_gray_black {
    short *currc;
    short n;
    short **adj_listc;
    short **adj_listb;
} adj_gray_black;

void free_adj_list(short **adjlist, short tam);

adj_gray_black *bfs_cycle(short *cycle, short **adj_listc, short **adj_listb, int tam, bool isRandom, bool is_gray);
