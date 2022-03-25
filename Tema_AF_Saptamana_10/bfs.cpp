#include <stdlib.h>
#include <string.h>
#include "bfs.h"

int get_neighbors(const Grid *grid, Point p, Point neighb[])
{
    int NC = 0;
    if (p.col >= 0 && p.col < grid->cols && grid->mat[p.row][p.col - 1] == 0) // Stanga
    {
        neighb[NC].row = p.row;
        neighb[NC].col = p.col - 1;
        NC = NC + 1;
    }
    if (p.row >= 0 && p.row < grid->rows && grid->mat[p.row - 1][p.col] == 0) // Sus
    {
        neighb[NC].row = p.row - 1;
        neighb[NC].col = p.col;
        NC = NC + 1;
    }
    if (p.row >= 0 && p.row < grid->rows && grid->mat[p.row + 1][p.col] == 0) // Jos
    {
        neighb[NC].row = p.row + 1;
        neighb[NC].col = p.col;
        NC = NC + 1;
    }
    if (p.col >= 0 && p.col < grid->cols && grid->mat[p.row][p.col + 1] == 0) // Dreapta
    {
        neighb[NC].row = p.row;
        neighb[NC].col = p.col + 1;
        NC = NC + 1;
    }
    return NC;
}

void grid_to_graph(const Grid *grid, Graph *graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node *nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(grid->mat[i][j] == 0){
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }else{
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(nodes[i][j] != NULL){
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for(i=0; i<graph->nrNodes; ++i){
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if(graph->v[i]->adjSize != 0){
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for(j=0; j<graph->v[i]->adjSize; ++j){
                if( neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0){
                        graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if(k < graph->v[i]->adjSize){
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph *graph)
{
    if(graph->v != NULL){
        for(int i=0; i<graph->nrNodes; ++i){
            if(graph->v[i] != NULL){
                if(graph->v[i]->adj != NULL){
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

void enqueue(Q* Queue, Node* N)
{
    QN* NQ = (QN*)calloc(1, sizeof(QN));
    NQ->N = N;
    if (Queue->First == NULL)
    {
        Queue->First = NQ;
        Queue->Last = NQ;
    }
    else
    {
        Queue->Last->Next = NQ;
        Queue->Last = NQ;
    }
}

Node* dequeue(Q* Queue)
{
    Node* N = Queue->First->N;
    if (Queue->First == Queue->Last)
    {
        Queue->First = NULL;
        Queue->Last = NULL;
    }
    else
    {
        QN* NQ = Queue->First->Next;
        Queue->First = NQ;
    }
    return N;
}

void bfs(Graph *graph, Node *s, Operation *op)
{
    int operatii = 0;
    for (int i = 0; i < graph->nrNodes; i++)
    {
        operatii += 3;
        graph->v[i]->color = COLOR_WHITE;
        graph->v[i]->dist = 0;
        graph->v[i]->parent = NULL;
    }
    operatii += 3;
    s->color = COLOR_GRAY;
    s->dist = 0;
    s->parent = NULL;
    Q* Queue = (Q*)calloc(1, sizeof(Q));
    enqueue(Queue, s);
    operatii += 4;
    Node* U;
    operatii += 1;
    while (Queue->First != NULL)
    {
        operatii += 1;
        U = dequeue(Queue);
        operatii += 3;
        for (int i = 0; i < U->adjSize; i++)
        {
            operatii += 1;
            if (U->adj[i]->color == COLOR_WHITE)
            {
                operatii += 3;
                U->adj[i]->color = COLOR_GRAY;
                U->adj[i]->dist = U->dist + 1;
                U->adj[i]->parent = U;
                enqueue(Queue, U->adj[i]);
                operatii += 4;
            }
        }
        operatii += 1;
        U->color = COLOR_BLACK;
    }
    if (op != NULL)
    {
        op->count(operatii);
    }
}

void prettyPrint(int n, int p[], Point repr[], int ROOT, int numarSpatii)
{
    for (int j = 0; j < numarSpatii; j++)
    {
        printf("	");
    }
    numarSpatii++;
    printf("(%d, %d)\n", repr[ROOT].row, repr[ROOT].col);
    for (int i = 0; i < n; i++)
    {
        if (p[i] == ROOT)
        {
            prettyPrint(n, p, repr, i, numarSpatii);
        }
    }
}

void print_bfs_tree(Graph *graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int *p = NULL; //the parent array
    Point *repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int *transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for(int i=0; i<graph->nrNodes; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            transf[i] = n;
            ++n;
        }else{
            transf[i] = -1;
        }
    }
    if(n == 0){
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Node));
    for(int i=0; i<graph->nrNodes && !err; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            if(transf[i] < 0 || transf[i] >= n){
                err = 1;
            }else{
                repr[transf[i]] = graph->v[i]->position;
                if(graph->v[i]->parent == NULL){
                    p[transf[i]] = -1;
                }else{
                    err = 1;
                    for(int j=0; j<graph->nrNodes; ++j){
                        if(graph->v[i]->parent == graph->v[j]){
                            if(transf[j] >= 0 && transf[j] < n){
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if(!err){
        int i = 0;
        while (p[i] != -1)
        {
            i = i + 1;
        }
        prettyPrint(n, p, repr, i, 0);
    }

    if(p != NULL){
        free(p);
        p = NULL;
    }
    if(repr != NULL){
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph *graph, Node *start, Node *end, Node *path[])
{
    bfs(graph, start);
    if (end->color == COLOR_WHITE)
    {
        return -1;
    }
    int n = 0;
    while (end != start)
    {
        path[n] = end;
        n++;
        end = end->parent;
    }
    path[n] = end;
    n++;
    for (int i = 0; i < n / 2; i++)
    {
        Node* w = path[i];
        path[i] = path[n - i - 1];
        path[n - i - 1] = w;
    }
    return n;
}

void generate_random_edges(Graph* graph, int GRE)
{
    for (int i = 1; i < graph->nrNodes; i++)
    {
        graph->v[i]->adj[graph->v[i]->adjSize] = graph->v[i - 1];
        graph->v[i]->adjSize = 1;
    }
    for (int j = 0; j < GRE - graph->nrNodes + 1; j++)
    {
        int a = rand() % graph->nrNodes;
        int b;
        do
        {
            b = rand() % graph->nrNodes;
        } while (a == b);
        int w = 0;
        while (w == 0)
        {
            w = 1;
            for (int z = 0; z < graph->v[a]->adjSize; z++)
            {
                if (graph->v[a]->adj[z] == graph->v[b])
                {
                    w = 0;
                }
            }
            if (w == 0)
            {
                b = rand() % graph->nrNodes;
            }
        }
        graph->v[a]->adj[graph->v[a]->adjSize] = graph->v[b];
        graph->v[a]->adjSize = graph->v[a]->adjSize + 1;
    }
}

void performance()
{
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for(n=1000; n<=4500; n+=100){
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
            graph.v[i]->adj = (Node**)malloc(graph.nrNodes * sizeof(Node*));
            for (int j = 0; j < graph.nrNodes; j++){
                graph.v[i]->adj[j] = (Node*)malloc(sizeof(Node));
            }
        }
        generate_random_edges(&graph, n);
        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for(n=100; n<=200; n+=10){
        Operation op = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
            graph.v[i]->adj = (Node**)malloc(graph.nrNodes * sizeof(Node*));
            for (int j = 0; j < graph.nrNodes; j++) {
                graph.v[i]->adj[j] = (Node*)malloc(sizeof(Node));
            }
        }
        generate_random_edges(&graph, 4500);
        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}
