#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

bool verbose;

// ---------------------------------------------------------------------------

double step(double x) {
  return x > 0. ? 1. : -1.;
}

double sigmoid(double x) {
  // This slope puts a fixed point (attractor) at approximately x = 0.5.
  double xcenter = 0.0, ymin = -1.0, ymax = 1.0, slope = 2.1972274554893376;
  double yscale = ymax - ymin;
  double ycenter = (ymax + ymin) / 2.0;
  double yoffset =- ycenter - (yscale / 2.0);
  return (yscale / (1.0 + exp(slope * (xcenter - x)))) + yoffset;
}

typedef struct {
  double value;
  double (*threshold_func)(double);
} Node;

typedef struct {
  int src;
  int dst;
  double weight;
} Edge;

typedef struct {
  Node *nodes;
  Edge *edges;
  int num_nodes;
  int num_edges;
  int num_in;
  int num_out;
} Genotype;

typedef struct {
  Genotype *genotype;
  double *activations;
} Organism;

// -- dot --------------------------------------------------------------------

void print_dot(Organism *o) {
  Genotype *g = o->genotype;
  double *activations = o->activations;

  printf("digraph g {\
{ rank=source edge [style=\"invis\"] n0 -> n1 }\
{ rank=sink edge [style=\"invis\"] n2 -> n3 }\n");
  for (int n = 0; n < g->num_nodes; n++)
    printf("  n%d [label=%lf]\n", n, activations[n]);
  for (int e = 0; e < g->num_edges; e++)
    printf("  n%d -> n%d [label=%lf];\n", g->edges[e].src, g->edges[e].dst,
      g->edges[e].weight);
  printf("; }\n");
}

// ---------------------------------------------------------------------------

#define EDGE_INFO_NEXT -1
#define EDGE_INFO_STOP -2

typedef struct edge_info {
  int node_idx;
  double weight;
} edge_info_t;

void init_edge_info(edge_info_t *edge_info, Genotype *g) {
  edge_info_t *ei = edge_info;
  for (int n=0; n<g->num_nodes; n++) {
    for (int e=0; e<g->num_edges; e++) {
      if (g->edges[e].dst == n) {
        ei->node_idx = e;
        ei->weight = g->edges[e].weight;
        if (verbose)
          printf("ei: %d (from %d, %4.2f)\n", n, e, ei->weight);
        ei++;
      }
    }
    ei->node_idx = EDGE_INFO_NEXT;
    ei++;
  }
  ei->node_idx = EDGE_INFO_STOP;
}

void init_in_values(double *current_values, Genotype *g) {
  for (int i=0; i<g->num_in; i++) {
    current_values[i] = g->nodes[i].value;
  }
}

void print_out_values(double *current_values, Genotype *g) {
  for (int i=g->num_in; i<g->num_in + g->num_out; i++) {
    printf("%4.2f ", current_values[i]);
  }
  printf("\n");
}

//void sa(Genotype *g, int iterations, double decay) {
void sa(Organism *o, int iterations, double decay) {
  Genotype *g = o->genotype;
  double *current_values = o->activations;

  edge_info_t edge_info[g->num_edges + g->num_nodes + 1];
  init_edge_info(edge_info, g);
  double next_values[g->num_nodes];
  bzero(next_values, sizeof(next_values));
  //double current_values[g->num_nodes];
  bzero(current_values, sizeof(current_values));
  init_in_values(current_values, g);

  while (iterations-- > 0) {
    double acc = 0.;
    double *current_value = current_values;
    double *next_value = next_values;
    edge_info_t *to_node = edge_info;
    Node *cur_node = g->nodes;
    while (to_node->node_idx != EDGE_INFO_STOP) {
      acc += current_values[to_node->node_idx] * to_node->weight;
      if (verbose)
        printf("n%ld acc=%4.2f\n", (cur_node - g->nodes), acc);
      if (to_node->node_idx == EDGE_INFO_NEXT) {
        *next_value++ = cur_node->threshold_func(decay * (*current_value++ + acc));
        cur_node++;
        acc = 0.;
      }
      to_node++;
    }
    memcpy(current_values, next_values, sizeof(next_values));
  }

  if (verbose)
    print_out_values(current_values, g);
}

// ---------------------------------------------------------------------------

double rand_edge_weight() {
  return rand() & 1 ? 1.0 : -1.0;
}

void init_random_genotype(Genotype *g, int num_edges, int num_nodes, int num_in,
        int num_out) {
  g->nodes = malloc(sizeof(Node) * num_nodes);
  g->edges = malloc(sizeof(Edge) * num_edges);
  g->num_nodes = num_nodes;
  g->num_edges = num_edges;
  g->num_in = num_in;
  g->num_out = num_out;
  for (int i=0; i<num_in; i++) {
    g->nodes[i].value =  (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1
  }
  for (int n=0; n<num_nodes; n++) {
    g->nodes[n].threshold_func = (rand() & 1) ? sigmoid : step;
  }
  for (int e=0; e<num_edges; e++) {
    g->edges[e].src = rand() % num_nodes;
    g->edges[e].dst = rand() % num_nodes;
    g->edges[e].weight = rand_edge_weight();
  }
}

void print_genotype(Genotype *g) {
  for (int e=0; e<g->num_edges; e++) {
    printf("e%d : %d -> %d\n", e, g->edges[e].src, g->edges[e].dst);
  }
}

void free_genotype(Genotype *g) {
  free(g->nodes);
  free(g->edges);
}

void init_organism(Organism *o, Genotype *g) {
  o->genotype = g;
  o->activations = calloc(sizeof(double), g->num_nodes);
}

void free_organism(Organism *o) {
  free(o->activations);
  free_genotype(o->genotype);
}

void test(int population, int num_edges, int num_nodes, int iterations,
          int epochs, double decay) {
  for (int e=0; e<epochs; e++) {
    printf("epoch %d\n", e);
    Genotype pop[population];
    for (int p=0; p<population; p++) {
      init_random_genotype(&pop[p], num_edges, num_nodes, 2, 2);
      if (verbose)
        print_genotype(&pop[p]);
    }
    Organism organisms[population];
    for (int p=0; p<population; p++) {
      init_organism(&organisms[p], &pop[p]);
      sa(&organisms[p], iterations, decay);
      print_dot(&organisms[p]);
    }
    for (int p=0; p<population; p++) {
      free_organism(&organisms[p]);
    }
  }
}

int main() {
  srand(0);
  //verbose = true;
  test(1, 30, 10, 20, 1, 1.);
  //test(100, 400, 100, 50, 100, 1.);
  return 0;
}
