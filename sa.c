// ideas:
// + world (contains pop, seed, options)
// + run world
// + generations per epoch
// - mutation/crossover
// - option to output dot (given organism/generation)
// - hill climb (from given organism/generation)
// - fitness function
// - change fitness function over time
// - option to output data on fitness over generations
// (if it's fast enough, maybe no need to save?)
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

int verbose = 0;
bool dot = false;

// -- graph ------------------------------------------------------------------

double step(double x) {
  return x > 0.0 ? 1.0 : -1.0;
}

double sigmoid(double x) {
  // This slope puts a fixed point (attractor) at approximately x = 0.5.
  double xcenter = 0.0, ymin = -1.0, ymax = 1.0, slope = 2.1972274554893376;
  double yscale = ymax - ymin;
  double ycenter = (ymax + ymin) / 2.0;
  double yoffset =- ycenter - (yscale / 2.0);
  double denom = (1.0 + exp(slope * (xcenter - x)));
  if (verbose > 1)
    printf("x = %lf; denom = %lf\n", x, denom);
  assert(denom != 0.0);
  return (yscale / denom) + yoffset;
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

// -- genotype ---------------------------------------------------------------

typedef struct {
  Node *nodes;
  Edge *edges;
  int num_nodes;
  int num_edges;
  int num_in;
  int num_out;
} Genotype;

double rand_edge_weight() {
  return rand() & 1 ? 1.0 : -1.0;
}

// float in range -1 to 1
double rand_value() {
  return (double)rand()/RAND_MAX*2.0-1.0;
}

void init_random_genotype(Genotype *g, int num_edges, int num_nodes, int num_in,
        int num_out) {
  assert(num_in >= 1);
  assert(num_out >= 1);
  assert(num_in + num_out <= num_nodes);
  g->nodes = malloc(sizeof(Node) * num_nodes);
  g->edges = malloc(sizeof(Edge) * num_edges);
  g->num_nodes = num_nodes;
  g->num_edges = num_edges;
  g->num_in = num_in;
  g->num_out = num_out;
  for (int i=0; i<num_in; i++) {
    g->nodes[i].value = rand_value();
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

// -- organism ---------------------------------------------------------------

typedef struct {
  Genotype *genotype;
  double *activations;
} Organism;

void init_organism(Organism *o, Genotype *g) {
  o->genotype = g;
  o->activations = calloc(sizeof(double), g->num_nodes);
}

void free_organism(Organism *o) {
  free(o->activations);
  free_genotype(o->genotype);
}

void print_organism_dot(Organism *o) {
  Genotype *g = o->genotype;
  double *activations = o->activations;

  printf("digraph g {\n");
  printf("  { rank=source edge [style=\"invis\"] ");
  for (int i = 0; i < g->num_in-1; i++)
    printf("n%d ->", i);
  printf(" n%d }\n", g->num_in-1);
  printf("  { rank=sink edge [style=\"invis\"] ");
  for (int o = 0; o < g->num_out-1; o++)
    printf("n%d ->", g->num_in + o);
  printf(" n%d }\n", g->num_in + g->num_out-1);
  for (int n = 0; n < g->num_nodes; n++)
    printf("  n%d [label=%.3lf]\n", n, activations[n]);
  for (int e = 0; e < g->num_edges; e++)
    printf("  n%d -> n%d [label=%.3lf];\n", g->edges[e].src, g->edges[e].dst,
      g->edges[e].weight);
  printf("}\n");
}

// -- spreading activation ---------------------------------------------------

#define EDGE_INFO_NEXT -1
#define EDGE_INFO_STOP -2

typedef struct {
  int node_idx;
  double weight;
} Edgeinfo;

void init_edge_info(Edgeinfo *edge_info, Genotype *g) {
  Edgeinfo *ei = edge_info;
  for (int n=0; n<g->num_nodes; n++) {
    for (int e=0; e<g->num_edges; e++) {
      if (g->edges[e].dst == n) {
        ei->node_idx = g->edges[e].src;
        ei->weight = g->edges[e].weight;
        if (verbose > 1)
          printf("ei: n%d (from %d, %4.2f)\n", n, g->edges[e].src, ei->weight);
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

void sa(Organism *o, int iterations, double decay) {
  Genotype *g = o->genotype;
  double *current_values = o->activations;
  bzero(current_values, sizeof(current_values));
  init_in_values(current_values, g);

  Edgeinfo edge_info[g->num_edges + g->num_nodes + 1];
  init_edge_info(edge_info, g);

  double next_values[g->num_nodes];
  bzero(next_values, sizeof(next_values));

  while (iterations-- > 0) {
    double acc = 0.0;
    double *current_value = current_values;
    double *next_value = next_values;
    Edgeinfo *to_node = edge_info;
    Node *cur_node = g->nodes;
    for (;;) {
      while (to_node->node_idx == EDGE_INFO_NEXT) {
        *next_value++ = cur_node->threshold_func(decay * (*current_value++ + acc));
        cur_node++;
        to_node++;
        acc = 0.0;
      }
      if (to_node->node_idx == EDGE_INFO_STOP) {
        break;
      }
      assert(to_node->node_idx >= 0 && to_node->node_idx < g->num_nodes);
      acc += current_values[to_node->node_idx] * to_node->weight;
      if (verbose > 1)
        printf("n%ld acc=%4.2lf to=%d cur=%4.2lf weight=%4.2lf\n", (cur_node - g->nodes), acc,
            to_node->node_idx, current_values[to_node->node_idx], to_node->weight);
      to_node++;
    }
    memcpy(current_values, next_values, sizeof(next_values));
  }

  if (verbose)
    print_out_values(current_values, g);
}

// -- world ------------------------------------------------------------------

typedef struct world {
  int random_seed;
  int num_organisms;
  int iterations_per_generation;
  int generations_per_epoch;
  int num_epochs;
  int num_in;
  int num_out;
  int num_nodes;
  int num_edges;
  double decay_rate;
  Genotype *genotypes;
  Organism *organisms;
} World;

World *create_world(int random_seed, int num_organisms, int iterations_per_generation,
                    int generations_per_epoch, int num_epochs,
                    int num_in, int num_out, int num_nodes, int num_edges, double decay_rate) {
  World *w = calloc(sizeof(World), 1);
  w->random_seed = random_seed;
  w->num_organisms = num_organisms;
  w->iterations_per_generation = iterations_per_generation;
  w->generations_per_epoch = generations_per_epoch;
  w->num_epochs = num_epochs;
  w->num_in = num_in;
  w->num_out = num_out;
  w->num_nodes = num_nodes;
  w->num_edges = num_edges;
  w->decay_rate = decay_rate;
  w->genotypes = calloc(sizeof(Genotype), num_organisms);
  w->organisms = calloc(sizeof(Organism), num_organisms);
  return w;
}

void init_random_population(World *w) {
  for (int n=0; n<w->num_organisms; n++) {
    init_random_genotype(&w->genotypes[n], w->num_edges, w->num_nodes, w->num_in, w->num_out);
    if (verbose > 1)
      print_genotype(&w->genotypes[n]);
    init_organism(&w->organisms[n], &w->genotypes[n]);
  }
}

void run_world(World *w) {
  srand(w->random_seed);
  for (int e=0; e<w->num_epochs; e++) {
    if (!dot)
      printf("epoch %d\n", e);
    for (int g=0; g<w->generations_per_epoch; g++) {
      if (!dot)
        printf("  generation %d\n", g);
      init_random_population(w);
      for (int n=0; n<w->num_organisms; n++) {
        sa(&w->organisms[n], w->iterations_per_generation, w->decay_rate);
        if (dot)
          print_organism_dot(&w->organisms[n]);
        free_organism(&w->organisms[n]);
      }
    }
  }
}

void free_world(World *w) {
  for (int i=0; i<w->num_organisms; i++)
    free_organism(&w->organisms[i]); // will free associated genotype
}

// ---------------------------------------------------------------------------

void quick_test() {
  verbose = 1;
  World *w = create_world(0, 1, 20, 1, 1, 2, 2, 30, 10, 1.0);
  run_world(w);
}

void dot_test() {
  dot = true;
  World *w = create_world(0, 1, 20, 1, 1, 2, 2, 30, 10, 1.0);
  run_world(w);
}

void long_test() {
  World *w = create_world(0, 50, 20, 20, 10, 2, 2, 200, 70, 1.0);
  run_world(w);
}

int main() {
  //quick_test();
  //dot_test();
  long_test();
  return 0;
}
