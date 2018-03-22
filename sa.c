// todo:
// + world (contains pop, seed, options)
// + run world
// + generations per epoch
// + fitness function
// - mutation
// - crossover
// - option to output dot (given organism/generation)
// - hill climb (from given organism/generation)
// - change fitness function over time
// - option to output data on fitness over generations
// - command line args (or separate mains)
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

double clamp(double x) {
  return x <= -1.0 ? -1.0 : (x >= 1.0 ? 1.0 : x);
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
  double fitness;
} Organism;

void init_organism(Organism *o, Genotype *g) {
  o->genotype = g;
  o->activations = calloc(sizeof(double), g->num_nodes);
  o->fitness = 0.0;
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

void init_in_activations(double *activations, Genotype *g) {
  for (int i=0; i<g->num_in; i++) {
    activations[i] = g->nodes[i].value;
  }
}

void print_out_activations(double *activations, Genotype *g) {
  printf("phenotype: ");
  for (int i=g->num_in; i<g->num_in + g->num_out; i++) {
    printf("%4.4f ", activations[i]);
  }
  printf("\n");
}

void print_all_activations(double *activations, int num_nodes) {
  for (int i=0; i<num_nodes; i++) {
    printf("%4.2f ", activations[i]);
  }
  printf("\n");
}

void sa(Organism *o, int timesteps, double decay) {
  Genotype *g = o->genotype;
  double *activations = o->activations;
  bzero(activations, sizeof(activations));
  init_in_activations(activations, g);

  Edgeinfo edge_info[g->num_edges + g->num_nodes + 1];
  init_edge_info(edge_info, g);

  double next_activations[g->num_nodes];
  bzero(next_activations, sizeof(next_activations));

  for (int timestep=1; timestep<=timesteps; timestep++) {
    double acc = 0.0;
    double *activation = activations;
    double *next_activation = next_activations;
    Edgeinfo *to_node = edge_info;
    Node *cur_node = g->nodes;
    for (;;) {
      while (to_node->node_idx == EDGE_INFO_NEXT) {
        *next_activation++ = cur_node->threshold_func(decay * (*activation++ + acc));
        cur_node++;
        to_node++;
        acc = 0.0;
      }
      if (to_node->node_idx == EDGE_INFO_STOP) {
        break;
      }
      assert(to_node->node_idx >= 0 && to_node->node_idx < g->num_nodes);
      acc += activations[to_node->node_idx] * to_node->weight;
      if (verbose > 1)
        printf("n%ld acc=%4.2lf to=%d cur=%4.2lf weight=%4.2lf\n", (cur_node - g->nodes), acc,
            to_node->node_idx, activations[to_node->node_idx], to_node->weight);
      to_node++;
    }
    if (verbose >= 9)
      printf("timestep: %d\n", timestep);
      print_all_activations(next_activations, g->num_nodes);
    memcpy(activations, next_activations, sizeof(next_activations));
  }

  if (verbose)
    print_out_activations(activations, g);
}

// -- world ------------------------------------------------------------------

typedef struct world_t {
  int random_seed;
  int num_organisms;
  int sa_timesteps;
  int generations_per_epoch;
  int num_epochs;
  int num_nodes;
  int num_edges;
  int num_in;
  int num_out;
  double decay_rate;
  Genotype *genotypes;
  Organism *organisms;
  double (*phenotype_fitness_func)(struct world_t *, Organism *);
  int generation;
  double c;
} World;

double phenotype_fitness(World *, Organism *);

World *create_world_full(int random_seed, int num_organisms, int sa_timesteps,
                         int generations_per_epoch, int num_epochs, int num_nodes, int num_edges,
                         int num_in, int num_out, double decay_rate) {
  World *w = calloc(sizeof(World), 1);
  w->random_seed = random_seed;
  w->num_organisms = num_organisms;
  w->sa_timesteps = sa_timesteps;
  w->generations_per_epoch = generations_per_epoch;
  w->num_epochs = num_epochs;
  w->num_nodes = num_nodes;
  w->num_edges = num_edges;
  w->num_in = num_in;
  w->num_out = num_out;
  w->decay_rate = decay_rate;
  w->genotypes = calloc(sizeof(Genotype), num_organisms);
  w->organisms = calloc(sizeof(Organism), num_organisms);
  w->phenotype_fitness_func = phenotype_fitness;
  w->c = 0.0;
  return w;
}

World *create_world(int random_seed, int num_organisms, int sa_timesteps,
                    int generations_per_epoch, int num_epochs, int num_nodes, int num_edges) {
  return create_world_full(random_seed, num_organisms, sa_timesteps,
      generations_per_epoch, num_epochs, num_nodes, num_edges, 2, 2, 1.0);
}

void init_random_population(World *w) {
  for (int n=0; n<w->num_organisms; n++) {
    init_random_genotype(&w->genotypes[n], w->num_edges, w->num_nodes, w->num_in, w->num_out);
    if (verbose > 1)
      print_genotype(&w->genotypes[n]);
    init_organism(&w->organisms[n], &w->genotypes[n]);
  }
}

void set_phenotypes_and_fitnesses(World *w) {
  if (!dot)
    printf("  generation %d\n", w->generation);
  for (int n=0; n<w->num_organisms; n++) {
    Organism *o = &w->organisms[n];
    sa(o, w->sa_timesteps, w->decay_rate);
    if (dot)
      print_organism_dot(o);
    o->fitness = w->phenotype_fitness_func(w, o);
    if (verbose)
      printf("fitness: %lf\n", o->fitness);
  }
}

Organism *mutate(Organism *);

void run_generation(World *w) {
  set_phenotypes_and_fitnesses(w);
  int i = 0; // why?
  Organism *o;
  for (o = w->organisms, i = 0; i < w->num_organisms; o++, i++) {
    Organism *new_o = mutate(o);
    free_organism(o);
    w->organisms[i] = *new_o;
  }
}

void run_world(World *w) {
  srand(w->random_seed);
  init_random_population(w);
  for (int e=0; e<w->num_epochs; e++) {
    if (!dot)
      printf("epoch %d\n", e);
    for (w->generation=0; w->generation<w->generations_per_epoch; w->generation++) {
      run_generation(w);
    }
  }
}

void free_world(World *w) {
  for (int i=0; i<w->num_organisms; i++)
    free_organism(&w->organisms[i]); // will free associated genotype
}

// -- fitness ----------------------------------------------------------------

double many_small_hills(double *phenotype) { // length is 2
  return cos(phenotype[0] * 20.0) * sin(phenotype[1] * 20.0);
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

double phenotype_fitness(World *w, Organism *o) {
  double phenotype[2] = {
    o->activations[o->genotype->num_in],
    o->activations[o->genotype->num_in+1]
  };
  return //many_small_hills(phenotype) +
    (1.0 - distance(w->c, w->c, phenotype[0], phenotype[1]));
}

// -- next generation via crossover and mutation -----------------------------

// select from previous population by fitness (no replacement)
// mutants = 70%
//   weighted choice
//     move-edge   1
//     add-node    1
//     remove-node 1
//     add-edge    1
//     remove-edge 1
//     turn-knob   10

Genotype *copy_genotype(Genotype *g) {
  Genotype *new_g = calloc(sizeof(Genotype), 1);
  new_g->num_nodes = g->num_nodes;
  new_g->num_edges = g->num_edges;
  new_g->num_in = g->num_in;
  new_g->num_out = g->num_out;
  new_g->nodes = calloc(sizeof(Node), new_g->num_nodes);
  memcpy(new_g->nodes, g->nodes, sizeof(Node) * new_g->num_nodes);
  new_g->edges = calloc(sizeof(Edge), new_g->num_edges);
  memcpy(new_g->edges, g->edges, sizeof(Edge) * new_g->num_edges);
  return new_g;
}

Organism *copy_organism(Organism *o) {
  Organism *new_o = calloc(sizeof(Organism), 1);
  new_o->genotype = copy_genotype(o->genotype);
  new_o->activations = calloc(sizeof(double), new_o->genotype->num_nodes);
  new_o->fitness = 0.0;
  return new_o;
}

Organism *add_edge(Organism *o) {
  return copy_organism(o);
}

Organism *move_edge(Organism *o) {
  if (o->genotype->num_edges == 0)
    return add_edge(o);
  Organism *new_o = copy_organism(o);
  Genotype *g = new_o->genotype;
  int selected_edge = rand() % g->num_edges;
  g->edges[selected_edge].src = rand() % g->num_nodes;
  g->edges[selected_edge].dst = rand() % g->num_nodes;
  return new_o;
}

// NEXT do the rest of these

Organism *add_node(Organism *o) {
  return copy_organism(o);
}

Organism *remove_node(Organism *o) {
  return copy_organism(o);
}

Organism *remove_edge(Organism *o) {
  return copy_organism(o);
}

Organism *turn_knob(Organism *o) {
  return copy_organism(o);
}

Organism *mutate(Organism *o) {
  int mutation_type = rand() % 16;
  switch (mutation_type) {
  case 0:
    return move_edge(o);
    break;
  case 1:
    return add_node(o);
    break;
  case 2:
    return remove_node(o);
    break;
  case 3:
    return add_edge(o);
    break;
  case 4:
    return remove_edge(o);
    break;
  default:
    return turn_knob(o);
    break;
  }
}

// crossovers = 30%

// ---------------------------------------------------------------------------

void sa_test() {
  Node nodes[6] = {
    { 0.2, clamp },
    { 0.4, clamp },
    { 0.0, clamp },
    { 0.0, clamp },
    { 0.0, clamp },
    { 0.0, clamp }
  };
  Edge edges[6] = {
    { 0, 4, 1.0 },
    { 1, 4, 1.0 },
    { 4, 2, 1.0 },
    { 4, 5, 1.0 },
    { 5, 0, -1.0 },
    { 5, 3, -1.0 }
  };
  Genotype genotype = { nodes, edges, 6, 6, 2, 2 };
  double activations[6];
  Organism o = { &genotype, activations, 0.0 };
  
  verbose = 9;
  sa(&o, 13, 1.0);
}

void quick_test() {
  verbose = 1;
  World *w = create_world(0, 2, 20, 1, 1, 10, 30);
  run_world(w);
}

void dot_test() {
  dot = true;
  World *w = create_world(0, 1, 20, 1, 1, 10, 30);
  run_world(w);
}

void long_test() {
  World *w = create_world(0, 50, 20, 20, 10, 70, 200);
  run_world(w);
}

int main() {
  quick_test();
  //dot_test();
  //long_test();
  //sa_test();
  return 0;
}
