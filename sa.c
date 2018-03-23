// todo:
// + world (contains pop, seed, options)
// + run world
// + generations per epoch
// + fitness function
// + mutation
// ~ crossover [push]
// + tournament_select
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
#include <time.h>

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
  bool in_use;
  double initial_activation;
  double final_activation;
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
  int num_nodes_in_use;
  int num_edges;
  int num_in;
  int num_out;
} Genotype;

double rand_edge_weight() {
  return rand() & 1 ? 1.0 : -1.0;
}

// float in range -1 to 1
double rand_value() {
  return (double) rand() / RAND_MAX * 2.0 - 1.0;
}

// MAYBE: shuffle and take without replace for deterministic select
int select_in_use_node(Genotype *g) {
  for (;;) {
    int index = rand() % g->num_nodes;
    if (g->nodes[index].in_use)
      return index;
  }
}

int select_in_use_removable_node(Genotype *g) {
  for (;;) {
    int unremovable = g->num_in + g->num_out;
    if (g->num_nodes <= unremovable)
      return -1;
    int index = unremovable + (rand() % (g->num_nodes - unremovable));
    if (g->nodes[index].in_use)
      return index;
  }
}

int take_first_unused_node(Genotype *g) {
  for (int n = 0; n < g->num_nodes; n++) {
    if (!g->nodes[n].in_use) {
      g->nodes[n].in_use = true;
      return n;
    }
  }
  return -1;
}

void init_random_genotype(Genotype *g, int num_edges, int num_nodes, int num_in,
        int num_out) {
  assert(num_in >= 1);
  assert(num_out >= 1);
  assert(num_in + num_out <= num_nodes);
  g->nodes = malloc(sizeof(Node) * num_nodes);
  g->edges = malloc(sizeof(Edge) * num_edges);
  g->num_nodes = num_nodes;
  g->num_nodes_in_use = num_nodes;
  g->num_edges = num_edges;
  g->num_in = num_in;
  g->num_out = num_out;
  for (int n = 0; n < num_nodes; n++) {
    if (n < num_in)
        g->nodes[n].initial_activation = rand_value();
    else
        g->nodes[n].initial_activation = 0.0;
    g->nodes[n].final_activation = 0.0;
    g->nodes[n].in_use = true;
    g->nodes[n].threshold_func = (rand() & 1) ? sigmoid : step;
  }
  for (int e = 0; e < num_edges; e++) {
    g->edges[e].src = rand() % num_nodes;
    g->edges[e].dst = rand() % num_nodes;
    g->edges[e].weight = rand_edge_weight();
  }
}

void print_genotype(Genotype *g) {
  for (int e = 0; e < g->num_edges; e++) {
    printf("e%d : %d -> %d\n", e, g->edges[e].src, g->edges[e].dst);
  }
}

void free_genotype(Genotype *g) {
  free(g->nodes);
  free(g->edges);
}

void print_phenotype(Genotype *g) {
  printf("phenotype: ");
  for (int p = g->num_in; p < g->num_in + g->num_out; p++) {
    printf("%4.4f ", g->nodes[p].final_activation);
  }
  printf("\n");
}

// -- organism ---------------------------------------------------------------

typedef struct {
  Genotype *genotype;
  double fitness;
} Organism;

void init_organism(Organism *o, Genotype *g) {
  o->genotype = g;
  o->fitness = 0.0;
}

void free_organism(Organism *o) {
  free_genotype(o->genotype);
}

void print_organism_dot(Organism *o) {
  Genotype *g = o->genotype;

  printf("digraph g {\n");
  printf("  { rank=source edge [style=\"invis\"] ");
  for (int i = 0; i < g->num_in - 1; i++)
    printf("n%d ->", i);
  printf(" n%d }\n", g->num_in - 1);
  printf("  { rank=sink edge [style=\"invis\"] ");
  for (int o = 0; o < g->num_out - 1; o++)
    printf("n%d ->", g->num_in + o);
  printf(" n%d }\n", g->num_in + g->num_out - 1);
  for (int n = 0; n < g->num_nodes; n++) {
    if (g->nodes[n].in_use)
      printf("  n%d [label=%.3lf]\n", n, g->nodes[n].final_activation);
  }
  for (int e = 0; e < g->num_edges; e++) {
    printf("  n%d -> n%d [label=%.3lf];\n", g->edges[e].src, g->edges[e].dst,
      g->edges[e].weight);
  }
  printf("}\n");
}

// -- spreading activation ---------------------------------------------------

void init_in_activations(Genotype *g, double *activations) {
  for (int i = 0; i < g->num_in; i++) {
    assert(g->nodes[i].in_use);
    activations[i] = g->nodes[i].initial_activation;
  }
}

void print_all_activations(Genotype *g, double *activations) {
  for (int n = 0; n < g->num_nodes; n++) {
    if (g->nodes[n].in_use)
      printf("%4.2f ", activations[n]);
  }
  printf("\n");
}

void sa(Organism *o, int timesteps, double decay) {
  Genotype *g = o->genotype;

  double activations[g->num_nodes];
  memset(activations, 0, sizeof(double) * g->num_nodes);
  init_in_activations(g, activations);

  double activation_deltas[g->num_nodes];

  if (verbose)
    print_all_activations(g, activations);

  for (int timestep = 1; timestep <= timesteps; timestep++) {
    memset(activation_deltas, 0, sizeof(activation_deltas));
    for (int e = 0; e < g->num_edges; e++) {
      Edge *edge = &g->edges[e];
      activation_deltas[edge->dst] += edge->weight * activations[edge->src];
    }
    for (int n = 0; n < g->num_nodes; n++) {
      Node *node = &g->nodes[n];
      if (node->in_use) {
        activations[n] =
            node->threshold_func(activations[n] + activation_deltas[n]);
      }
    }
    if (verbose >= 9) {
      printf("timestep: %d\n", timestep);
      print_all_activations(g, activations);
    }
  }

  for (int n = 0; n < g->num_nodes; n++) {
    Node *node = &g->nodes[n];
    if (node->in_use)
      node->final_activation = activations[n];
  }

  if (verbose)
    print_phenotype(g);
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
  int num_candidates;
} World;

double phenotype_fitness(World *, Organism *);

World *create_world_full(int random_seed, int num_organisms, int sa_timesteps,
                         int generations_per_epoch, int num_epochs, int num_nodes, int num_edges,
                         int num_in, int num_out, double decay_rate) {
  World *w = calloc(1, sizeof(World));
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
  w->genotypes = calloc(num_organisms, sizeof(Genotype));
  w->organisms = calloc(num_organisms, sizeof(Organism));
  w->phenotype_fitness_func = phenotype_fitness;
  w->c = 0.0;
  w->num_candidates = 5;
  return w;
}

World *create_world(int random_seed, int num_organisms, int sa_timesteps,
                    int generations_per_epoch, int num_epochs, int num_nodes, int num_edges) {
  return create_world_full(random_seed, num_organisms, sa_timesteps,
      generations_per_epoch, num_epochs, num_nodes, num_edges, 2, 2, 1.0);
}

void init_random_population(World *w) {
  for (int n = 0; n < w->num_organisms; n++) {
    init_random_genotype(&w->genotypes[n], w->num_edges, w->num_nodes, w->num_in, w->num_out);
    if (verbose > 1)
      print_genotype(&w->genotypes[n]);
    init_organism(&w->organisms[n], &w->genotypes[n]);
  }
}

void set_phenotypes_and_fitnesses(World *w) {
  if (!dot)
    printf("  generation %d\t", w->generation);
  for (int n = 0; n < w->num_organisms; n++) {
    Organism *o = &w->organisms[n];
    sa(o, w->sa_timesteps, w->decay_rate);
    if (dot)
      print_organism_dot(o);
    o->fitness = w->phenotype_fitness_func(w, o);
    if (verbose)
      printf("fitness: %lf\n", o->fitness);
  }
}

void mutate(Organism *);

int tournament_select(World *w) {
  int pool[w->num_candidates];
  for (int n = 0; n < w->num_candidates; n++) {
    pool[n] = rand() % w->num_organisms;
  }
  double max_fitness = -1e20;
  int max_fitness_index = -1;
  for (int n = 0; n < w->num_candidates; n++) {
    if (w->organisms[pool[n]].fitness > max_fitness) {
      max_fitness = w->organisms[pool[n]].fitness;
      max_fitness_index = pool[n];
    }
  }
  assert(max_fitness_index > -1);
  return max_fitness_index;
}

void copy_organism(Organism *, Organism *);

void run_generation(World *w) {
  set_phenotypes_and_fitnesses(w);
  Organism *new_population = calloc(w->num_organisms, sizeof(Organism));
  for (int p = 0; p < w->num_organisms; p++) {
    int selected_organism = tournament_select(w);
    copy_organism(&new_population[p], &w->organisms[selected_organism]);
    mutate(&new_population[p]);
  }
  for (int p = 0; p < w->num_organisms; p++)
    free_organism(&w->organisms[p]);
  w->organisms = new_population;
}

void print_best_fitness(World *w) {
  double max_fitness = -1e20;
  int max_fitness_index = -1;
  for (int n = 0; n < w->num_organisms; n++) {
    if (w->organisms[n].fitness > max_fitness) {
      max_fitness = w->organisms[n].fitness;
      max_fitness_index = n;
    }
  }
  assert(max_fitness_index > -1);
  printf("    best fitness: %d %lf\n", max_fitness_index, max_fitness);
}

void run_world(World *w) {
  srand(time(NULL)); //w->random_seed);
  init_random_population(w);
  for (int e=0; e<w->num_epochs; e++) {
    if (!dot)
      printf("epoch %d\n", e);
    for (w->generation = 0; w->generation < w->generations_per_epoch; w->generation++) {
      run_generation(w);
      print_best_fitness(w);
    }
  }
}

void free_world(World *w) {
  for (int i = 0; i < w->num_organisms; i++)
    free_organism(&w->organisms[i]); // will free associated genotype
}

// -- fitness ----------------------------------------------------------------

double many_small_hills(double *phenotype) { // length is 2
  return cos(phenotype[0] * 20.0) * sin(phenotype[1] * 20.0);
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double phenotype_fitness(World *w, Organism *o) {
  double phenotype[2] = {
    o->genotype->nodes[o->genotype->num_in].final_activation,
    o->genotype->nodes[o->genotype->num_in + 1].final_activation
  };
  return many_small_hills(phenotype) +
    (1.0 - distance(w->c, w->c, phenotype[0], phenotype[1]));
}

// -- next generation via crossover and mutation -----------------------------

Genotype *copy_genotype(Genotype *g) {
  Genotype *new_g = calloc(sizeof(Genotype), 1);
  new_g->num_nodes = g->num_nodes;
  new_g->num_nodes_in_use = g->num_nodes_in_use;
  new_g->num_edges = g->num_edges;
  new_g->num_in = g->num_in;
  new_g->num_out = g->num_out;
  new_g->nodes = malloc(new_g->num_nodes * sizeof(Node));
  memcpy(new_g->nodes, g->nodes, sizeof(Node) * new_g->num_nodes);
  new_g->edges = malloc(new_g->num_edges * sizeof(Edge));
  memcpy(new_g->edges, g->edges, sizeof(Edge) * new_g->num_edges);
  return new_g;
}

void copy_organism(Organism *new_o, Organism *old_o) {
  new_o->genotype = copy_genotype(old_o->genotype);
  new_o->fitness = old_o->fitness;
}

void mut_add_edge(Organism *o) {
  Genotype *g = o->genotype;
  g->num_edges++;
  g->edges = realloc(g->edges, sizeof(Edge) * g->num_edges);
  int add_index = g->num_edges-1;
  g->edges[add_index].src = select_in_use_node(g);
  g->edges[add_index].dst = select_in_use_node(g);
  g->edges[add_index].weight = rand_edge_weight();
}

void remove_edge(Genotype *g, int e) {
  if (e < g->num_edges - 1)
    memmove(&g->edges[e], &g->edges[e + 1], sizeof(Edge) * (g->num_edges - e - 1));
  g->num_edges--;
}

void mut_remove_edge(Organism *o) {
  Genotype *g = o->genotype;
  if (g->num_edges > 0) {
    int selected_edge = rand() % g->num_edges;
    remove_edge(g, selected_edge);
  }
}

void mut_move_edge(Organism *o) {
  Genotype *g = o->genotype;
  if (g->num_edges == 0)
    return mut_add_edge(o);
  int selected_edge = rand() % g->num_edges;
  if (rand() & 1)
    g->edges[selected_edge].src = select_in_use_node(g);
  else
    g->edges[selected_edge].dst = select_in_use_node(g);
}

void mut_add_node(Organism *o) {
  Genotype *g = o->genotype;
  int add_index = take_first_unused_node(g);
  g->num_nodes_in_use++;
  if (add_index == -1) {
    g->num_nodes++;
    g->nodes = realloc(g->nodes, sizeof(Node) * g->num_nodes);
    add_index = g->num_nodes - 1;
  }
  g->nodes[add_index].initial_activation = 0.0;
  g->nodes[add_index].final_activation = 0.0;
  g->nodes[add_index].threshold_func = (rand() & 1) ? sigmoid : step;
  g->nodes[add_index].in_use = true;
}

void mut_remove_node(Organism *o) {
  Genotype *g = o->genotype;
  int unremovable = g->num_in + g->num_out;
  if (g->num_nodes_in_use <= unremovable)
    return;
  int selected_node = select_in_use_removable_node(g);
  if (selected_node == -1)
    return;
  // mark unused
  g->nodes[selected_node].in_use = false;
  g->num_nodes_in_use--;
  // redirect any affected edges
  for (Edge *e = g->edges; (e - g->edges) < g->num_edges; e++) {
    if (e->src == selected_node || e->dst == selected_node)
      remove_edge(g, (e - g->edges));
  }
}

void mut_turn_knob(Organism *o) {
  int genotype_index = rand() & o->genotype->num_in;
  double nudge = rand_value() * 0.02;
  Node *node_to_change = &o->genotype->nodes[genotype_index];
  node_to_change->initial_activation =
      clamp(node_to_change->initial_activation + nudge);
}

void mutate(Organism *o) {
  int mutation_type = rand() % 16;
  switch (mutation_type) {
  case 0:
    mut_move_edge(o);
    break;
  case 1:
    mut_add_node(o);
    break;
  case 2:
    mut_remove_node(o);
    break;
  case 3:
    mut_add_edge(o);
    break;
  case 4:
    mut_remove_edge(o);
    break;
  default:
    mut_turn_knob(o);
    break;
  }
}

// ---------------------------------------------------------------------------

void sa_test() {
  Node nodes[6] = {
    { true, 0.2, 0.0, clamp },
    { true, 0.4, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp }
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
  Organism o = { &genotype, 0.0 };
  
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
  World *w = create_world(0, 100, 20, 20, 10, 70, 200);
  run_world(w);
}

int main() {
  //quick_test();
  //dot_test();
  long_test();
  //sa_test();
  return 0;
}
