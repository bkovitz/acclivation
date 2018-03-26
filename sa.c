// todo:
// - clean up options
//   - move any hardcoded ones to world
// - improve organism sanity_check
// - clean up output
//   - option to output dot (given organism/generation/epoch)
//   - option to output dot for every best fitness per generation
//   - option to dump virtual fitness func for given organism/generation/epoch
//   - minimal standard output
// - hill climb (from given organism/generation/epoch)

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>

#define MAXS 128
#define array_len(a) (sizeof(a) / sizeof(a[0]))

int verbose = 0;
bool debug = false;
bool dot = false;

// -- graph ------------------------------------------------------------------

double step(double x) {
  return x > 0.0 ? 1.0 : -1.0;
}

double clamp(double x) {
  return x <= -1.0 ? -1.0 : (x >= 1.0 ? 1.0 : x);
}

// float in range 0 to 1
double rand_float() {
  return (double) rand() / RAND_MAX;
}

// http://www.bearcave.com/misl/misl_tech/wavelets/hurst/random.html
double sample_normal(const double sigma)
{
   double x, y, r2;

   do {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x = -1 + 2 * rand_float();
      y = -1 + 2 * rand_float();
      /* see if it is in the unit circle */
      r2 = x * x + y * y;
   } while (r2 > 1.0 || r2 == 0);
   /* Box-Muller transform */
   return sigma * y * sqrt (-2.0 * log (r2) / r2);
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

bool coin_flip() {
  return rand() & 1;
}

int rand_int(int lb, int ub) {
  return (rand() % (ub - lb + 1)) + lb;
}

double rand_edge_weight() {
  return rand() & 1 ? 1.0 : -1.0;
}

// float in range -1 to 1
double rand_activation() {
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
    //if (n < num_in)
        //g->nodes[n].initial_activation = rand_activation();
        g->nodes[n].initial_activation = rand_int(-100, +100) * 0.01;
    //else
        //g->nodes[n].initial_activation = 0.0;
    g->nodes[n].final_activation = 0.0;
    g->nodes[n].in_use = true;
    g->nodes[n].threshold_func = clamp;
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

#define UNWRITTEN -1000.0

void init_activations(Genotype *g, double *activations) {
  // "in" activations get their initial activations from the Genotype
  for (int i = 0; i < g->num_in; i++) {
    activations[i] = g->nodes[i].initial_activation;
  }
  // all other nodes get initialized to UNWRITTEN
  //for (int o = g->num_in; o < g->num_nodes; o++) {
    //activations[o] = UNWRITTEN;
  //}
  // out nodes get initialized to UNWRITTEN
  for (int o = g->num_in; o < g->num_in + g->num_out; o++) {
    activations[o] = UNWRITTEN;
  }
  // all others get their initial activations from the Genotype
  for (int n = g->num_in + g->num_out; n < g->num_nodes; n++) {
    activations[n] = g->nodes[n].initial_activation;
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
  init_activations(g, activations);

  double activation_deltas[g->num_nodes];

  if (verbose)
    print_all_activations(g, activations);

  for (int timestep = 1; timestep <= timesteps; timestep++) {
    for (int n = 0; n < g->num_nodes; n++) {
      activation_deltas[n] = UNWRITTEN;
    }
    for (int e = 0; e < g->num_edges; e++) {
      Edge *edge = &g->edges[e];
      if (activations[edge->src] != UNWRITTEN) {
        if (activation_deltas[edge->dst] == UNWRITTEN)
          activation_deltas[edge->dst] = 0.0;
        activation_deltas[edge->dst] += edge->weight * activations[edge->src];
      }
    }
    for (int n = 0; n < g->num_nodes; n++) {
      Node *node = &g->nodes[n];
      if (node->in_use) {
        if (activation_deltas[n] != UNWRITTEN) {
          if (activations[n] == UNWRITTEN)
            activations[n] = 0.0;
          activations[n] =
                node->threshold_func(activations[n] + decay * activation_deltas[n]);
        }
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
  double (*phenotype_fitness_func)(struct world_t *, Genotype *);
  int epoch;
  int generation;
  double c1, c2, c3;
  double ridge_radius;
  int mutation_type_ub;
  double extra_mutation_rate;
  double crossover_freq;
  enum { NO_EDGES_ACROSS_PARENTS,
         INHERIT_SRC_EDGES_FROM_MOMMY,
         INHERIT_SRC_EDGES_FROM_BOTH_PARENTS,
         INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS } edge_inheritance;
  int num_candidates;
  double knob_constant;
  enum { KNOB_DISCRETE, KNOB_NORMAL } knob_type;
  bool dump_fitness_nbhd;
  int dump_fitness_epoch;
  int dump_fitness_generation;
} World;

double phenotype_fitness(World *, Genotype *);

World *create_world(int num_organisms) {
  World *w = calloc(1, sizeof(World));
  w->random_seed = 0;
  w->num_organisms = num_organisms;
  w->sa_timesteps = 20;
  w->generations_per_epoch = 20;
  w->num_epochs = 100;
  w->num_nodes = 4;
  w->num_edges = 0;
  w->num_in = 2;
  w->num_out = 2;
  w->decay_rate = 0.01;
  w->genotypes = calloc(num_organisms, sizeof(Genotype)); // THIS IS CRAZY!
  w->organisms = calloc(num_organisms, sizeof(Organism));
  w->phenotype_fitness_func = phenotype_fitness;
  w->epoch = 0;
  w->generation = 0;
  w->c1 = 0.5;
  w->c2 = 2.0;
  w->c3 = 0.45;
  //w->c2 = 1.0;
  //w->c3 = 0.0;
  w->ridge_radius = 0.05;
  w->mutation_type_ub = 16;
  w->extra_mutation_rate = 0.1;
  w->crossover_freq = 0.1;
  w->edge_inheritance = INHERIT_SRC_EDGES_FROM_MOMMY;
  w->num_candidates = 7;
  w->knob_constant = 0.02;
  w->knob_type = KNOB_DISCRETE;
  w->dump_fitness_nbhd = false;
  w->dump_fitness_epoch = -1;
  w->dump_fitness_generation = -1;
  return w;
}

void set_phenotypes_and_fitnesses(World *w) {
  for (int n = 0; n < w->num_organisms; n++) {
    Organism *o = &w->organisms[n];
    sa(o, w->sa_timesteps, w->decay_rate);
//    if (dot)
//      print_organism_dot(o);
    o->fitness = w->phenotype_fitness_func(w, o->genotype);
  }
}

void init_random_population(World *w) {
  for (int n = 0; n < w->num_organisms; n++) {
    init_random_genotype(&w->genotypes[n], w->num_edges, w->num_nodes, w->num_in, w->num_out);
    if (verbose > 1)
      print_genotype(&w->genotypes[n]);
    init_organism(&w->organisms[n], &w->genotypes[n]);
  }
  set_phenotypes_and_fitnesses(w);
}

void mutate(World *w, Organism *);

int tournament_select(World *w) {
  int pool[w->num_candidates];
  for (int n = 0; n < w->num_candidates; n++) {
    pool[n] = rand() % w->num_organisms;
  }
  double max_fitness = -1e20;
  int best_organism_index = -1;
  for (int n = 0; n < w->num_candidates; n++) {
    if (w->organisms[pool[n]].fitness > max_fitness) {
      max_fitness = w->organisms[pool[n]].fitness;
      best_organism_index = pool[n];
    }
  }
  assert(best_organism_index > -1);
  return best_organism_index;
}

void copy_organism(Organism *, Organism *);
void crossover(World *, Organism *, Organism *, Organism *);

void sanity_check(World *w) {
  for (int p = 0; p < w->num_organisms; p++) {
    Organism *o = &w->organisms[p];
    Genotype *g = o->genotype;
    assert(g);
    // check in/out nodes in use
    for (int i = 0; i < g->num_in + g->num_out; i++)
      assert(g->nodes[i].in_use);
    // check edges
    for (int e = 0; e < g->num_edges; e++) {
      Edge *edge = &g->edges[e];
      assert(edge->src >= 0);
      assert(edge->src < g->num_nodes);
      assert(edge->dst >= 0);
      assert(edge->dst < g->num_nodes);
      assert(g->nodes[edge->src].in_use);
      assert(g->nodes[edge->dst].in_use);
    }
  }
}

void dump_fitness_nbhd(World *w);

void run_generation(World *w) {
  Organism *old_population = w->organisms;
  Organism *new_population = calloc(w->num_organisms, sizeof(Organism));
  for (int p = 0; p < w->num_organisms; p++) {
    if (rand_float() <= w->crossover_freq) {
      int mommy = tournament_select(w);
      int daddy = tournament_select(w);
      Organism *baby = &new_population[p];
      baby->genotype = calloc(1, sizeof(Genotype)); // SUPER UGLY
      crossover(w, baby, &old_population[mommy], &old_population[daddy]);
    } else {
      int selected_organism = tournament_select(w);
      copy_organism(&new_population[p], &w->organisms[selected_organism]);
      mutate(w, &new_population[p]);
    }
  }
  for (int p = 0; p < w->num_organisms; p++)
    free_organism(&w->organisms[p]);
  w->organisms = new_population;
  set_phenotypes_and_fitnesses(w);
  if (debug)
    sanity_check(w);
}

int find_best_organism(World *w) {
  double max_fitness = -1e20;
  int max_fitness_index = -1;
  for (int n = 0; n < w->num_organisms; n++) {
    if (w->organisms[n].fitness > max_fitness) {
      max_fitness = w->organisms[n].fitness;
      max_fitness_index = n;
    }
  }
  assert(max_fitness_index > -1);
  return max_fitness_index;
}

void print_best_fitness(World *w) {
  int best_organism_index = find_best_organism(w);
  Organism *o = &w->organisms[best_organism_index];
  double max_fitness = o->fitness;
  Genotype *g = o->genotype;
  printf("    best fitness=%.16lf  index=%d nodes=%d edges=%d g-vector=[%lf %lf] phenotype=[%.16lf %.16lf]\n",
    max_fitness, best_organism_index,
    g->num_nodes_in_use,
    g->num_edges,
    g->nodes[0].initial_activation,
    g->nodes[1].initial_activation,
    g->nodes[2].final_activation,
    g->nodes[3].final_activation);
  print_organism_dot(o);
}

void print_generation_results(World *w) {
  if (!dot) {
    printf("  generation %d\n", w->generation);
    print_best_fitness(w);
  }
}

void dump_fitness_nbhd(World *w) {
  int best_organism_index = find_best_organism(w);
  double m = 2;
  Organism *original = &w->organisms[best_organism_index];
  Organism o;
  copy_organism(&o, original);
  Genotype *g = o.genotype;
  printf("neighborhood:\n");
  for (double dx = -m * w->knob_constant; dx <= m * w->knob_constant; dx += w->knob_constant) {
    for (double dy = -m * w->knob_constant; dy <= m * w->knob_constant; dy += w->knob_constant) {
      g->nodes[0].initial_activation = original->genotype->nodes[0].initial_activation + dx;
      g->nodes[1].initial_activation = original->genotype->nodes[1].initial_activation + dy;
      sa(&o, w->sa_timesteps, w->decay_rate);
      o.fitness = w->phenotype_fitness_func(w, o.genotype);
      printf("  %lf %lf %.16lf %.16lf %lf\n",
        dx,
        dy,
        g->nodes[2].final_activation,
        g->nodes[3].final_activation,
        o.fitness);
    }
  }
}

void free_world(World *w) {
  for (int i = 0; i < w->num_organisms; i++)
    free_organism(&w->organisms[i]); // will free associated genotype
}

void change_fitness_constants(World *w) {
  w->c1 = rand_activation();
}

void run_epoch(World *w, int e) {
  change_fitness_constants(w);
  if (!dot)
    printf("\nepoch %d (c1=%lf, c2=%lf, c3=%lf)\n", e, w->c1, w->c2, w->c3);
  w->epoch = e;
  w->generation = 0;
  set_phenotypes_and_fitnesses(w);
  print_generation_results(w);
  for (w->generation = 1;
       w->generation <= w->generations_per_epoch;
       w->generation++) {
    run_generation(w);
    print_generation_results(w);
//    if (w->dump_fitness_nbhd
//        && w->generation == w->dump_fitness_generation
//        && w->epoch == w->dump_fitness_epoch) {
    if (w->epoch % 10 == 0 && w->generation == w->generations_per_epoch) {
      dump_fitness_nbhd(w);
    }
  }
  fflush(stdout);
}

void dump_virtual_fitness_func(World *w) {
  int best_organism_index = find_best_organism(w);
  Organism o;
  copy_organism(&o, &w->organisms[best_organism_index]);
  double delta = 0.1;
  for (double g1 = -1.0; g1 <= 1.0; g1 += delta) {
    for (double g2 = -1.0; g2 <= 1.0; g2 += delta) {
      o.genotype->nodes[0].initial_activation = g1;
      o.genotype->nodes[1].initial_activation = g2;
      sa(&o, w->sa_timesteps, w->decay_rate);
      o.fitness = w->phenotype_fitness_func(w, o.genotype);
      printf("%lf %lf %lf %lf %lf\n",
        g1,
        g2,
        o.genotype->nodes[2].final_activation,
        o.genotype->nodes[3].final_activation,
        o.fitness);
    }
  }
}

void dump_phenotype_fitness_func(World *w) {
  double delta = 0.05;
  Genotype g;
  init_random_genotype(&g, 0, 4, 2, 2);
  for (double p1 = -1.0; p1 <= 1.0; p1 += delta) {
    for (double p2 = -1.0; p2 <= 1.0; p2 += delta) {
      g.nodes[2].final_activation = p1;
      g.nodes[3].final_activation = p2;
      double fitness = w->phenotype_fitness_func(w, &g);
      printf("%lf %lf %lf\n", p1, p2, fitness);
    }
  }
}

void run_world(World *w) {
  printf("--------------------------------------------------------------------------------\n");
  printf("w->random_seed=%d;\n", w->random_seed);
  printf("w->ridge_radius=%lf;\n", w->ridge_radius);
  printf("w->c2=%lf; w->c3=%lf;\n", w->c2, w->c3);
  printf("w->decay_rate=%lf;\n", w->decay_rate);
  printf("w->mutation_type_ub=%d;\n", w->mutation_type_ub);
  printf("w->extra_mutation_rate=%lf;\n", w->extra_mutation_rate);
  printf("w->crossover_freq=%lf;\n", w->crossover_freq);
  printf("w->edge_inheritance=");
  switch (w->edge_inheritance) {
    case NO_EDGES_ACROSS_PARENTS:
      printf("NO_EDGES_ACROSS_PARENTS");
      break;
    case INHERIT_SRC_EDGES_FROM_MOMMY:
      printf("INHERIT_SRC_EDGES_FROM_MOMMY");
      break;
    case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
      printf("INHERIT_SRC_EDGES_FROM_BOTH_PARENTS");
      break;
    case INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS:
      printf("INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS");
      break;
  }
  puts(";\n");
  printf("w->num_organisms=%d;\n", w->num_organisms);
  printf("w->num_candidates=%d;\n", w->num_candidates);
  printf("w->generations_per_epoch=%d;\n", w->generations_per_epoch);
  printf("w->sa_timesteps=%d;\n", w->sa_timesteps);

  srand(w->random_seed);
  init_random_population(w);
  for (int e = 1; e <= w->num_epochs; e++) {
    run_epoch(w, e);
  }
  dump_virtual_fitness_func(w);
  //free_world();
}

// -- fitness ----------------------------------------------------------------

double many_small_hills(double *phenotype) { // length is 2
  return cos(phenotype[0] * 20.0) * sin(phenotype[1] * 20.0);
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double invv(double target, double radius, double x) {
  double dist = fabs(target - x);
  if (dist >= radius)
    return 0.0;
  else
    return (radius - dist) / radius;
}

double along_ridge(World *w, double x, double y) {
  return invv(0.0, w->ridge_radius, fabs(y - (w->c2 * x + w->c3)));
}

double phenotype_fitness(World *w, Genotype *g) {
  const double sqrt8 = sqrt(8.0);
  double phenotype[2] = {
    g->nodes[g->num_in].final_activation,
    g->nodes[g->num_in + 1].final_activation
  };
  double peak_x = w->c1;
  double peak_y = w->c2 * w->c1 + w->c3;
  if (phenotype[0] != UNWRITTEN && phenotype[1] != UNWRITTEN) {
    return //many_small_hills(phenotype) +
      //(5 * (sqrt8 - distance(w->c1, w->c1, phenotype[0], phenotype[1])));
      along_ridge(w, phenotype[0], phenotype[1]) *
      (5.0 * (sqrt8 - distance(peak_x, peak_y, phenotype[0], phenotype[1])));
  } else {
    return -10.0;
  }
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

bool has_node(Genotype *g, int n) {
  return n < g->num_nodes && g->nodes[n].in_use;
}

bool has_edge(Genotype *g, int src, int dst) {
  for (int e = 0; e < g->num_edges; e++) {
    if (g->edges[e].src == src && g->edges[e].dst == dst) {
      return true;
    }
  }
  return false;
}

void add_edge(Genotype *g, int src, int dst, double weight) {
  assert(has_node(g, src));
  assert(has_node(g, dst));
  g->num_edges++;
  g->edges = realloc(g->edges, sizeof(Edge) * g->num_edges);
  int e = g->num_edges - 1;
  g->edges[e].src = src;
  g->edges[e].dst = dst;
  g->edges[e].weight = weight;
}

void mut_add_edge(Organism *o) {
  Genotype *g = o->genotype;
  add_edge(g, select_in_use_node(g), select_in_use_node(g), rand_edge_weight());
//  g->num_edges++;
//  g->edges = realloc(g->edges, sizeof(Edge) * g->num_edges);
//  int add_index = g->num_edges - 1;
//  g->edges[add_index].src = select_in_use_node(g);
//  g->edges[add_index].dst = select_in_use_node(g);
//  g->edges[add_index].weight = rand_edge_weight();
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
  //g->nodes[add_index].initial_activation = 0.0;
  g->nodes[add_index].initial_activation = rand_activation();
  g->nodes[add_index].final_activation = 0.0;
  g->nodes[add_index].threshold_func = clamp;
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
  assert(selected_node >= g->num_in + g->num_out);
  // mark unused
  g->nodes[selected_node].in_use = false;
  g->num_nodes_in_use--;
  // redirect any affected edges
  for (Edge *e = g->edges; (e - g->edges) < g->num_edges; ) {
    if (e->src == selected_node || e->dst == selected_node)
      // can't increment e because memory may shift
      remove_edge(g, (e - g->edges));
    else
      // can always increment because even when memory doesn't shift, we'll be on last element
      e++;
  }
}

void mut_turn_knob(World *w, Organism *o) {
  int genotype_index = rand() % o->genotype->num_in;
  double nudge;
  switch (w->knob_type) {
    case KNOB_DISCRETE:
      nudge = (rand() & 1) ? w->knob_constant : -w->knob_constant;
      //nudge = rand_activation() * w->knob_constant;
      break;
    case KNOB_NORMAL:
      nudge = sample_normal(w->knob_constant);
      break;
  }
  Node *node_to_change = &o->genotype->nodes[genotype_index];
  node_to_change->initial_activation =
      clamp(node_to_change->initial_activation + nudge);
}

void mutate(World *w, Organism *o) {
  int num_mutations = 1 + (int)(w->extra_mutation_rate * rand_float() * (o->genotype->num_nodes + o->genotype->num_edges));
  for (int i = 0; i < num_mutations; i++) {
    int mutation_type = rand_int(0, w->mutation_type_ub); //rand() % 16;
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
      mut_turn_knob(w, o);
      break;
    }
  }
}

int count_internal_edges(Genotype *g, int start, int end) {
  int num_internal_edges = 0;
  for (int e = 0; e < g->num_edges; e++) {
    Edge *edge = &g->edges[e];
    if (edge->src >= start && edge->src < end &&
        edge->dst >= start && edge->dst < end) {
      num_internal_edges++;
    }
  }
  return num_internal_edges;
}

void crossover(World *w, Organism *b, Organism *m, Organism *d) {
  Genotype *mommy = m->genotype;
  Genotype *daddy = d->genotype;
  Genotype *baby = b->genotype;

  double crossover_frac = rand_float();
  int mommy_crossover_point = mommy->num_nodes * crossover_frac;
  int daddy_crossover_point = daddy->num_nodes * crossover_frac;

  int num_from_mommy = mommy_crossover_point;
  int num_from_daddy = daddy->num_nodes - daddy_crossover_point;
  baby->num_nodes = num_from_mommy + num_from_daddy;
  baby->nodes = malloc(sizeof(Node) * baby->num_nodes);
  assert(baby->num_nodes >= mommy->num_in + mommy->num_out);
  int n = 0;
  for (int m = 0; m < mommy_crossover_point; m++)
    baby->nodes[n++] = mommy->nodes[m];
  for (int d = daddy_crossover_point; d < daddy->num_nodes; d++)
    baby->nodes[n++] = daddy->nodes[d];
  assert(n == baby->num_nodes);

  assert(mommy->num_in == daddy->num_in);
  assert(mommy->num_out == daddy->num_out);
  baby->num_in = mommy->num_in;
  baby->num_out = mommy->num_out;

  int in_use = 0;
  for (int n = 0; n < baby->num_nodes; n++) {
    if (n < baby->num_in + baby->num_out)
      baby->nodes[n].in_use = true;
    if (baby->nodes[n].in_use)
      in_use++;
  }
  baby->num_nodes_in_use = in_use;

//  baby->num_edges = count_internal_edges(mommy, 0, mommy_crossover_point)
//    + count_internal_edges(daddy, daddy_crossover_point, daddy->num_nodes);
  //baby->edges = malloc(sizeof(Edge) * baby->num_edges);
  baby->edges = NULL;
  baby->num_edges = 0;
  //int e = 0;
  for (int m = 0; m < mommy->num_edges; m++) {
    Edge *edge = &mommy->edges[m];
    switch (w->edge_inheritance) {
      case NO_EDGES_ACROSS_PARENTS:
        if (edge->src < mommy_crossover_point && edge->dst < mommy_crossover_point) {
          add_edge(baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_SRC_EDGES_FROM_MOMMY:
      case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
        if (edge->src < mommy_crossover_point &&
            has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS:
        if (coin_flip() &&
            has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(baby, edge->src, edge->dst, edge->weight);
        }
        break;
//    else if ((edge->src < baby->num_nodes && baby->nodes[edge->src].in_use) &&
//               (edge->dst < baby->num_nodes && baby->nodes[edge->dst].in_use) ) { //&&
               //!has_edge(baby, edge->src, edge->dst)) {
      //baby->edges[e++] = *edge;
//      add_edge(baby, edge->src, edge->dst, edge->weight);
    }
  }
  for (int d = 0; d < daddy->num_edges; d++) {
    Edge *edge = &daddy->edges[d];
    int bsrc = (edge->src - daddy_crossover_point) + mommy_crossover_point;
    int bdst = (edge->dst - daddy_crossover_point) + mommy_crossover_point;
    switch (w->edge_inheritance) {
      case NO_EDGES_ACROSS_PARENTS:
        if (edge->src >= daddy_crossover_point && edge->dst >= daddy_crossover_point) {
          //baby->edges[e].src = (edge->src - daddy_crossover_point) + mommy_crossover_point;
          //baby->edges[e].dst = (edge->dst - daddy_crossover_point) + mommy_crossover_point;
          //baby->edges[e].weight = edge->weight;
          //e++;
          add_edge(baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
        if (edge->src >= daddy_crossover_point &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS:
        if (coin_flip() &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(baby, bsrc, bdst, edge->weight);
        }
        break;
    }
  }
  //assert(e == baby->num_edges);
}

// -- measuring acclivation --------------------------------------------------

void set_random_genotype(Organism *o) {
  for (int i = 0; i < o->genotype->num_in; i++)
    o->genotype->nodes[i].initial_activation = rand_activation();
}

double climb_hill(World *w, Organism *o) {
  Genotype *g = o->genotype;
  double last_fitness = -1e10;
  double starting_fitness = o->fitness;

  while (o->fitness > last_fitness) {
    last_fitness = o->fitness;
    mut_turn_knob(w, o);
    sa(o, w->sa_timesteps, w->decay_rate);
    o->fitness = w->phenotype_fitness_func(w, o->genotype);
  }
  return o->fitness - starting_fitness;
}

double get_acclivation(World *w, Organism *test_o) {
  Organism o;
  double total_ending_fitness = 0.0;
  int num_hill_climbers = 3; // w->num_hill_climbers

  for (int i = 0; i < num_hill_climbers; i++) {
    copy_organism(&o, test_o);
    set_random_genotype(&o);
    total_ending_fitness += climb_hill(w, &o);
    free_organism(&o);
  }

  return total_ending_fitness / num_hill_climbers;
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

void dump_phenotype_fitness() {
  World *w = create_world(40);
  dump_phenotype_fitness_func(w);
}

void quick_test(int seed) {
  verbose = 1;
  World *w = create_world(2);
  w->random_seed = seed;
  w->generations_per_epoch = 1;
  w->num_epochs = 1;
  w->num_nodes = 10;
  w->num_edges = 30;
  run_world(w);
}

void dot_test(int seed) {
  dot = true;
  World *w = create_world(1);
  w->random_seed = seed;
  w->generations_per_epoch = 1;
  w->num_epochs = 1;
  w->num_nodes = 10;
  w->num_edges = 30;
  run_world(w);
}

void long_test_start_small(int seed) {
  World *w = create_world(40);
  w->random_seed = seed;
  w->num_epochs = 50;
  w->edge_inheritance = INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS;
  w->c2 = 2.0;
  w->c3 = 0.45;
  w->ridge_radius = 0.2;
  w->crossover_freq = 0.3;
  w->dump_fitness_nbhd = true;
  w->dump_fitness_epoch = 5;
  w->dump_fitness_generation = 20;
  run_world(w);
}

typedef struct {
  int id;
  int random_seed;
  double ridge_radius;
  enum {YX_RIDGE, OBLIQUE_RIDGE} ridge_type;
  int edge_inheritance;
} PARAMS;

void init_params(PARAMS *p, int random_seed) {
  p->id = 0;
  p->random_seed = random_seed;
  p->ridge_radius = 0.05;
  p->ridge_type = OBLIQUE_RIDGE;
  p->edge_inheritance = INHERIT_SRC_EDGES_FROM_MOMMY;
}

void reopen_stdout_from_param(PARAMS *p) {
  char filename[MAXS];

  snprintf(filename, sizeof(filename), "out%d", p->id);
  if (freopen(filename, "w", stdout) == NULL) {
    perror(filename);
    exit(errno);
  }
}

World *create_world_from_param(PARAMS *p) {
  World *w = create_world(40);
  w->random_seed = p->random_seed;
  w->num_epochs = 400;

  w->ridge_radius = p->ridge_radius;
  switch (p->ridge_type) {
    case YX_RIDGE:
      w->c2 = 1.0;
      w->c3 = 0.0;
      break;
    case OBLIQUE_RIDGE:
      w->c2 = 2.0;
      w->c3 = 0.45;
      break;
  }
  w->edge_inheritance = p->edge_inheritance;

  return w;
}

void parameter_sweep(int seed) {
  static double ridge_radii[] = {0.05, 0.2};
  static int ridge_types[] = {YX_RIDGE, OBLIQUE_RIDGE};
  static int edge_inheritances[] = {
    NO_EDGES_ACROSS_PARENTS,
    INHERIT_SRC_EDGES_FROM_MOMMY,
    INHERIT_SRC_EDGES_FROM_BOTH_PARENTS,
    INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS
  };

  PARAMS p;
  init_params(&p, seed);

  for (int rr = 0; rr < array_len(ridge_radii); rr++) {
    for (int rt = 0; rt < array_len(ridge_types); rt++) {
      for (int ei = 0; ei < array_len(edge_inheritances); ei++) {
        p.ridge_radius = ridge_radii[rr];
        p.ridge_type = ridge_types[rt];
        p.edge_inheritance = edge_inheritances[ei];
        World *w = create_world_from_param(&p);
        reopen_stdout_from_param(&p);
        run_world(w);
        fflush(stdout);
        //free_world(w);
        p.id++;
      }
    }
  }
}

void one_long_epoch(int seed) {
  World *w = create_world(40);
  w->random_seed = seed;
  w->generations_per_epoch = 5000;
  w->num_epochs = 1;
  run_world(w);
}

void good_run_oblique() {
  World *w = create_world(40);
  w->random_seed=203540935;
  w->ridge_radius=0.200000;
  w->c2=2.000000; w->c3=0.450000;
  w->decay_rate=0.010000;
  w->mutation_type_ub=15;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;
  w->edge_inheritance=INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS;

  //w->num_organisms=40;
  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;

  w->num_epochs = 400;
  run_world(w);
}

void good_run_oblique2() {
  World *w = create_world(40);
  w->random_seed=203540935;
  w->ridge_radius=0.200000;
  w->c2=2.000000; w->c3=0.450000;
  w->decay_rate=0.010000;
  w->mutation_type_ub=16;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;
  w->edge_inheritance=INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS;

  //w->num_organisms=40;
  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;
  run_world(w);
}

int get_seed(char **argv, int argc) {
  if (argc > 1) {
    return atoi(argv[1]);
  } else {
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    return tm.tv_nsec;
  }
}

void acclivation_test(int seed) {
  World *w = create_world(40);
  w->random_seed = seed;
  w->num_epochs = 20;
  run_world(w);
  double acclivation = get_acclivation(w, &w->organisms[rand() % 40]);
  printf("acclivation = %lf\n", acclivation);
}

int main(int argc, char **argv) {
  int seed = get_seed(argv, argc);
  //sa_test();
  //quick_test(seed);
  //dot_test(seed);
  //long_test(seed);
  //long_test_start_small(seed); // the main test
  //parameter_sweep(seed);
  //good_run_oblique();
  good_run_oblique2();
  //one_long_epoch(seed);
  //dump_virt_test(seed);
  //dump_phenotype_fitness();
  //acclivation_test(seed);
  return 0;
}
