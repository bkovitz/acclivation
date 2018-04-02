#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>

#define M_PI 3.14159265358979323846

#define MAXS 128
#define array_len(a) (sizeof(a) / sizeof(a[0]))

int verbose = 0;
bool quiet = false;
bool debug = false;
bool dot = false;

double max(double x, double y) {
  if (x >= y)
    return x;
  else
    return y;
}

// ----------------------------------------------------------------------

double make_random_seed() {
  struct timespec tm;
  clock_gettime(CLOCK_REALTIME, &tm);
  return tm.tv_nsec;
}

bool coin_flip() {
  return rand() & 1;
}

int rand_int(int lb, int ub) {
  return (rand() % (ub - lb + 1)) + lb;
}

double rand_double(double lb, double ub) {
  return (double)rand() / RAND_MAX * (ub - lb) + lb;
}

// float in range -1 to 1
double rand_activation() {
  return (double) rand() / RAND_MAX * 2.0 - 1.0;
}

// -- activation types -------------------------------------------------------

typedef enum {
  SUM_INCOMING,
  MULT_INCOMING,
  MIN_INCOMING
} ACTIVATION_TYPE;

const char *activation_type_string(ACTIVATION_TYPE activation_type) {
  switch (activation_type) {
  case SUM_INCOMING:
    return "+";
  case MULT_INCOMING:
    return "*";
  case MIN_INCOMING:
    return "m";
  default:
    assert(false);
  }
}

typedef enum {
  ONLY_SUM_INCOMING,
  SUM_AND_MULT_INCOMING,
  SUM_AND_MIN_INCOMING
} ACTIVATION_TYPES;

const char *activation_types_string(ACTIVATION_TYPES activation_types) {
  switch (activation_types) {
  case ONLY_SUM_INCOMING:
    return "ONLY_SUM_INCOMING";
  case SUM_AND_MULT_INCOMING:
    return "SUM_AND_MULT_INCOMING";
  case SUM_AND_MIN_INCOMING:
    return "SUM_AND_MIN_INCOMING";
  default:
    assert(false);
  }
}

typedef enum {
  EDGE_WEIGHTS_ONLY_PLUS_1,
  EDGE_WEIGHTS_POS_OR_NEG
} EDGE_WEIGHTS;

const char *edge_weights_string(EDGE_WEIGHTS edge_weights) {
  switch (edge_weights) {
    case EDGE_WEIGHTS_ONLY_PLUS_1:
      return "EDGE_WEIGHTS_ONLY_PLUS_1";
    case EDGE_WEIGHTS_POS_OR_NEG:
      return "EDGE_WEIGHTS_POS_OR_NEG";
    default:
      assert(false);
  }
}

typedef enum {
  PASS_THROUGH,  // output = activation
  STEEP_SIGMOID
} OUTPUT_TYPE;


typedef enum {
  ONLY_PASS_THROUGH,
  PASS_THROUGH_AND_STEEP_SIGMOID
} OUTPUT_TYPES;

// -- accumulating data ------------------------------------------------------

typedef struct {
  double *array;
  int len;
} DATA;

DATA *create_data() {
  DATA *data = malloc(sizeof(DATA));
  data->array = NULL;
  data->len = 0;
  return data;
}

void free_data(DATA *data) {
  free(data->array);
  free(data);
}

DATA *copy_data(DATA* data) {
  DATA *result = create_data();
  result->array = malloc(data->len * sizeof(double));
  memmove(result->array, data->array, data->len * sizeof(double));
  result->len = data->len;
  return result;
}

void add_datum(DATA *data, double datum) {
  data->array = realloc(data->array, (data->len + 1) * sizeof(double));
  data->array[data->len++] = datum;
}

void reset_data(DATA *data) {
  if (data->array) {
    free(data->array);
    data->array = NULL;
    data->len = 0;
  }
  assert(data->len == 0);
}

int compar_double(const void *pa, const void *pb) {
  double a = *(double *)pa, b = *(double *)pb;
  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else
    return 0;
}

DATA *make_sorted_data(DATA *data) {
  DATA *result = copy_data(data);
  qsort(result->array, result->len, sizeof(result->array[0]),
    compar_double);
  return result;
}

void print_data(DATA *data) {
  if (data->len == 0) {
    puts("NONE");
  } else {
    printf("%lf", data->array[0]);
    for (int i = 1; i < data->len; i++) {
      printf(" %lf", data->array[i]);
    }
    putchar('\n');
  }
}

double sum_data(DATA *data) {
  double sum = 0.0;
  for (int i = 0; i < data->len; i++) {
    sum += data->array[i];
  }
  return sum;
}

double average_data(DATA *data) {
  if (data->len == 0)
    return 0.0;  // "safe" mean: 0.0 if no data
  else
    return sum_data(data) / data->len;
}

void print_stats(DATA *data) {
  if (data->len == 0) {
    puts("fd_mean=NA fd_median=NA fd_sd=NA fd_min=NA fd_max=NA");
  } else {
    DATA *s = make_sorted_data(data);

    double mean = average_data(s);

    double median;
    if ((s->len & 1) == 0)
      median = (s->array[s->len / 2 - 1] + s->array[s->len / 2]) / 2.0;
    else
      median = s->array[s->len / 2];

    DATA *sq_deviations = create_data();
    for (int i = 0; i < s->len; i++) {
      double deviation = s->array[i] - mean;
      add_datum(sq_deviations, deviation * deviation);
    }
    double sd;
    if (sq_deviations->len == 1)
      sd = sqrt(sq_deviations->array[0]);
    else                                  // Bessel's correction
      sd = sqrt(sum_data(sq_deviations) / (sq_deviations->len - 1));
    
    printf("fd_mean=% .6lf  fd_median=% .6lf  fd_sd=% .6lf  fd_min=% .6lf  fd_max=% .6lf\n",
      mean, median, sd, s->array[0], s->array[s->len - 1]);
  }
}

// -- graph ------------------------------------------------------------------

double step(double x) {
  return x > 0.0 ? 1.0 : -1.0;
}

double two_step(double x) {
  if (x >= 0.2)
    return 1.0;
  else if (x <= -0.2)
    return -1.0;
  else
    return x;
}

double clamp(double x) {
  return x <= -1.0 ? -1.0 : (x >= 1.0 ? 1.0 : x);
}

double clamp2(double x, double lb, double ub) {
  if (x < lb)
    return lb;
  else if (x > ub)
    return ub;
  else
    return x;
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
  double yoffset = ycenter - (yscale / 2.0);
  double denom = (1.0 + exp(slope * (xcenter - x)));
  if (verbose > 1)
    printf("x = %lf; denom = %lf\n", x, denom);
  assert(denom != 0.0);
  return (yscale / denom) + yoffset;
}

/*double steep_sigmoid(double x) {
  //double xcenter = 0.0, ymin = -1.0, ymax = 1.0, slope = 4.0;
  double xcenter = 0.5, ymin = 0.0, ymax = 1.0, slope = 2.0;
  double yscale = ymax - ymin;
  double ycenter = (ymax + ymin) / 2.0;
  double yoffset = ycenter - (yscale / 2.0);
  double denom = (1.0 + exp(slope * (xcenter - x)));
  if (verbose > 1)
    printf("x = %lf; denom = %lf\n", x, denom);
  assert(denom != 0.0);
  return (yscale / denom) + yoffset;
}*/

double steep_sigmoid(double x, double xcenter) {
  //double xcenter = 0.0, ymin = -1.0, ymax = 1.0, slope = 4.0;
  //double ymin = 0.0, ymax = 1.0;
  double ymin = 1.0, ymax = 0;
  double slope = 2.0;
  double yscale = ymax - ymin;
  double ycenter = (ymax + ymin) / 2.0;
  double yoffset = ycenter - (yscale / 2.0);
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
  ACTIVATION_TYPE activation_type; // activation level as a function of inputs
  OUTPUT_TYPE output_type;  // output level as a function of activation level
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

void print_genotype(Genotype *g) {
  for (int e = 0; e < g->num_edges; e++) {
    printf("e%d : %d -> %d\n", e, g->edges[e].src, g->edges[e].dst);
  }
}

void free_genotype(Genotype *g) {
  if (g != NULL) {
    free(g->nodes);
    free(g->edges);
    free(g);
  }
}

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
  bool from_turned_knob;
} Organism;

/*void init_organism(Organism *o, Genotype *g) {
  o->genotype = g;
  o->fitness = 0.0;
}*/

void free_organism(Organism *o) {
  if (o != NULL) {
    free_genotype(o->genotype);
    free(o);
  }
}

void free_organisms(Organism **organisms, int num_organisms) {
  if (organisms != NULL) {
    for (int i = 0; i < num_organisms; i++)
      free_organism(organisms[i]);
    free(organisms);
  }
}

/*void copy_organism(Organism *new_o, Organism *old_o) {
  new_o->genotype = copy_genotype(old_o->genotype);
  new_o->fitness = old_o->fitness;
}*/

Organism *copy_organism(Organism *o) {
  Organism *new_o = calloc(1, sizeof(Organism));
  new_o->genotype = copy_genotype(o->genotype);
  new_o->fitness = o->fitness;
  new_o->from_turned_knob = false;
  return new_o;
}

void print_organism_dot(Organism *o, FILE *f) {
  Genotype *g = o->genotype;

  fprintf(f, "digraph g {\n");
  fprintf(f, "  { rank=source edge [style=\"invis\"] ");
  for (int i = 0; i < g->num_in - 1; i++)
    fprintf(f, "n%d ->", i);
  fprintf(f, " n%d }\n", g->num_in - 1);
  fprintf(f, "  { rank=sink edge [style=\"invis\"] ");
  for (int o = 0; o < g->num_out - 1; o++)
    fprintf(f, "n%d ->", g->num_in + o);
  fprintf(f, " n%d }\n", g->num_in + g->num_out - 1);
  for (int n = 0; n < g->num_nodes; n++) {
    if (g->nodes[n].in_use) {
      switch (g->nodes[n].output_type) {
      case PASS_THROUGH:
        fprintf(f, "  n%d [label=\"%.3lf %s P\"]\n", n, g->nodes[n].final_activation,
               activation_type_string(g->nodes[n].activation_type));
        break;
      case STEEP_SIGMOID:
        fprintf(f, "  n%d [label=\"%.3lf %s S %.3lf\"]\n", n, g->nodes[n].final_activation,
               activation_type_string(g->nodes[n].activation_type),
               steep_sigmoid(g->nodes[n].final_activation,
                 g->nodes[n].initial_activation));
        break;
      }
    }
  }
  for (int e = 0; e < g->num_edges; e++) {
    fprintf(f, "  n%d -> n%d [label=%.3lf];\n", g->edges[e].src, g->edges[e].dst,
      g->edges[e].weight);
  }
  fprintf(f, "}\n");
}

// -- ancestor logging -------------------------------------------------------

typedef struct {
  bool enabled;
  char *path;
  FILE *f;
} ANCESTOR_LOG;

ANCESTOR_LOG *create_ancestor_log() {
  ANCESTOR_LOG *log = malloc(sizeof(ANCESTOR_LOG));
  log->enabled = true;
  log->path = "ancestors";
  log->f = NULL;
  return log;
}

void free_ancestor_log(ANCESTOR_LOG *log) {
  free(log);
}

void open_ancestor_log(ANCESTOR_LOG *log) {
  if (log->enabled) {
    assert(log->path);
    log->f = fopen(log->path, "w");
  }
}

void close_ancestor_log(ANCESTOR_LOG *log) {
  if (log->enabled) {
    assert(log->f);
    fclose(log->f);
  }
}

typedef enum {
  LINE,
  CIRCLE
} RIDGE_TYPE;

// -- world ------------------------------------------------------------------

typedef struct world_t {
  int seed;
  int num_organisms;
  int sa_timesteps;
  int generations_per_epoch;
  int num_epochs;
  int num_nodes;
  int num_edges;
  int num_in;
  int num_out;
  double decay;
  double spreading_rate;
  ACTIVATION_TYPES activation_types;
  EDGE_WEIGHTS edge_weights;
  OUTPUT_TYPES output_types;
  bool multi_edges;
  bool allow_move_edge;
  //Genotype *genotypes;
  Organism **organisms;
  double (*phenotype_fitness_func)(struct world_t *, Genotype *);
  double distance_weight;
  bool bumps;
  int epoch;
  int generation;
  double c1, c2, c3;
  double c1_lb, c1_ub;
  double peak_x, peak_y;
  enum { JUMPY_PEAK_MOVEMENT, GRADUAL_PEAK_MOVEMENT } peak_movement;
  double max_dist;
  RIDGE_TYPE ridge_type;
  double ridge_radius;
  int mutation_type_ub;
  double extra_mutation_rate;
  double crossover_freq;
  enum { NO_EDGES_ACROSS_PARENTS,
         INHERIT_SRC_EDGES_FROM_MOMMY,
         INHERIT_SRC_EDGES_FROM_BOTH_PARENTS,
         INHERIT_HALF_OF_CROSSOVER_EDGES,
         INHERIT_HALF_OF_ALL_EDGES,
         INHERIT_ALL_EDGES } edge_inheritance;
  int num_candidates;
  double knob_constant;
  enum { KNOB_DISCRETE, KNOB_NORMAL } knob_type;
  bool dump_fitness_nbhd;
  int dump_fitness_epoch;
  int dump_fitness_generation;
  DATA *epoch_fitness_deltas;
  ANCESTOR_LOG *log;
  int num_hill_climbers;
  int num_generations_measured;
  int num_fitness_increases_from_knob_turn;
} World;

double phenotype_fitness(World *, Genotype *);

//World *create_world(int num_organisms) {
World *create_world() {
  World *w = calloc(1, sizeof(World));
  w->seed = make_random_seed();
  w->num_organisms = 80;
  w->sa_timesteps = 20;
  w->generations_per_epoch = 20;
  w->num_epochs = 100;
  w->num_nodes = 4; //8;
  w->num_edges = 0; //16;
  w->num_in = 2;
  w->num_out = 2;
  w->decay = 0.8;
  w->spreading_rate = 0.2;
  w->activation_types = SUM_AND_MULT_INCOMING; //SUM_AND_MIN_INCOMING;
  w->edge_weights = EDGE_WEIGHTS_ONLY_PLUS_1; //EDGE_WEIGHTS_POS_OR_NEG;
  w->multi_edges = true;
  w->allow_move_edge = false;
  w->output_types = ONLY_PASS_THROUGH;
  //w->genotypes = NULL; //calloc(num_organisms, sizeof(Genotype)); // THIS IS CRAZY!
  w->organisms = NULL; //calloc(num_organisms, sizeof(Organism));
  w->phenotype_fitness_func = phenotype_fitness;
  w->distance_weight = 10.0;
  w->bumps = true;
  w->epoch = 0;
  w->generation = 0;
  w->c1 = 0.5;
  w->c2 = 1.0;
  w->c3 = 0.0;
  w->c1_lb = 0.2;
  w->c1_ub = 0.8;
  w->peak_movement = JUMPY_PEAK_MOVEMENT;
  w->peak_x = 0.0; // dependent variable
  w->peak_y = 0.0; // dependent variable
  w->max_dist = 0.0; // dependent variable
  //w->c2 = 1.0;
  //w->c3 = 0.0;
  w->ridge_type = LINE;
  w->ridge_radius = 0.2;
  w->mutation_type_ub = 10;
  w->extra_mutation_rate = 0.0; //0.1;
  w->crossover_freq = 0.02;
  w->edge_inheritance = INHERIT_SRC_EDGES_FROM_MOMMY;
  w->num_candidates = 7;
  w->knob_constant = 0.02;
  w->knob_type = KNOB_DISCRETE;
  w->dump_fitness_nbhd = false;
  w->dump_fitness_epoch = -1;
  w->dump_fitness_generation = -1;

  w->epoch_fitness_deltas = create_data();
  w->log = create_ancestor_log();
  w->num_hill_climbers = 30;

  w->num_generations_measured = 0;
  w->num_fitness_increases_from_knob_turn = 0;
  return w;
}

double rand_edge_weight(World *w) {
  switch (w->edge_weights) {
    case EDGE_WEIGHTS_POS_OR_NEG:
      return coin_flip() ? 1.0 : -1.0;
    case EDGE_WEIGHTS_ONLY_PLUS_1:
      return 1.0;
    default:
      assert(false);
  }
}

/*void init_random_genotype(World *w, Genotype *g, int num_edges, int num_nodes,
        int num_in, int num_out) {
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
    switch (w->activation_types) {
      case ONLY_SUM_INCOMING:
        g->nodes[n].activation_type = SUM_INCOMING;
        break;
      case SUM_AND_MULT_INCOMING:
        g->nodes[n].activation_type = coin_flip() ? SUM_INCOMING : MULT_INCOMING;
    }
  }
  for (int e = 0; e < num_edges; e++) {
    g->edges[e].src = rand() % num_nodes;
    g->edges[e].dst = rand() % num_nodes;
    g->edges[e].weight = rand_edge_weight(w);
  }
}*/

void init_random_node(World *w, Node *n) {
  //if (n < num_in)
      //g->nodes[n].initial_activation = rand_activation();
      n->initial_activation = rand_int(-100, +100) * 0.01;
  //else
      //g->nodes[n].initial_activation = 0.0;
  n->final_activation = 0.0;
  n->in_use = true;
  n->threshold_func = clamp;
  switch (w->activation_types) {
    case ONLY_SUM_INCOMING:
      n->activation_type = SUM_INCOMING;
      break;
    case SUM_AND_MULT_INCOMING:
      n->activation_type = coin_flip() ? SUM_INCOMING : MULT_INCOMING;
      break;
    case SUM_AND_MIN_INCOMING:
      n->activation_type = coin_flip() ? SUM_INCOMING : MIN_INCOMING;
      break;
  }
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      n->output_type = PASS_THROUGH;
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      n->output_type = coin_flip() ? PASS_THROUGH : STEEP_SIGMOID;
      break;
  }
}

void add_edge(World *, Genotype *, int, int, double);

Genotype *create_random_genotype(World *w) {
  Genotype *g = calloc(1, sizeof(Genotype));
  g->nodes = malloc(sizeof(Node) * w->num_nodes);
  //g->edges = malloc(sizeof(Edge) * w->num_edges);
  g->edges = NULL;
  g->num_nodes = w->num_nodes;
  g->num_nodes_in_use = w->num_nodes;
  g->num_edges = 0; //w->num_edges;
  g->num_in = w->num_in;
  g->num_out = w->num_out;
  for (int n = 0; n < g->num_nodes; n++) {
    init_random_node(w, &g->nodes[n]);
    /* //if (n < num_in)
        //g->nodes[n].initial_activation = rand_activation();
        g->nodes[n].initial_activation = rand_int(-100, +100) * 0.01;
    //else
        //g->nodes[n].initial_activation = 0.0;
    g->nodes[n].final_activation = 0.0;
    g->nodes[n].in_use = true;
    g->nodes[n].threshold_func = clamp;
    switch (w->activation_types) {
      case ONLY_SUM_INCOMING:
        g->nodes[n].activation_type = SUM_INCOMING;
        break;
      case SUM_AND_MULT_INCOMING:
        g->nodes[n].activation_type = coin_flip() ? SUM_INCOMING : MULT_INCOMING;
    }*/
  }
  for (int e = 0; e < g->num_edges; e++) {
    int src = rand() % g->num_nodes;
    int dst = rand() % g->num_nodes;
    double weight = rand_edge_weight(w);
    add_edge(w, g, src, dst, weight);
//    g->edges[e].src = rand() % g->num_nodes;
//    g->edges[e].dst = rand() % g->num_nodes;
//    g->edges[e].weight = rand_edge_weight(w);
  }
  return g;
}

Organism *create_random_organism(World *w) {
  assert(w->num_in >= 1);
  assert(w->num_out >= 1);
  assert(w->num_in + w->num_out <= w->num_nodes);

  Organism *o = calloc(1, sizeof(Organism));
  o->genotype = create_random_genotype(w);
  o->fitness = 0.0;
  o->from_turned_knob = false;
  return o;
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
      printf("%.16f ", activations[n]);
  }
  printf("\n");
}

double src_output(Node *src, int index, double *activations) {
  switch (src->output_type) {
    case PASS_THROUGH:
      return activations[index];
    case STEEP_SIGMOID:
      //return steep_sigmoid(activations[index], src->initial_activation);
      return steep_sigmoid(activations[index], activations[index]);
    default:
      assert(false);
  }
}

void sa(Organism *o, int timesteps, double decay, double spreading_rate) {
  Genotype *g = o->genotype;

  double activations[g->num_nodes];
  memset(activations, 0, sizeof(double) * g->num_nodes);
  init_activations(g, activations);

  double incoming_activations[g->num_nodes];

  if (verbose)
    print_all_activations(g, activations);

  for (int timestep = 1; timestep <= timesteps; timestep++) {
    for (int n = 0; n < g->num_nodes; n++) {
      incoming_activations[n] = UNWRITTEN;
    }
    for (int e = 0; e < g->num_edges; e++) {
      Edge *edge = &g->edges[e];
      double incoming_output = src_output(&g->nodes[edge->src], edge->src,
          activations);
      if (activations[edge->src] != UNWRITTEN) {
        switch (g->nodes[edge->dst].activation_type) {
          case SUM_INCOMING:
            if (incoming_activations[edge->dst] == UNWRITTEN)
              incoming_activations[edge->dst] = 0.0;
            incoming_activations[edge->dst] +=
                  //edge->weight * activations[edge->src];
                  edge->weight * incoming_output; //src_output(src_node);
            break;
          case MULT_INCOMING:
            if (incoming_activations[edge->dst] == UNWRITTEN)
              incoming_activations[edge->dst] = 1.0;
            incoming_activations[edge->dst] *=
                  //edge->weight * activations[edge->src];
                  edge->weight * incoming_output;
            break;
          case MIN_INCOMING:
//            if (incoming_activations[edge->dst] == UNWRITTEN)
//              incoming_activations[edge->dst] = activations[edge->src];
//            else if (activations[edge->src] < incoming_activations[edge->dst])
//              incoming_activations[edge->dst] = activations[edge->src];
            if (incoming_activations[edge->dst] == UNWRITTEN)
              incoming_activations[edge->dst] = incoming_output;
            else if (incoming_output < incoming_activations[edge->dst])
              incoming_activations[edge->dst] = incoming_output;
            break;
        }
      }
    }
    for (int n = 0; n < g->num_nodes; n++) {
      Node *node = &g->nodes[n];
      if (node->in_use) {
        if (incoming_activations[n] != UNWRITTEN) {
          switch (node->activation_type) {
            case SUM_INCOMING:
              if (activations[n] == UNWRITTEN)
                activations[n] = 0.0;
              activations[n] = node->threshold_func(
                  activations[n] + spreading_rate
                                   * incoming_activations[n]
                                   * pow(decay, timestep - 1));
              break;
            case MULT_INCOMING:
              // ignore previous activation
              activations[n] = sigmoid(incoming_activations[n]);
              break;
            case MIN_INCOMING:
              // ignore previous activation
              activations[n] = incoming_activations[n];
          }
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

// ------------------------------------------------------------------------

void set_phenotypes_and_fitnesses(World *w) {
  for (int n = 0; n < w->num_organisms; n++) {
    Organism *o = w->organisms[n];
    sa(o, w->sa_timesteps, w->decay, w->spreading_rate);
//    if (dot)
//      print_organism_dot(o, stdout);
    o->fitness = w->phenotype_fitness_func(w, o->genotype);
  }
}

void init_random_population(World *w) {
  /*if (w->organisms != NULL) {
    for (int i = 0; i < w->num_organisms; i++)
      free_organism(w->organisms[i]);
  }*/
  free_organisms(w->organisms, w->num_organisms);
  w->organisms = calloc(w->num_organisms, sizeof(Organism *));

  for (int n = 0; n < w->num_organisms; n++) {
    w->organisms[n] = create_random_organism(w);
    //init_random_genotype(w, &w->genotypes[n],
        //w->num_edges, w->num_nodes, w->num_in, w->num_out);
    //if (verbose > 1)
      //print_genotype(&w->genotypes[n]);
    //init_organism(&w->organisms[n], &w->genotypes[n]);
  }
  set_phenotypes_and_fitnesses(w);
}

Organism *mutate(World *w, Organism *);

int tournament_select(World *w) {
  int pool[w->num_candidates];
  for (int n = 0; n < w->num_candidates; n++) {
    pool[n] = rand() % w->num_organisms;
  }
  double max_fitness = -1e20;
  int best_organism_index = -1;
  for (int n = 0; n < w->num_candidates; n++) {
    if (w->organisms[pool[n]]->fitness > max_fitness) {
      max_fitness = w->organisms[pool[n]]->fitness;
      best_organism_index = pool[n];
    }
  }
  assert(best_organism_index > -1);
  return best_organism_index;
}

Organism *copy_organism(Organism *);
Organism *crossover(World *, Organism *, Organism *);

void sanity_check(World *w) {
  for (int p = 0; p < w->num_organisms; p++) {
    Organism *o = w->organisms[p];
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

void log_organisms(World *w) {
  if (w->log->enabled) {
    for (int i = 0; i < w->num_organisms; i++) {
      Organism *o = w->organisms[i];
      Genotype *g = o->genotype;
      fprintf(w->log->f, "organism [%d,%d,%d] fitness=%20.16lf  nodes=%2d  edges=%2d  g-vector=[% lf % lf] phenotype=[% .16lf % .16lf]\n",
        w->epoch,
        w->generation,
        i,
        o->fitness,
        g->num_nodes_in_use,
        g->num_edges,
        g->nodes[0].initial_activation,
        g->nodes[1].initial_activation,
        g->nodes[2].final_activation,
        g->nodes[3].final_activation);
      print_organism_dot(o, w->log->f);
    }
  }
}

int prev_generation(World *w) {
  return w->generation > 1 ? (w->generation - 1) : (w->generations_per_epoch - 1);
}

int maybe_prev_epoch(World *w) {
  return w->generation > 1 ? w->epoch : (w->epoch - 1);
}

void log_mutation_start(World *w, int parent, int child) {
  if (w->log->enabled) {
    fprintf(w->log->f, "mutation %d,%d,%d %d,%d,%d [ ",
      maybe_prev_epoch(w),
      prev_generation(w),
      parent,
      w->epoch,
      w->generation,
      child);
  }
}

void log_mutation(World *w, char *type) {
  if (w->log->enabled) {
    fprintf(w->log->f, "%s ", type);
  }
}

void log_mutation_end(World *w) {
  if (w->log->enabled) {
    fprintf(w->log->f, "]\n");
  }
}

void log_crossover(World *w, int mommy, int daddy, int child) {
  if (w->log->enabled) {
    fprintf(w->log->f, "crossover %d,%d,%d %d,%d,%d %d,%d,%d\n",
      maybe_prev_epoch(w),
      prev_generation(w),
      mommy,
      maybe_prev_epoch(w),
      prev_generation(w),
      daddy,
      w->epoch,
      w->generation,
      child);
  }
}

void dump_fitness_nbhd(World *w);

void run_generation(World *w) {
  Organism **old_population = w->organisms;
  Organism **new_population = calloc(w->num_organisms, sizeof(Organism *));
  for (int p = 0; p < w->num_organisms; p++) {
    if (rand_float() <= w->crossover_freq) {
      int mommy = tournament_select(w);
      int daddy = tournament_select(w);
      //Organism *baby = &new_population[p];
      //baby->genotype = calloc(1, sizeof(Genotype)); // SUPER UGLY
      //crossover(w, baby, &old_population[mommy], &old_population[daddy]);
      new_population[p] = crossover(w, old_population[mommy], old_population[daddy]);
      log_crossover(w, mommy, daddy, p);
    } else {
      int selected_organism = tournament_select(w);
      //copy_organism(&new_population[p], &w->organisms[selected_organism]);
      log_mutation_start(w, selected_organism, p);
      //mutate(w, &new_population[p]);
      new_population[p] = mutate(w, old_population[selected_organism]);
      log_mutation_end(w);
    }
  }
  free_organisms(w->organisms, w->num_organisms);
  /*for (int p = 0; p < w->num_organisms; p++)
    free_organism(w->organisms[p]);
  free(old_population);*/
  w->organisms = new_population;
  set_phenotypes_and_fitnesses(w);
  if (debug)
    sanity_check(w);
  if (!quiet)
    log_organisms(w);
}

int find_best_organism_in(Organism **organisms, int num_organisms) {
  double max_fitness = -1e20;
  int max_fitness_index = -1;
  for (int n = 0; n < num_organisms; n++) {
    if (organisms[n]->fitness > max_fitness) {
      max_fitness = organisms[n]->fitness;
      max_fitness_index = n;
    }
  }
  assert(max_fitness_index > -1);
  return max_fitness_index;
}

int find_best_organism(World *w) {
  return find_best_organism_in(w->organisms, w->num_organisms);
}
/*int find_best_organism(World *w) {
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
}*/

double find_best_fitness(World *w) {
  int best_organism_index = find_best_organism(w);
  return w->organisms[best_organism_index]->fitness;
}

void print_best_fitness(World *w) {
  int best_organism_index = find_best_organism(w);
  Organism *o = w->organisms[best_organism_index];
  double max_fitness = o->fitness;
  Genotype *g = o->genotype;
  printf("    best fitness=%20.16lf  index=%2d  nodes=%2d  edges=%2d  g-vector=[% lf % lf] phenotype=[% .16lf % .16lf]\n",
    max_fitness, best_organism_index,
    g->num_nodes_in_use,
    g->num_edges,
    g->nodes[0].initial_activation,
    g->nodes[1].initial_activation,
    g->nodes[2].final_activation,
    g->nodes[3].final_activation);
  print_organism_dot(o, stdout);
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
  Organism *original = w->organisms[best_organism_index];
  //Organism o;
  //copy_organism(&o, original);
  Organism *o = copy_organism(original);
  Genotype *g = o->genotype;
  printf("neighborhood:\n");
  for (double dx = -m * w->knob_constant; dx <= m * w->knob_constant; dx += w->knob_constant) {
    for (double dy = -m * w->knob_constant; dy <= m * w->knob_constant; dy += w->knob_constant) {
      g->nodes[0].initial_activation = original->genotype->nodes[0].initial_activation + dx;
      g->nodes[1].initial_activation = original->genotype->nodes[1].initial_activation + dy;
      sa(o, w->sa_timesteps, w->decay, w->spreading_rate);
      o->fitness = w->phenotype_fitness_func(w, o->genotype);
      printf("  % lf % lf % .16lf % .16lf % lf\n",
        dx,
        dy,
        g->nodes[2].final_activation,
        g->nodes[3].final_activation,
        o->fitness);
    }
  }
}

// -- measuring acclivity ----------------------------------------------------

void set_random_gvector(Organism *o) {
  for (int i = 0; i < o->genotype->num_in; i++)
    o->genotype->nodes[i].initial_activation = rand_activation();
}

void nudge_candidate(Organism *o, int node_index, double nudge_amount) {
  Node *node_to_change = &o->genotype->nodes[node_index];
  node_to_change->initial_activation =
    clamp(node_to_change->initial_activation + nudge_amount);
}

typedef struct {
  double fitness_delta;
  double ending_fitness;
} HILL_CLIMBING_RESULT;

HILL_CLIMBING_RESULT climb_hill(World *w, Organism *o) {
  sa(o, w->sa_timesteps, w->decay, w->spreading_rate);
  o->fitness = w->phenotype_fitness_func(w, o->genotype);

  const double starting_fitness = o->fitness;
  double last_fitness = starting_fitness;
  double nudge_amount = w->knob_constant;

  int num_candidates = w->num_in << 1;
  Organism *candidates[num_candidates];

  int num_neutral_steps = 0;

  //printf("starting_fitness = %lf\n", starting_fitness);
  for (;;) {
    // create organisms that take a step in each possible direction
    for (int i = 0; i < num_candidates; i++) {
      candidates[i] = copy_organism(o);
      Organism *candidate = candidates[i];
      nudge_candidate(candidate, i >> 1, (i & 1) ? nudge_amount : -nudge_amount);
      sa(candidate, w->sa_timesteps, w->decay, w->spreading_rate);
      candidate->fitness = w->phenotype_fitness_func(w, candidate->genotype);
      //printf("candidate %d fitness = %lf\n", i, candidate->fitness);
      //print_phenotype(candidate->genotype);
    }
    // take the best positive step, or bail out
    int best_candidate_index = find_best_organism_in(candidates, num_candidates);
    if (candidates[best_candidate_index]->fitness > last_fitness) {
      free_organism(o);
      o = copy_organism(candidates[best_candidate_index]);
      last_fitness = o->fitness;
      //printf("fitness, last fitness-> %lf, %lf\n", o->fitness, last_fitness);
      for (int i = 0; i < num_candidates; i++)
        free_organism(candidates[i]);
      num_neutral_steps = 0;
    } else if (candidates[best_candidate_index]->fitness == last_fitness) {
      // neutral plateau
      if (++num_neutral_steps >= 200) {
        //printf("crazy plateau\n");
        break;
      }
      free_organism(o);
      o = copy_organism(candidates[rand_int(0, num_candidates - 1)]);
      for (int i = 0; i < num_candidates; i++)
        free_organism(candidates[i]);
    } else {
      // reached a peak
      break;
    }
  }
  HILL_CLIMBING_RESULT result = { last_fitness - starting_fitness, last_fitness };
  free_organism(o);
  return result;
}

HILL_CLIMBING_RESULT measure_acclivity(World *w, Organism *test_o) {
  Organism *o;
  HILL_CLIMBING_RESULT total = { 0.0, 0.0 };

  srand(0);
  for (int i = 0; i < w->num_hill_climbers; i++) {
    o = copy_organism(test_o);
    set_random_gvector(o);
    //printf("gvector-> %lf, %lf\n", o->genotype->nodes[0].initial_activation, o->genotype->nodes[1].initial_activation);
    HILL_CLIMBING_RESULT result = climb_hill(w, o);
    total.fitness_delta += result.fitness_delta;
    total.ending_fitness += result.ending_fitness;
    //printf("->delta %lf, abs %lf\n", result.fitness_delta, result.ending_fitness);
    //free_organism(o);
  }

  total.fitness_delta /= w->num_hill_climbers;
  total.ending_fitness /= w->num_hill_climbers;
  return total;
}

HILL_CLIMBING_RESULT phenotype_acclivity(World *w) {
  Node nodes[] = {
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp }
  };
  Edge edges[] = {
    { 0, 2, 1.0 },
    { 1, 3, 1.0 }
  };
  Genotype genotype = { nodes, edges, 4, 4, 2, 2, 2 };
  Organism null_organism = { &genotype, 0.0 };
  
  return measure_acclivity(w, &null_organism);
}

// ----------------------------------------------------------------------

void free_world(World *w) {
  free_organisms(w->organisms, w->num_organisms);
  //for (int i = 0; i < w->num_organisms; i++)
    //free_organism(w->organisms[i]); // will free associated genotype
  free_data(w->epoch_fitness_deltas);
  free_ancestor_log(w->log);
  free(w);
}

double distance(double x1, double y1, double x2, double y2);

void change_fitness_constants(World *w) {
  //w->c1 = rand_activation();
  const double sqrt2 = sqrt(2.0);
  switch (w->ridge_type) {
  case LINE:
    switch (w->peak_movement) {
    case JUMPY_PEAK_MOVEMENT:
      w->c1 = rand_double(w->c1_lb, w->c1_ub);
      break;
    case GRADUAL_PEAK_MOVEMENT:
      w->c1 += sample_normal(0.4);
      if (w->c1 < w->c1_lb || w->c1 > w->c1_ub)
        w->c1 = rand_double(w->c1_lb, w->c1_ub);
      break;
    default:
      assert(false);
    }
    // dependent variables
    w->peak_x = w->c1;
    w->peak_y = clamp2(w->c2 * w->c1 + w->c3, w->c1_lb, w->c1_ub);
    //w->max_dist = sqrt2; // wrong
    w->max_dist = distance(w->peak_x, w->peak_y, -1.0, -1.0);
    double d = distance(w->peak_x, w->peak_y, -1.0, +1.0);
    if (d > w->max_dist)
      w->max_dist = d;
    d = distance(w->peak_x, w->peak_y, +1.0, -1.0);
    if (d > w->max_dist)
      w->max_dist = d;
    d = distance(w->peak_x, w->peak_y, +1.0, +1.0);
    if (d > w->max_dist)
      w->max_dist = d;
    printf("peak=(%lf,%lf) max_dist=%lf\n", w->peak_x, w->peak_y, w->max_dist);
    break;
  case CIRCLE:
    switch (w->peak_movement) {
    case JUMPY_PEAK_MOVEMENT:
      w->c1 = rand_double(0.0, 2 * M_PI);
      break;
    case GRADUAL_PEAK_MOVEMENT:
      {
        //double delta = sample_normal(2 * M_PI / 10.0);
        double delta = rand_double(2 * M_PI / 10, 2 * M_PI / 5.0);
        printf("c1 delta = %lf\n", delta);
        w->c1 += delta;
      }
      break;
    default:
      assert(false);
    }
    // dependent variables
    w->peak_x = 0.5 * cos(w->c1);
    w->peak_y = 0.5 * sin(w->c1);
    w->max_dist = 0.5 + sqrt2; // radius + edge-of-circle-to-corner
    break;
  }
}

void check_knob_turn(World *w, double last_fitness) {
  int best_organism_index = find_best_organism(w);
  Organism *best = w->organisms[best_organism_index];
  double fitness_delta = best->fitness - last_fitness;
  w->num_generations_measured++;
  if (fitness_delta > 0.0) {
    if (best->from_turned_knob) {
      w->num_fitness_increases_from_knob_turn++;
      //printf("    from_knob_turn: delta: %lf\n", fitness_delta);
    }
  }
}

void run_epoch(World *w, int e) {
  double epoch_start_fitness;
  double last_fitness;
  change_fitness_constants(w);
  if (!quiet && !dot)
    printf("\nepoch %d (c1=%lf, c2=%lf, c3=%lf)\n", e, w->c1, w->c2, w->c3);
  w->epoch = e;
  w->generation = 0;
  set_phenotypes_and_fitnesses(w);
  if (!quiet)
    print_generation_results(w);
  epoch_start_fitness = max(find_best_fitness(w), 0.0);
  for (w->generation = 1;
       w->generation <= w->generations_per_epoch;
       w->generation++) {
    last_fitness = max(find_best_fitness(w), 0.0);
    run_generation(w);
    if (!quiet) {
      print_generation_results(w);
//    if (w->dump_fitness_nbhd
//        && w->generation == w->dump_fitness_generation
//        && w->epoch == w->dump_fitness_epoch)
      if (w->epoch % 10 == 0 && w->generation == w->generations_per_epoch)
        dump_fitness_nbhd(w);
    }
    check_knob_turn(w, last_fitness);
  }
  add_datum(w->epoch_fitness_deltas,
    ((find_best_fitness(w) - epoch_start_fitness) /
     (11.0 - epoch_start_fitness)));
  fflush(stdout);
}

void dump_virtual_fitness_func(World *w) {
  int best_organism_index = find_best_organism(w);
  //Organism o;
  //copy_organism(&o, &w->organisms[best_organism_index]);
  Organism *o = copy_organism(w->organisms[best_organism_index]);
  double delta = 0.02;
  puts("BEGIN VFUNC");
  for (double g1 = -1.0; g1 <= 1.0; g1 += delta) {
    for (double g2 = -1.0; g2 <= 1.0; g2 += delta) {
      o->genotype->nodes[0].initial_activation = g1;
      o->genotype->nodes[1].initial_activation = g2;
      sa(o, w->sa_timesteps, w->decay, w->spreading_rate);
      o->fitness = w->phenotype_fitness_func(w, o->genotype);
      printf("%lf %lf %lf %lf %lf\n",
        g1,
        g2,
        o->genotype->nodes[2].final_activation,
        o->genotype->nodes[3].final_activation,
        o->fitness);
    }
  }
  puts("END VFUNC");
}

void dump_phenotype_fitness_func(World *w) {
  double delta = 0.02;
  //Genotype g;
  //init_random_genotype(w, &g, 0, 4, 2, 2);
  Genotype *g = create_random_genotype(w);
  puts("BEGIN PHFUNC");
  for (double p1 = -1.0; p1 <= 1.0; p1 += delta) {
    for (double p2 = -1.0; p2 <= 1.0; p2 += delta) {
      g->nodes[2].final_activation = p1;
      g->nodes[3].final_activation = p2;
      double fitness = w->phenotype_fitness_func(w, g);
      printf("%lf %lf %lf\n", p1, p2, fitness);
    }
  }
  puts("END PHFUNC");
}

void print_world_params(World *w) {
  printf("w->seed=%d;\n", w->seed);
  printf("w->ridge_type=");
  switch (w->ridge_type) {
    case LINE:
      printf("LINE;\n");
      break;
    case CIRCLE:
      printf("CIRCLE;\n");
      break;
    default:
      assert(false);
      break;
  }
  printf("w->ridge_radius=%lf;\n", w->ridge_radius);
  printf("w->c2=%lf; w->c3=%lf;\n", w->c2, w->c3);
  printf("w->c1_lb=%lf; w->c1_ub=%lf;\n", w->c1_lb, w->c1_ub);
  printf("w->peak_movement=");
  switch (w->peak_movement) {
  case JUMPY_PEAK_MOVEMENT:
    printf("JUMPY_PEAK_MOVEMENT;\n");
    break;
  case GRADUAL_PEAK_MOVEMENT:
    printf("GRADUAL_PEAK_MOVEMENT;\n");
    break;
  default:
    assert(false);
  }
  printf("w->decay=%lf;\n", w->decay);
  printf("w->spreading_rate=%lf;\n", w->spreading_rate);
  printf("w->distance_weight=%lf;\n", w->distance_weight);
  printf("w->bumps=%s;\n", w->bumps ? "true" : "false");
  printf("w->knob_type=%d\n;", w->knob_type);
  printf("w->knob_constant=%lf\n;", w->knob_constant);
  printf("w->mutation_type_ub=%d;\n", w->mutation_type_ub);
  printf("w->extra_mutation_rate=%lf;\n", w->extra_mutation_rate);
  printf("w->crossover_freq=%lf;\n", w->crossover_freq);
  putchar('\n');
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
    case INHERIT_HALF_OF_CROSSOVER_EDGES:
      printf("INHERIT_HALF_OF_CROSSOVER_EDGES");
      break;
    case INHERIT_HALF_OF_ALL_EDGES:
      printf("INHERIT_HALF_OF_ALL_EDGES");
      break;
    case INHERIT_ALL_EDGES:
      printf("INHERIT_ALL_EDGES");
      break;
  }
  puts(";");
  printf("w->edge_weights=%s;\n", edge_weights_string(w->edge_weights));
  printf("w->multi_edges=%s;\n", w->multi_edges ? "true" : "false");
  printf("w->allow_move_edge=%s;\n", w->allow_move_edge ? "true" : "false");
  printf("w->activation_types=%s;\n",
      activation_types_string(w->activation_types));
  putchar('\n');
  printf("w->output_types=");
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      puts("ONLY_PASS_THROUGH;");
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      puts("PASS_THROUGH_AND_STEEP_SIGMOID;");
      break;
    default:
      assert(false);
  }
  printf("w->num_organisms=%d;\n", w->num_organisms);
  printf("w->num_candidates=%d;\n", w->num_candidates);
  printf("w->generations_per_epoch=%d;\n", w->generations_per_epoch);
  printf("w->num_epochs=%d;\n", w->num_epochs);
  printf("w->sa_timesteps=%d;\n", w->sa_timesteps);
  printf("w->num_hill_climbers=%d;\n", w->num_hill_climbers);
  printf("w->num_nodes=%d;\n", w->num_nodes);
  printf("w->num_edges=%d;\n", w->num_edges);
}

void print_acclivity_measures_of_best(World *w) {
  int best_organism_index = find_best_organism(w);
  //print_phenotype(w->organisms[best_organism_index]->genotype);
  HILL_CLIMBING_RESULT gvector_result = measure_acclivity(w, w->organisms[best_organism_index]);
  printf("gvector_fitness_delta = %lf, gvector_absolute_fitness = %lf\n", gvector_result.fitness_delta, gvector_result.ending_fitness);
  HILL_CLIMBING_RESULT phenotype_result = phenotype_acclivity(w);
  printf("acclivity: ph_fitness_delta = %lf, ph_absolute_fitness = %lf\n", phenotype_result.fitness_delta, phenotype_result.ending_fitness);
}

void print_knob_fitness_numbers(World *w) {
  printf("pos_knob_turns = %lf\n", w->num_fitness_increases_from_knob_turn /
      (double) w->num_generations_measured);
}

void run_world(World *w) {
  printf("--------------------------------------------------------------------------------\n");
  print_world_params(w);

  open_ancestor_log(w->log);

  srand(w->seed);
  init_random_population(w);
  for (int e = 1; e <= w->num_epochs; e++) {
    run_epoch(w, e);
  }
  if (!quiet) {
    dump_virtual_fitness_func(w);
    dump_phenotype_fitness_func(w);
  }
  printf("epoch fitness deltas: ");
  print_stats(w->epoch_fitness_deltas);
  //print_data(w->epoch_fitness_deltas);
  print_acclivity_measures_of_best(w);
  print_knob_fitness_numbers(w);

  close_ancestor_log(w->log);
  //free_world();
}

// -- fitness ----------------------------------------------------------------

double many_small_hills(double *phenotype) { // length is 2
  return cos(phenotype[0] * 30.0) * sin(phenotype[1] * 30.0);
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
  //printf("x = %lf; y = %lf; fabs(%lf) = %lf\n", x, y, y - (w->c2 * x + w->c3), fabs(y - (w->c2 * x + w->c3)));
  switch (w->ridge_type) {
  case LINE:
    return invv(0.0, w->ridge_radius, fabs(y - (w->c2 * x + w->c3)));
  case CIRCLE:
    return invv(0.0, w->ridge_radius, fabs((x*x + y*y) - (0.5*0.5)));
  default:
    assert(false);
  }
}

double require_valid_region(World *w, double x, double y) {
  switch (w->ridge_type) {
  case LINE:
    if (x >= w->c1_lb && x <= w->c1_ub &&
        //y >= w->c2 * w->c1_lb + w->c3 && y <= w->c2 * w->c1_ub + w->c3
        y >= w->c1_lb && y <= w->c1_ub
        )
      return 1.0;
    else
      return 0.0;
  case CIRCLE:
    return 1.0;
  default:
    assert(false);
  }
}

double phenotype_fitness(World *w, Genotype *g) {
  //const double sqrt8 = sqrt(8.0);
  //const double sqrt2 = sqrt(2.0);
  double phenotype[2] = {
    g->nodes[g->num_in].final_activation,
    g->nodes[g->num_in + 1].final_activation
  };
  //double peak_x = w->c1;
  //double peak_y = w->c2 * w->c1 + w->c3;
  if (verbose >= 2) {
    printf("require_valid_region(w, %lf, %lf) = %lf\n", phenotype[0], phenotype[1], require_valid_region(w, phenotype[0], phenotype[1]));
    printf("along_ridge(%lf, %lf) = %lf\n", phenotype[0], phenotype[1], along_ridge(w, phenotype[0], phenotype[1]));
  }
  double fitness = 0.0;
  if (phenotype[0] != UNWRITTEN && phenotype[1] != UNWRITTEN) {
    double dist = distance(w->peak_x, w->peak_y, phenotype[0], phenotype[1]);
    //double scaled_dist = (sqrt8 - dist) / sqrt8;  // 0.0 to 1.0; 1.0 is right on it
    //double scaled_dist = (sqrt2 - dist) / sqrt2;  // 0.0 to 1.0; 1.0 is right on it
    double scaled_dist = (w->max_dist - dist) / w->max_dist;
    fitness = //many_small_hills(phenotype) +
      //(5 * (sqrt8 - distance(w->c1, w->c1, phenotype[0], phenotype[1])));
      require_valid_region(w, phenotype[0], phenotype[1]) *
      along_ridge(w, phenotype[0], phenotype[1]) *
      (w->distance_weight * scaled_dist * scaled_dist);
    if (w->bumps) {
      double bump_amt = many_small_hills(phenotype);
//      if (bump_amt <= 0.0)
//        fitness = 0;
//      else
        fitness += bump_amt;
    }
  }
  if (fitness < 0.0) {
    return 0.0; //-10.0;
  } else {
    return fitness;
  }
}

// -- next generation via crossover and mutation -----------------------------

bool has_node(Genotype *g, int n) {
  return n >= 0 && n < g->num_nodes && g->nodes[n].in_use;
}

bool has_edge(Genotype *g, int src, int dst) {
  for (int e = 0; e < g->num_edges; e++) {
    if (g->edges[e].src == src && g->edges[e].dst == dst) {
      return true;
    }
  }
  return false;
}

void add_edge(World *w, Genotype *g, int src, int dst, double weight) {
  assert(has_node(g, src));
  assert(src >= 0);
  assert(has_node(g, dst));
  assert(dst >= 0);
  if (w->multi_edges || !has_edge(g, src, dst)) {
    g->num_edges++;
    g->edges = realloc(g->edges, sizeof(Edge) * g->num_edges);
    int e = g->num_edges - 1;
    g->edges[e].src = src;
    g->edges[e].dst = dst;
    g->edges[e].weight = weight;
  }
}

void mut_add_edge(World *w, Organism *o) {
  Genotype *g = o->genotype;
  add_edge(w, g, select_in_use_node(g), select_in_use_node(g), rand_edge_weight(w));
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

void mut_move_edge(World *w, Organism *o) {
  Genotype *g = o->genotype;
  if (g->num_edges == 0)
    return mut_add_edge(w, o);
  int selected_edge = rand() % g->num_edges;
  if (rand() & 1)
    g->edges[selected_edge].src = select_in_use_node(g);
  else
    g->edges[selected_edge].dst = select_in_use_node(g);
}

void mut_add_node(World *w, Organism *o) {
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
  switch (w->activation_types) {
    case ONLY_SUM_INCOMING:
      g->nodes[add_index].activation_type = SUM_INCOMING;
      break;
    case SUM_AND_MULT_INCOMING:
      g->nodes[add_index].activation_type = coin_flip() ? SUM_INCOMING : MULT_INCOMING;
      break;
    case SUM_AND_MIN_INCOMING:
      g->nodes[add_index].activation_type = coin_flip() ? SUM_INCOMING : MIN_INCOMING;
      break;
  }
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      g->nodes[add_index].output_type = PASS_THROUGH;
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      g->nodes[add_index].output_type = coin_flip() ? PASS_THROUGH : STEEP_SIGMOID;
      break;
  }
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
  o->from_turned_knob = true;
}

void mut_alter_activation_type(World *w, Organism *o) {
  switch (w->activation_types) {
  case ONLY_SUM_INCOMING:
    mut_turn_knob(w, o);
    break;
  case SUM_AND_MULT_INCOMING:
    {
      int n = select_in_use_node(o->genotype);
      Node *node_to_change = &o->genotype->nodes[n];
      switch (node_to_change->activation_type) {
        case SUM_INCOMING:
          node_to_change->activation_type = MULT_INCOMING;
          break;
        case MULT_INCOMING:
          node_to_change->activation_type = SUM_INCOMING;
          break;
        default:
          // shouldn't reach
          assert(false);
          break;
      }
    }
    break;
  case SUM_AND_MIN_INCOMING:
    {
      int n = select_in_use_node(o->genotype);
      Node *node_to_change = &o->genotype->nodes[n];
      switch (node_to_change->activation_type) {
        case SUM_INCOMING:
          node_to_change->activation_type = MIN_INCOMING;
          break;
        case MIN_INCOMING:
          node_to_change->activation_type = SUM_INCOMING;
          break;
        default:
          // shouldn't reach
          assert(false);
          break;
      }
    }
    break;
  }
}

void mut_alter_output_type(World *w, Organism *o) {
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        switch (node_to_change->output_type) {
          case PASS_THROUGH:
            node_to_change->output_type = STEEP_SIGMOID;
            break;
          case STEEP_SIGMOID:
            node_to_change->output_type = PASS_THROUGH;
            break;
          default:
            assert(false);
        }
      }
      break;
  }
}

Organism *mutate(World *w, Organism *old_o) {
  Organism *o = copy_organism(old_o);
  int num_mutations = 1 + (int)(w->extra_mutation_rate * rand_float() * (o->genotype->num_nodes + o->genotype->num_edges));
  for (int i = 0; i < num_mutations; i++) {
    int mutation_type = rand_int(0, w->mutation_type_ub); //rand() % 16;
    switch (mutation_type) {
    case 0:
      mut_add_node(w, o);
      log_mutation(w, "add_node");
      break;
    case 1:
      mut_remove_node(o);
      log_mutation(w, "remove_node");
      break;
    case 2:
      mut_add_edge(w, o);
      log_mutation(w, "add_edge");
      break;
    case 3:
      mut_remove_edge(o);
      log_mutation(w, "remove_edge");
      break;
    case 4:
      mut_alter_activation_type(w, o);
      log_mutation(w, "alter_act_type");
      break;
    case 5:
      mut_alter_output_type(w, o);
      log_mutation(w, "alter_out_type");
      break;
    case 6:
      if (w->allow_move_edge) {
        mut_move_edge(w, o);
        log_mutation(w, "move_edge");
        break;
      }
    default:
      mut_turn_knob(w, o);
      log_mutation(w, "turn_knob");
      break;
    }
  }
  return o;
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

//void crossover(World *w, Organism *b, Organism *m, Organism *d) {
Organism *crossover(World *w, Organism *m, Organism *d) {
  Genotype *mommy = m->genotype;
  Genotype *daddy = d->genotype;
  //Genotype *baby = b->genotype;
  Organism *baby_o = calloc(1, sizeof(Organism));
  Genotype *baby = calloc(1, sizeof(Genotype));
  baby_o->genotype = baby;
  baby_o->fitness = 0.0;
  baby_o->from_turned_knob = false;

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
          add_edge(w, baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_SRC_EDGES_FROM_MOMMY:
      case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
        if (edge->src < mommy_crossover_point &&
            has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(w, baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_CROSSOVER_EDGES:
        if ((edge->src < mommy_crossover_point || edge->dst < mommy_crossover_point) &&
            coin_flip() &&
            has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(w, baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_ALL_EDGES:
        if (coin_flip() &&
            has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(w, baby, edge->src, edge->dst, edge->weight);
        }
        break;
      case INHERIT_ALL_EDGES:
        if (has_node(baby, edge->src) && has_node(baby, edge->dst)) {
          add_edge(w, baby, edge->src, edge->dst, edge->weight);
        }
        break;
//    else if ((edge->src < baby->num_nodes && baby->nodes[edge->src].in_use) &&
//               (edge->dst < baby->num_nodes && baby->nodes[edge->dst].in_use) ) { //&&
               //!has_edge(baby, edge->src, edge->dst)) {
      //baby->edges[e++] = *edge;
//      add_edge(w, baby, edge->src, edge->dst, edge->weight);
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
          add_edge(w, baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_SRC_EDGES_FROM_MOMMY:
        // nothing
        break;
      case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
        if (edge->src >= daddy_crossover_point &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(w, baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_CROSSOVER_EDGES:
        if ((edge->src >= daddy_crossover_point || edge->dst >= daddy_crossover_point) &&
            coin_flip() &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(w, baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_HALF_OF_ALL_EDGES:
        if (coin_flip() &&
            !has_edge(mommy, bsrc, bdst) &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(w, baby, bsrc, bdst, edge->weight);
        }
        break;
      case INHERIT_ALL_EDGES:
        if (!has_edge(mommy, bsrc, bdst) &&
            has_node(baby, bsrc) && has_node(baby, bdst)) {
          add_edge(w, baby, bsrc, bdst, edge->weight);
        }
        break;
//        if (edge->src >= daddy_crossover_point &&
//            coin_flip() &&
//            has_node(baby, bsrc) && has_node(baby, bdst % baby->num_nodes)) {
//          add_edge(w, baby, bsrc, bdst % baby->num_nodes, edge->weight);
//        }
//        break;
    }
  }
  //assert(e == baby->num_edges);
  return baby_o;
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
  Genotype genotype = { nodes, edges, 6, 6, 6, 2, 2 };
  Organism o = { &genotype, 0.0 };
  
  verbose = 9;
  sa(&o, 13, 1.0, 1.0);
}

void sa_test2() {
  Node nodes[] = {
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, 0.0, 0.0, clamp },
    { true, -1.0, 0.0, clamp }
  };
  Edge edges[] = {
    { 4, 2, 1 },
    { 2, 2, 1 },
    { 2, 3, -1 },
    { 2, 3, -1 },
    { 2, 3, -1 },
    { 3, 3, 1 },
    { 3, 2, 1 }
  };
  Genotype genotype = { nodes, edges, 5, 5, 7, 2, 2 };
  Organism o = { &genotype, 0.0 };
  
  verbose = 9;
  sa(&o, 20, 1.0, 0.01);
}

void dump_phenotype_fitness() {
  World *w = create_world();
  w->num_organisms = 40;
  w->ridge_radius=0.200000;
  w->c2=1.000000; w->c3=0.000000;
  print_world_params(w);
  dump_phenotype_fitness_func(w);
}

void quick_test(int seed) {
  verbose = 1;
  World *w = create_world();
  w->num_organisms = 2;
  w->seed = seed;
  w->generations_per_epoch = 1;
  w->num_epochs = 1;
  w->num_nodes = 10;
  w->num_edges = 30;
  run_world(w);
}

void dot_test(int seed) {
  dot = true;
  World *w = create_world();
  w->num_organisms = 1;
  w->seed = seed;
  w->generations_per_epoch = 1;
  w->num_epochs = 1;
  w->num_nodes = 10;
  w->num_edges = 30;
  run_world(w);
}

void horizontal_ridge(World *w) {
  w->c2 = 0.0;
  w->c3 = 0.0;
}

void diagonal_ridge(World *w) {
  w->c2 = 1.0;
  w->c3 = 0.0;
}

void oblique_ridge(World *w) {
  w->c2 = 2.0;
  w->c3 = +0.45;
}

void easier_oblique_ridge(World *w) {
  w->c2 = 1.3;
  w->c3 = -0.45;
}

void long_test_start_small(int seed) {
  World *w = create_world();
  w->seed = seed;
  w->num_epochs = 200;
  w->sa_timesteps = 20;
  //w->generations_per_epoch = 20;
  //w->num_candidates = 5;
  w->edge_inheritance = INHERIT_ALL_EDGES;
  w->edge_weights = EDGE_WEIGHTS_ONLY_PLUS_1;
  //w->edge_weights = EDGE_WEIGHTS_POS_OR_NEG;
  //w->activation_types = ONLY_SUM_INCOMING;
  w->activation_types = SUM_AND_MULT_INCOMING;
  diagonal_ridge(w);
  w->bumps = true;
  w->ridge_radius = 0.05;
  w->crossover_freq = 0.3;
  //w->c1_lb = -1.0; w->c1_ub = +1.0;
  //w->mutation_type_ub = 30;
  //w->sa_timesteps = 20;
  //w->distance_weight = 5.0;
  w->dump_fitness_nbhd = true;
  w->dump_fitness_epoch = 5;
  w->dump_fitness_generation = 20;
  w->log->enabled = false;
  w->log->path = "./ancestors";
  run_world(w);
}

typedef struct {
  int id;
  int seed;
  double ridge_radius;
  enum {YX_RIDGE, OBLIQUE_RIDGE} ridge_type;
  int edge_inheritance;
} PARAMS;

void init_params(PARAMS *p, int seed) {
  p->id = 0;
  p->seed = seed;
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
  World *w = create_world();
  w->num_organisms = 40;
  w->seed = p->seed;
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
    INHERIT_HALF_OF_CROSSOVER_EDGES
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
  World *w = create_world();
  w->num_organisms = 40;
  w->seed = seed;
  w->generations_per_epoch = 5000;
  w->num_epochs = 1;
  run_world(w);
}

void good_run_oblique() {
  World *w = create_world();
  w->num_organisms = 40;
  w->seed=203540935;
  w->ridge_radius=0.200000;
  w->c2=2.000000; w->c3=0.450000;
  w->spreading_rate=0.010000;
  w->mutation_type_ub=15;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;
  //w->edge_inheritance=INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS;

  //w->num_organisms=40;
  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;

  w->num_epochs = 400;
  run_world(w);
}

void good_run_oblique2() {
  World *w = create_world();
  w->num_organisms = 40;
  w->seed=203540935;
  w->ridge_radius=0.200000;
  w->c2=2.000000; w->c3=0.450000;
  w->spreading_rate=0.010000;
  w->mutation_type_ub=16;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;
  //w->edge_inheritance=INHERIT_HALF_OF_EDGES_FROM_BOTH_PARENTS;

  //w->num_organisms=40;
  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;
  run_world(w);
}

void good_run_with_bumps() {
  World *w = create_world();
  w->num_organisms = 40;
  w->seed=23992348;
  w->ridge_radius=0.200000;
  w->c2=1.000000; w->c3=0.000000;
  w->spreading_rate=0.010000;
  w->distance_weight=10.000000;
  w->bumps=true;
  w->mutation_type_ub=16;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;

  w->edge_inheritance=INHERIT_ALL_EDGES;
  w->edge_weights=EDGE_WEIGHTS_ONLY_PLUS_1;
  w->activation_types=SUM_AND_MULT_INCOMING;

  w->num_organisms=40;
  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;
  run_world(w);
}

void tom() {
  World *w = create_world();
  w->seed=520664716;
  w->num_epochs = 200;
  w->ridge_radius=0.200000;
  w->c2=1.000000; w->c3=0.000000;
  w->decay=0.900000;
  w->spreading_rate=0.200000;
  w->distance_weight=10.000000;
  w->bumps=true;
  w->mutation_type_ub=16;
  w->extra_mutation_rate=0.100000;
  w->crossover_freq=0.300000;

  w->edge_inheritance=INHERIT_ALL_EDGES;
  w->edge_weights=EDGE_WEIGHTS_ONLY_PLUS_1;
  w->activation_types=SUM_AND_MULT_INCOMING;

  w->num_candidates=7;
  w->generations_per_epoch=20;
  w->sa_timesteps=20;
  run_world(w);
}

int get_seed(char **argv, int argc) {
  if (argc > 1) {
    return atoi(argv[1]);
  } else {
    return make_random_seed();
  }
}

void acclivation_test(int seed) {
  World *w = create_world();
  w->seed = seed;
  w->num_epochs = 20;
  w->sa_timesteps = 20;
  //w->generations_per_epoch = 20;
  //w->num_candidates = 5;
  w->edge_inheritance = INHERIT_ALL_EDGES;
  w->edge_weights = EDGE_WEIGHTS_ONLY_PLUS_1;
  //w->edge_weights = EDGE_WEIGHTS_POS_OR_NEG;
  //w->activation_types = ONLY_SUM_INCOMING;
  w->activation_types = SUM_AND_MULT_INCOMING;
  diagonal_ridge(w);
  w->bumps = true;
  w->ridge_radius = 0.2;
  w->crossover_freq = 0.3;
  //w->mutation_type_ub = 30;
  //w->sa_timesteps = 20;
  //w->distance_weight = 5.0;
  w->dump_fitness_nbhd = true;
  w->dump_fitness_epoch = 5;
  w->dump_fitness_generation = 20;
  w->log->enabled = false;
  w->log->path = "./ancestors";
  run_world(w);
  int best_organism_index = find_best_organism(w);
  print_phenotype(w->organisms[best_organism_index]->genotype);
  HILL_CLIMBING_RESULT gvector_result = measure_acclivity(w, w->organisms[best_organism_index]);
  printf("gvector: average delta = %lf, average fitness = %lf\n", gvector_result.fitness_delta, gvector_result.ending_fitness);
  HILL_CLIMBING_RESULT phenotype_result = phenotype_acclivity(w);
  printf("phenotype: average delta = %lf, average fitness = %lf\n", phenotype_result.fitness_delta, phenotype_result.ending_fitness);
}

void run_from_options(int argc, char **argv) {
  World *w = create_world();
  int option_index = 0;
  static struct option long_options[] = {
    { "seed", required_argument, 0, 0 },
    { "num_organisms", required_argument, 0, 0 },
    { "sa_timesteps", required_argument, 0, 0 },
    { "generations_per_epoch", required_argument, 0, 0 },
    { "num_epochs", required_argument, 0, 0 },
    { "num_nodes", required_argument, 0, 0 },
    { "num_edges", required_argument, 0, 0 },
    { "num_in", required_argument, 0, 0 },
    { "num_out", required_argument, 0, 0 },
    { "decay", required_argument, 0, 0 },
    { "spreading_rate", required_argument, 0, 0 },
    { "activation_types", required_argument, 0, 0 },
    { "edge_weights", required_argument, 0, 0 },
    { "distance_weight", required_argument, 0, 0 },
    { "bumps", required_argument, 0, 0 },
    { "c2", required_argument, 0, 0 },
    { "c3", required_argument, 0, 0 },
    { "c1_lb", required_argument, 0, 0 },
    { "c1_ub", required_argument, 0, 0 },
    { "ridge_radius", required_argument, 0, 0 },
    { "mutation_type_ub", required_argument, 0, 0 },
    { "extra_mutation_rate", required_argument, 0, 0 },
    { "crossover_freq", required_argument, 0, 0 },
    { "edge_inheritance", required_argument, 0, 0 },
    { "num_candidates", required_argument, 0, 0 },
    { "knob_constant", required_argument, 0, 0 },
    { "knob_type", required_argument, 0, 0 },
    { "num_hill_climbers", required_argument, 0, 0 },
    { "quiet", no_argument, 0, 0 },
    { "multi_edges", required_argument, 0, 0 },
    { "allow_move_edge", required_argument, 0, 0 },
    { "ridge_type", required_argument, 0, 0 },
    { "peak_movement", required_argument, 0, 0 },
    { "output_types", required_argument, 0, 0 },
    { NULL, 0, 0, 0 },
  };
  int c;

  for (;;) {
    c = getopt_long(argc, argv, "", long_options, &option_index);
    if (c == -1)
      break;
    if (c == 0) {
      switch (option_index) {
      case 0:
        w->seed = atoi(optarg);
        break;
      case 1:
        w->num_organisms = atoi(optarg);
        break;
      case 2:
        w->sa_timesteps = atoi(optarg);
        break;
      case 3:
        w->generations_per_epoch = atoi(optarg);
        break;
      case 4:
        w->num_epochs = atoi(optarg);
        break;
      case 5:
        w->num_nodes = atoi(optarg);
        break;
      case 6:
        w->num_edges = atoi(optarg);
        break;
      case 7:
        w->num_in = atoi(optarg);
        break;
      case 8:
        w->num_out = atoi(optarg);
        break;
      case 9:
        w->decay = atof(optarg);
        break;
      case 10:
        w->spreading_rate = atof(optarg);
        break;
      case 11:
        w->activation_types = atoi(optarg);
        break;
      case 12:
        w->edge_weights = atoi(optarg);
        break;
      case 13:
        w->distance_weight = atof(optarg);
        break;
      case 14:
        w->bumps = atoi(optarg);
        break;
      case 15:
        w->c2 = atof(optarg);
        break;
      case 16:
        w->c3 = atof(optarg);
        break;
      case 17:
        w->c1_lb = atof(optarg);
        break;
      case 18:
        w->c1_ub = atof(optarg);
        break;
      case 19:
        w->ridge_radius = atof(optarg);
        break;
      case 20:
        w->mutation_type_ub = atoi(optarg);
        break;
      case 21:
        w->extra_mutation_rate = atof(optarg);
        break;
      case 22:
        w->crossover_freq = atof(optarg);
        break;
      case 23:
        w->edge_inheritance = atoi(optarg);
        break;
      case 24:
        w->num_candidates = atoi(optarg);
        break;
      case 25:
        w->knob_constant = atof(optarg);
        break;
      case 26:
        w->knob_type = atoi(optarg);
        break;
      case 27:
        w->num_hill_climbers = atoi(optarg);
        break;
      case 28:
        quiet = true;
        break;
      case 29:
        w->multi_edges = atoi(optarg);
        break;
      case 30:
        w->allow_move_edge = atoi(optarg);
        break;
      case 31:
        w->ridge_type = atoi(optarg);
        break;
      case 32:
        w->peak_movement = atoi(optarg);
        break;
      case 33:
        w->output_types = atoi(optarg);
        break;
      default:
        printf("Internal error\n");
        exit(3);
      }
    } else {
      // unrecognized option
      exit(2);
    }
   }
   if (optind < argc) {
     printf("unknown args --\n");
     while (optind < argc)
       printf("  %s\n", argv[optind++]);
     exit(1);
   }

   run_world(w);
   free_world(w);
}

_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Wunused-variable\"")
int main(int argc, char **argv) {
  int seed = get_seed(argv, argc);
  run_from_options(argc, argv);
  //sa_test();
  //sa_test2();
  //quick_test(seed);
  //dot_test(seed);
  //long_test(seed);
  //long_test_start_small(seed);  //(677953487); // the main test
  //parameter_sweep(seed);
  //good_run_oblique();
  //good_run_oblique2();
  //one_long_epoch(seed);
  //dump_virt_test(seed);
  //dump_phenotype_fitness();
  //acclivation_test(374815447 /* seed */);
  //good_run_with_bumps();
  //tom();
  return 0;
}
_Pragma("GCC diagnostic pop")
