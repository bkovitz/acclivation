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
#include "sds.h"

#include "coordset.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAXS 128
#define array_len(a) (sizeof(a) / sizeof(a[0]))

#define UNWRITTEN -1000.0  // activation/output level that means that
                           // no activation/output level has been written
                           // to the node yet.
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

unsigned int make_random_seed() {
//  struct timespec tm;
//  clock_gettime(CLOCK_REALTIME, &tm);
//  return tm.tv_nsec;
  // https://developer.apple.com/library/content/documentation/Security/Conceptual/cryptoservices/RandomNumberGenerationAPIs/RandomNumberGenerationAPIs.html
  FILE *fp = fopen("/dev/random", "r");

  if (!fp) {
    perror("make_random_seed");
    exit(-1);
  }

  unsigned int value = 0;
  for (int i = 0; i < sizeof(value); i++) {
    value <<= 8;
    value |= fgetc(fp);
  }

  fclose(fp);
  return value;
}

bool coin_flip() {
  return rand() & 1;
}

// Lower bound and upper bound are both inclusive.
int rand_int(int lb, int ub) {
  return (rand() % (ub - lb + 1)) + lb;
}

double rand_double(double lb, double ub) {
  return (double)rand() / RAND_MAX * (ub - lb) + lb;
}

// float in range -1.0 to 1.0
double rand_activation() {
  return rand_double(-1.0, +1.0);
}

// -- Node parameters --------------------------------------------------------

// how inputs accumulate during one timestep
typedef enum {
  SUM_INCOMING,
  MULT_INCOMING,
  MIN_INCOMING
} INPUT_ACC;

// allowed INPUT_ACCs: bit 0 is SUM_INCOMING, bit 1 is MULT_INCOMING,
// and bit 2 is MIN_INCOMING
typedef unsigned int INPUT_ACCS;

const char *input_accs_string(INPUT_ACCS ia) {
  switch (ia) {
    case 0x01:
      return "0x01 /* SUM_INCOMING only */ ";
    case 0x02:
      return "0x02 /* MULT_INCOMING only */ ";
    case 0x03:
      return "0x03 /* SUM_INCOMING and MULT_INCOMING */ ";
    case 0x04:
      return "0x04 /* MIN_INCOMING only */ ";
    case 0x05:
      return "0x05 /* SUM_INCOMING and MIN_INCOMING */ ";
    case 0x06:
      return "0x06 /* MULT_INCOMING and MIN_INCOMING */ ";
    case 0x07:
      return "0x07 /* SUM_INCOMING, MULT_INCOMING, MIN_INCOMING */ ";
  }
  assert(false);
}

// How activation is calculated from accumulated inputs
typedef enum {
  SIGMOID,
  CLAMP_ONLY
  //SIGMOID_SUM
} ACTIVATION_TYPE;

// bits enable corresponding ACTIVATION_TYPEs
typedef unsigned int ACTIVATION_TYPES;

const char *activation_types_string(ACTIVATION_TYPES at) {
  switch (at) {
    case 0x01:
      return "0x01 /* SIGMOID only */ ";
    case 0x02:
      return "0x02 /* CLAMP_ONLY only */ ";
    case 0x03:
      return "0x03 /* SIGMOID and CLAMP_ONLY */ ";
  }
  assert(false);
}

// How/whether non-gvector nodes have an initial activation
typedef enum {
  NO_INITIAL_ACTIVATION,
  HAS_INITIAL_ACTIVATION
} INITIAL_ACTIVATION_TYPE;

const char *initial_activation_type_string(INITIAL_ACTIVATION_TYPE ia) {
  switch (ia) {
    case NO_INITIAL_ACTIVATION:
      return "NO_INITIAL_ACTIVATION";
    case HAS_INITIAL_ACTIVATION:
      return "HAS_INITIAL_ACTIVATION";
    default:
      assert(false);
  }
}

const char *input_acc_string(INPUT_ACC input_acc) {
  switch (input_acc) {
  case SUM_INCOMING:
    return "+";
  case MULT_INCOMING:
    return "*";
  case MIN_INCOMING:
    return "m";
  }
  assert(false);
}

//typedef enum {
//  ONLY_SUM_INCOMING,
//  SUM_AND_MULT_INCOMING,
//  SUM_AND_MIN_INCOMING,
//  ONLY_SIGMOID_SUM
//} ACTIVATION_TYPES;

//const char *activation_types_string(ACTIVATION_TYPES activation_types) {
//  switch (activation_types) {
//  case ONLY_SUM_INCOMING:
//    return "ONLY_SUM_INCOMING";
//  case SUM_AND_MULT_INCOMING:
//    return "SUM_AND_MULT_INCOMING";
//  case SUM_AND_MIN_INCOMING:
//    return "SUM_AND_MIN_INCOMING";
//  case ONLY_SIGMOID_SUM:
//    return "ONLY_SIGMOID_SUM";
//  default:
//    assert(false);
//  }
//}

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
  STEEP_SIGMOID,
  TWO_STEP
} OUTPUT_TYPE;

const char *output_type_string(OUTPUT_TYPE output_type) {
  switch (output_type) {
  case PASS_THROUGH:
    return "PASS_THROUGH";
  case STEEP_SIGMOID:
    return "STEEP_SIGMOID";
  case TWO_STEP:
    return "TWO_STEP";
  }
  assert(false);
}

typedef enum {
  ONLY_PASS_THROUGH,
  PASS_THROUGH_AND_STEEP_SIGMOID,
  PASS_THROUGH_AND_TWO_STEP
} OUTPUT_TYPES;

// How a Node's 'control' is updated each timestep
typedef enum {
  CONSTANT,         // 'control' doesn't change
  NUDGED_BY_INPUT   // 'control' is altered by first input
} CONTROL_UPDATE;

const char *control_update_string(CONTROL_UPDATE cu) {
  switch (cu) {
    case CONSTANT:
      return "CONSTANT";
    case NUDGED_BY_INPUT:
      return "NUDGED_BY_INPUT";
    default:
      assert(false);
  }
}

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

//double two_step(double x) {
//  if (x >= 0.2)
//    return 1.0;
//  else if (x <= -0.2)
//    return -1.0;
//  else
//    return x;
//}

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

// float in range 0.0 to 1.0
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
  if (verbose >= 5)
    printf("x = %lf; denom = %lf\n", x, denom);
  assert(denom != 0.0);
  return (yscale / denom) + yoffset;
}

double steep_sigmoid(double x, double xcenter, double slope) {
  //double xcenter = 0.0, ymin = -1.0, ymax = 1.0, slope = 4.0;
  //double ymin = 0.0, ymax = 1.0;
  double ymin = -1.0, ymax = +1.0;
  //double slope = 2.0;
  double yscale = ymax - ymin;
  double ycenter = (ymax + ymin) / 2.0;
  double yoffset = ycenter - (yscale / 2.0);
  double denom = (1.0 + exp(slope * (xcenter - x)));
  if (verbose >= 5)
    printf("x = %lf; denom = %lf\n", x, denom);
  assert(denom != 0.0);
  return (yscale / denom) + yoffset;
}

typedef struct {
  bool in_use;
  double initial_activation;
  double final_output;
  double final_activation;
  double control;  // a control parameter for the output_type function
  INPUT_ACC input_acc; // how inputs accumulate during one timestep
  ACTIVATION_TYPE activation_type; // activation as fn of accumulated inputs
  OUTPUT_TYPE output_type;  // output level as a function of activation level
  CONTROL_UPDATE control_update; // how 'control' is updated on each timestep
  INITIAL_ACTIVATION_TYPE initial_activation_type;
} Node;

typedef struct {
  int src;
  int dst;
  double weight;
} Edge;

// -- genotype ---------------------------------------------------------------

typedef struct {
  Node *nodes;  // array of num_nodes Nodes
  Edge *edges;  // array of num_edges Edges
  int num_nodes;
  int num_nodes_in_use;
  int num_edges;
  int num_in;
  int num_out;
} Genotype;

Node *get_node_i(Genotype *g, int i) {
  if (i >= g->num_nodes)
    return NULL;
  return &g->nodes[i];
}

Edge *get_edge_i(Genotype *g, int i) {
  if (i >= g->num_edges)
    return NULL;
  return &g->edges[i];
}

int get_edge_by_value(Genotype *g, int src, int dst, double weight) {
  for (int i = 0; i < g->num_edges; i++) {
    Edge *e = &g->edges[i];
    if (e->src == src && e->dst == dst && fabs(e->weight - weight) <= 0.001)
      return i;
  }
  return -1;
}

void set_gvector(Genotype *g, double x, double y) {
  g->nodes[0].initial_activation = x;
  g->nodes[1].initial_activation = y;
}

// MAYBE: shuffle and take without replace for deterministic select
int select_in_use_node(Genotype *g) {
  for ( ; ; ) {
    int index = rand() % g->num_nodes;
    if (g->nodes[index].in_use)
      return index;
  }
}

int select_in_use_removable_node(Genotype *g) {
  for ( ; ; ) {
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
    //printf("%4.4f ", g->nodes[p].final_activation);
    printf("%4.4f ", g->nodes[p].final_output);
  }
  printf("\n");
}

// -- organism ---------------------------------------------------------------

typedef enum {
  RANDOM,
  CROSSOVER,
  MUTATION
} BIRTH_TYPE;

typedef struct {
  int epoch;
  int generation;
  int org_index;
} ORGANISM_ID;

typedef struct {
  ORGANISM_ID mom;
  ORGANISM_ID dad;
  double crossover_frac;
} CROSSOVER_INFO;

typedef enum {
  MUT_ADD_NODE,
  MUT_REMOVE_NODE,
  MUT_ADD_EDGE,
  MUT_REMOVE_EDGE,
  MUT_ALTER_ACTIVATION_TYPE,
  MUT_ALTER_INPUT_ACC,
  MUT_ALTER_OUTPUT_TYPE,
  MUT_TURN_CONTROL,
  MUT_MOVE_EDGE,
  MUT_TURN_KNOB
} MUTATION_TYPE;

typedef struct {
  int index;
  double initial_activation;
  double control;
  INPUT_ACC input_acc;
  ACTIVATION_TYPE activation_type;
  OUTPUT_TYPE output_type;
  CONTROL_UPDATE control_update;
  INITIAL_ACTIVATION_TYPE initial_activation_type;
} NODE_MUTATION_RECORD;

typedef struct {
  int index;
  int src;
  int dst;
  double weight;
} EDGE_MUTATION_RECORD;

typedef struct {
  int index;
  int old_src;
  int src;
  int old_dst;
  int dst;
  double weight;
} MOVE_EDGE_MUTATION_RECORD;

typedef struct {
  int index;
  double old_value;
  double nudge;
  double value;
} KNOB_TURN_RECORD;

typedef struct {
  int index;
  ACTIVATION_TYPE old_type;
  ACTIVATION_TYPE type;
} ALTER_ACTIVATION_TYPE_RECORD;

typedef struct {
  int index;
  INPUT_ACC old_acc;
  INPUT_ACC acc;
} ALTER_INPUT_ACC_RECORD;

typedef struct {
  int index;
  OUTPUT_TYPE old_type;
  OUTPUT_TYPE type;
} ALTER_OUTPUT_TYPE_RECORD;

typedef struct {
  MUTATION_TYPE type;
  union {
    NODE_MUTATION_RECORD node_mutation;
    EDGE_MUTATION_RECORD edge_mutation;
    MOVE_EDGE_MUTATION_RECORD move_edge_mutation;
    KNOB_TURN_RECORD knob_turn;
    ALTER_ACTIVATION_TYPE_RECORD alter_activation_type;
    ALTER_INPUT_ACC_RECORD alter_input_acc;
    ALTER_OUTPUT_TYPE_RECORD alter_output_type;
  };
} MUTATION_RECORD;

typedef struct {
  ORGANISM_ID parent;
  int num_mutations;
  MUTATION_RECORD *mutations;
} MUTATION_INFO;

typedef struct {
  BIRTH_TYPE type;
  union {
    CROSSOVER_INFO crossover_info;
    MUTATION_INFO mutation_info;
  };
} BIRTH_INFO;

typedef struct {
  Genotype *genotype;
  double fitness;
  bool from_turned_knob;
  BIRTH_INFO birth_info;
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

Organism *copy_organism(Organism *o) {
  Organism *new_o = calloc(1, sizeof(Organism));
  new_o->genotype = copy_genotype(o->genotype);
  new_o->fitness = o->fitness;
  new_o->from_turned_knob = false;
  return new_o;
}

void print_id(ORGANISM_ID id) {
  printf("%d.%d.%d", id.epoch, id.generation, id.org_index);
}

void verify(bool);

void print_node_mutation_record(NODE_MUTATION_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    initial_activation=%lf\n", rec->initial_activation);
  printf("    control=%lf\n", rec->control);
  printf("    input_acc=%s\n", input_acc_string(rec->input_acc));
  printf("    activation_type=%s\n",
      initial_activation_type_string(rec->activation_type));
  printf("    output_type=%s\n", output_type_string(rec->output_type));
  printf("    control_update=%s\n", control_update_string(rec->control_update));
  printf("    initial_activation_type=%s\n",
      initial_activation_type_string(rec->initial_activation_type));
}

void print_edge_mutation_record(EDGE_MUTATION_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    src=%d\n", rec->src);
  printf("    dst=%d\n", rec->dst);
  printf("    weight=%lf\n", rec->weight);
}

void print_move_edge_mutation_record(MOVE_EDGE_MUTATION_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  if (rec->old_src != rec->src)
    printf("    src=%d (moved from %d)\n", rec->src, rec->old_src);
  else
    printf("    src=%d\n", rec->src);
  if (rec->old_dst != rec->dst)
    printf("    dst=%d (moved from %d)\n", rec->dst, rec->old_dst);
  else
    printf("    dst=%d\n", rec->dst);
  printf("    weight=%lf\n", rec->weight);
}

void print_alter_activation_type_record(ALTER_ACTIVATION_TYPE_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    activation_type=%s (was %s)\n",
      initial_activation_type_string(rec->old_type),
      initial_activation_type_string(rec->type));
}

void print_alter_input_acc_record(ALTER_INPUT_ACC_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    input_acc=%s (was %s)\n",
      input_acc_string(rec->old_acc),
      input_acc_string(rec->acc));
}

void print_alter_output_type_record(ALTER_OUTPUT_TYPE_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    output_type=%s (was %s)\n",
      output_type_string(rec->old_type),
      output_type_string(rec->type));
}

void print_knob_turn_record(KNOB_TURN_RECORD *rec) {
  printf("    index=%d\n", rec->index);
  printf("    old_value=%lf\n", rec->old_value);
  printf("    new_value=%lf\n", rec->value);
  printf("    nudge=%lf\n", rec->nudge);
}

void print_birth_record(Organism *o) {
  switch (o->birth_info.type) {
  case RANDOM:
    printf("  part of initial population\n");
    break;
  case MUTATION:
    {
      MUTATION_INFO *mi = &o->birth_info.mutation_info;
      printf("  from mutation [num_mutations=%d", mi->num_mutations);
      printf(" parent=");
      print_id(mi->parent);
      printf("]\n");
      for (int i = 0; i < mi->num_mutations; i++) {
        printf("  #%d: ", i + 1);
        MUTATION_RECORD *mr = &mi->mutations[i];
        switch (mr->type) {
        case MUT_ADD_NODE:
          printf("add node\n");
          print_node_mutation_record(&mr->node_mutation);
          break;
        case MUT_REMOVE_NODE:
          printf("remove node\n");
          print_node_mutation_record(&mr->node_mutation);
          break;
        case MUT_ADD_EDGE:
          printf("add edge\n");
          print_edge_mutation_record(&mr->edge_mutation);
          break;
        case MUT_REMOVE_EDGE:
          printf("remove edge\n");
          print_edge_mutation_record(&mr->edge_mutation);
          break;
        case MUT_ALTER_ACTIVATION_TYPE:
          printf("alter activation type\n");
          print_alter_activation_type_record(&mr->alter_activation_type);
          break;
        case MUT_ALTER_INPUT_ACC:
          printf("alter input acc\n");
          print_alter_input_acc_record(&mr->alter_input_acc);
          break;
        case MUT_ALTER_OUTPUT_TYPE:
          printf("alter output type\n");
          print_alter_output_type_record(&mr->alter_output_type);
          break;
        case MUT_TURN_CONTROL:
          printf("turn control\n");
          print_knob_turn_record(&mr->knob_turn);
          break;
        case MUT_MOVE_EDGE:
          printf("move edge\n");
          print_move_edge_mutation_record(&mr->move_edge_mutation);
          break;
        case MUT_TURN_KNOB:
          printf("turn knob\n");
          print_knob_turn_record(&mr->knob_turn);
          break;
        default:
          verify(false);
          break;
        }
      }
    }
    break;
  case CROSSOVER:
    {
      CROSSOVER_INFO *ci = &o->birth_info.crossover_info;
      printf("  from crossover [frac=%lf", ci->crossover_frac);
      printf(" mom=");
      print_id(ci->mom);
      printf(" dad=");
      print_id(ci->dad);
      printf("]\n");
      break;
    }
  default:
    verify(false);
    break;
  }
}

void print_organism(Organism *o) {
  Genotype *g = o->genotype;
  printf("  fitness=%.16lf\n  num_nodes=%d\n  num_edges=%d\n  g-vector=[%lf %lf]\n  phenotype=[%.16lf %.16lf]\n",
    o->fitness,
    g->num_nodes_in_use,
    g->num_edges,
    g->nodes[0].initial_activation,
    g->nodes[1].initial_activation,
    g->nodes[2].final_output,
    g->nodes[3].final_output);
  print_birth_record(o);
}

// -- world ------------------------------------------------------------------

typedef enum {
  LINE,
  CIRCLE
} RIDGE_TYPE;

typedef struct world_t {
  unsigned int seed;
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
  INITIAL_ACTIVATION_TYPE initial_activation_type;
  INPUT_ACCS input_accs;
  EDGE_WEIGHTS edge_weights;
  OUTPUT_TYPES output_types;
  CONTROL_UPDATE control_update; // all nodes get a copy of this
  bool multi_edges;
  bool allow_move_edge;
  Organism **organisms;
  double (*genotype_fitness_func)(struct world_t *, Genotype *);
  double distance_weight;
  double dist_exponent;
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
  bool reward_coverage;
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
  double control_increment;
  bool dump_fitness_nbhd;
  int dump_fitness_epoch;
  int dump_fitness_generation;
  DATA *epoch_fitness_deltas;
  FILE *log;
  int num_hill_climbers;
  int num_generations_measured;
  int num_fitness_increases_from_knob_turn;
  int param_set;
  int run;
  COORD_SET ridge_coords;
  bool invu; // smooth version of invv?
} World;

double phenotype_fitness(World *, const double *phenotype);
double genotype_fitness(World *, Genotype *);
double coverage_reward(World *w, Genotype *g);

// Returns an empty world with all default parameters.
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
  w->input_accs = 0x01; // SUM_INCOMING only
  w->activation_types = 0x01;  // SIGMOID only
  w->initial_activation_type = NO_INITIAL_ACTIVATION;
  w->edge_weights = EDGE_WEIGHTS_ONLY_PLUS_1; //EDGE_WEIGHTS_POS_OR_NEG;
  w->multi_edges = true;
  w->allow_move_edge = false;
  w->output_types = ONLY_PASS_THROUGH;
  w->control_update = CONSTANT;
  w->organisms = NULL; //calloc(num_organisms, sizeof(Organism));
  w->genotype_fitness_func = genotype_fitness;
  w->distance_weight = 10.0;
  w->dist_exponent = 1.0;
  w->bumps = true;
  w->epoch = 0;
  w->generation = 0;
  w->c1 = 0.5;
  w->c2 = 1.0;
  w->c3 = 0.0;
  w->c1_lb = -1.0;
  w->c1_ub = +1.0;
  w->peak_movement = JUMPY_PEAK_MOVEMENT;
  w->peak_x = 0.0; // dependent variable
  w->peak_y = 0.0; // dependent variable
  w->max_dist = 0.0; // dependent variable
  w->ridge_type = LINE;
  w->ridge_radius = 0.2;
  w->reward_coverage = false;
  w->mutation_type_ub = 10;
  w->extra_mutation_rate = 0.0; //0.1;
  w->crossover_freq = 0.02;
  w->edge_inheritance = INHERIT_SRC_EDGES_FROM_MOMMY;
  w->num_candidates = 7;
  w->knob_constant = 0.02;
  w->knob_type = KNOB_DISCRETE;
  w->control_increment = 0.2;
  w->dump_fitness_nbhd = false;
  w->dump_fitness_epoch = -1;
  w->dump_fitness_generation = -1;
  w->param_set = -1; // No param_set number provided on command line
  w->run = -1; // No run number provided on command line

  w->epoch_fitness_deltas = create_data();
  w->log = NULL;
  w->num_hill_climbers = 30;

  w->num_generations_measured = 0;
  w->num_fitness_increases_from_knob_turn = 0;

  w->ridge_coords = NULL;

  return w;
}

bool is_gvector_index(World *w, int index) {
  return index < w->num_in;
}

bool is_phenotype_index(World *w, int index) {
  return index >= w->num_in &&
         index < w->num_in + w->num_out;
}

double rand_initial_activation(World *w) {
  int half_range = 1.0 / w->knob_constant;
  return rand_int(-half_range, +half_range) * w->knob_constant;
}

void init_initial_activation(World *w, Node *node, int index) {
  switch (node->initial_activation_type) {
    case HAS_INITIAL_ACTIVATION:
      if (is_phenotype_index(w, index))
        node->initial_activation = UNWRITTEN;
      else
        node->initial_activation = rand_initial_activation(w);
      break;
    case NO_INITIAL_ACTIVATION:
      if (is_gvector_index(w, index))
        node->initial_activation = rand_initial_activation(w);
      else
        node->initial_activation = UNWRITTEN;
      break;
    default:
      assert(false);
  }
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

double node_output(Node *node, double activation);

void init_random_node(World *w, Node *n, int index) {
  n->in_use = true;
  switch (w->input_accs) {
    case 0x01:
      n->input_acc = SUM_INCOMING;
      break;
    case 0x02:
      n->input_acc = MULT_INCOMING;
      break;
    case 0x03:
      n->input_acc = coin_flip() ? SUM_INCOMING : MULT_INCOMING;
      break;
    case 0x04:
      n->input_acc = MIN_INCOMING;
      break;
    case 0x05:
      n->input_acc = coin_flip() ? SUM_INCOMING : MIN_INCOMING;
      break;
    case 0x06:
      n->input_acc = coin_flip() ? MULT_INCOMING : MIN_INCOMING;
      break;
    case 0x07:
      n->input_acc = rand_int(0, MIN_INCOMING);
      break;
    default:
      assert(false); // invalid input_accs
  }

  switch (w->activation_types) {
    case 0x01:
      n->activation_type = SIGMOID;
      break;
    case 0x02:
      n->activation_type = CLAMP_ONLY;
      break;
    case 0x03:
      n->activation_type = coin_flip() ? SIGMOID : CLAMP_ONLY;
      break;
    default:
      assert(false); // invalid activation_types
  }

  n->control_update = w->control_update;
  n->control = rand_activation();
  n->initial_activation_type = w->initial_activation_type;
  init_initial_activation(w, n, index);
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:  // SIGMOID_SUM should have PASS_THROUGH
      n->output_type = PASS_THROUGH;
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      n->output_type = coin_flip() ? PASS_THROUGH : STEEP_SIGMOID;
      break;
    case PASS_THROUGH_AND_TWO_STEP:
      n->output_type = coin_flip() ? PASS_THROUGH : TWO_STEP;
      break;
  }
  n->final_activation = UNWRITTEN;
  n->final_output = UNWRITTEN;
}

void add_edge(World *, Genotype *, int, int, double);
double node_output(Node *node, double activation);

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
    init_random_node(w, &g->nodes[n], n);
  }

  for (int e = 0; e < w->num_edges; e++) {
    int src = rand() % g->num_nodes;
    int dst = rand() % g->num_nodes;
    double weight = rand_edge_weight(w);
    add_edge(w, g, src, dst, weight);
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
  o->from_turned_knob = false; // FIXME: convert to use birth info
  o->birth_info.type = RANDOM;
  return o;
}

// TODO Incorporate input_acc
const char *node_type_string(const Node *node) {
  switch (node->activation_type) {
    case SIGMOID:
      return "S";
    case CLAMP_ONLY:
      return "C";
  }
  assert(false);
}

bool is_affected_by_control(Node *node) {
  return node->activation_type == SIGMOID ||
         node->output_type == STEEP_SIGMOID ||
         node->output_type == TWO_STEP;
}

void print_dot(World *w, Organism *o, FILE *f) {
  Genotype *g = o->genotype;

  sds digraph_name = sdsnew("digraph ");
  if (w->param_set >= 0)
    digraph_name = sdscatprintf(digraph_name, "p%d", w->param_set);
  if (w->run >= 0)
    digraph_name = sdscatprintf(digraph_name, "r%d", w->run);
  digraph_name = sdscatprintf(digraph_name, "e%dg%d", w->epoch, w->generation);
  digraph_name = sdscat(digraph_name, " {\n");
  fputs(digraph_name, f);
  sdsfree(digraph_name);

  fprintf(f, "  { rank=source edge [style=\"invis\"] ");
  for (int i = 0; i < g->num_in - 1; i++)
    fprintf(f, "n%d ->", i);
  fprintf(f, " n%d }\n", g->num_in - 1);
  fprintf(f, "  { rank=sink edge [style=\"invis\"] ");
  for (int o = 0; o < g->num_out - 1; o++)
    fprintf(f, "n%d ->", g->num_in + o);
  fprintf(f, " n%d }\n", g->num_in + g->num_out - 1);

  for (int n = 0; n < g->num_nodes; n++) {
    Node *node = &g->nodes[n];
    if (node->in_use) {
      fprintf(f, "  n%d [label=\"n%d (%s", n, n, node_type_string(node));
      fputs(input_acc_string(node->input_acc), f);
      switch (node->output_type) {
        case STEEP_SIGMOID:
          fputc('S', f);
          break;
        case TWO_STEP:
          fputc('T', f);
          break;
        default:
          // nothing
          break;
      }
      fputc(')', f);

      if (node->initial_activation != UNWRITTEN) {
        fprintf(f, " i=%.3lf", node->initial_activation);
      }
      if (node->final_activation != UNWRITTEN) {
        fprintf(f, " %.3lf", node->final_activation);
      }

      switch (node->output_type) {
        case STEEP_SIGMOID:
        case TWO_STEP:
          fprintf(f, " o=%.3lf", node->final_output);
          break;
        default:
          // nothing
          break;
      }

      if (is_affected_by_control(node))
        fprintf(f, " c=%.3lf", node->control);

      fprintf(f, "\"]\n");
    }
  }

  for (int e = 0; e < g->num_edges; e++) {
    fprintf(f, "  n%d -> n%d [label=%.3lf];\n",
        g->edges[e].src, g->edges[e].dst, g->edges[e].weight);
  }
  fprintf(f, "}\n");
}

// -- spreading activation ---------------------------------------------------

void init_activations(Genotype *g, double *activations) {
  for (int n = 0; n < g->num_nodes; n++) {
    activations[n] = g->nodes[n].initial_activation;
  }
}

void print_all_activations(FILE *f, Genotype *g, double *activations) {
  for (int n = 0; n < g->num_nodes; n++) {
    if (g->nodes[n].in_use)
      fprintf(f, "%.16f ", activations[n]);
  }
  fprintf(f, "\n");
}

double node_output(Node *node, double activation) {
  if (activation == UNWRITTEN) {
    return UNWRITTEN;
  } else {
    switch (node->output_type) {
      case PASS_THROUGH:
        return activation;
      case STEEP_SIGMOID:
        //return steep_sigmoid(activation, 0.0, node->control);
        return steep_sigmoid(activation, node->control / 4.0, node->control);
      case TWO_STEP:
        if (node->control > 0.0) {
          if (activation > node->control)
            return 1.0;
          else if (activation > 0.0)
            return 0.0;
          else
            return -1.0;
        } else {
          if (activation < node->control)
            return -1.0;
          else if (activation < 0.0)
            return  0.0;
          else
            return 1.0;
        }
      default:
        assert(false);
    }
  }
}

double src_output(Node *src, int index, double *activations) {
  return node_output(src, activations[index]);
}

double control_value(double input) {
  return exp(2.0 * input);
}

//void sa(Genotype *g, int timesteps, double decay, double spreading_rate) {
void sa(World *w, Genotype *g, FILE *outf) {
  int timesteps = w->sa_timesteps;
  double decay = w->decay;
  double spreading_rate = w->spreading_rate;
  //Genotype *g = o->genotype;

  // TODO refactor/simplify
  double activations[g->num_nodes];
  memset(activations, 0, sizeof(double) * g->num_nodes);
  init_activations(g, activations);
  
  double incoming_activations[g->num_nodes];
  double next_controls[g->num_nodes];

  if (outf) {
    fputs("initial: ", outf);
    print_all_activations(outf, g, activations);
  }

  for (int timestep = 1; timestep <= timesteps; timestep++) {
    // Initialize incoming_activations and next_controls for this timestep
    for (int n = 0; n < g->num_nodes; n++) {
      incoming_activations[n] = UNWRITTEN;
      next_controls[n] = UNWRITTEN;
    }

    // Load incoming_activations[] and next_controls[] with values coming
    // in to nodes via their edges
    for (int e = 0; e < g->num_edges; e++) {
      Edge *edge = &g->edges[e];
      if (activations[edge->src] != UNWRITTEN) {
        double incoming_output = src_output(&g->nodes[edge->src], edge->src,
            activations);
        switch (g->nodes[edge->dst].control_update) {
          case NUDGED_BY_INPUT:
            if (next_controls[edge->dst] == UNWRITTEN) {
              next_controls[edge->dst] = clamp(
                  g->nodes[edge->dst].control +
                  incoming_output / spreading_rate);
              break;
            }
            // fall through
          case CONSTANT:
            switch (g->nodes[edge->dst].input_acc) {
              case SUM_INCOMING:
                if (incoming_activations[edge->dst] == UNWRITTEN)
                  incoming_activations[edge->dst] = 0.0;
                incoming_activations[edge->dst] +=
                      edge->weight * incoming_output;
                break;
              case MULT_INCOMING:
                if (incoming_activations[edge->dst] == UNWRITTEN)
                  incoming_activations[edge->dst] = 1.0;
                incoming_activations[edge->dst] *=
                      edge->weight * incoming_output;
                break;
              case MIN_INCOMING:
                if (incoming_activations[edge->dst] == UNWRITTEN ||
                    incoming_output < incoming_activations[edge->dst]) {
                  incoming_activations[edge->dst] = incoming_output;
                }
                break;
            }
            break;
          default:
            assert(false);
        }
      }
    }

    // Update each Node's 'control' from next_controls
    for (int n = 0; n < g->num_nodes; n++) {
      Node *node = &g->nodes[n];
      if (node->in_use &&
          node->control_update == NUDGED_BY_INPUT &&
          next_controls[n] != UNWRITTEN) {
        node->control = next_controls[n];
      }
    }

    if (outf) {
      fputs("incoming: ", outf);
      print_all_activations(outf, g, incoming_activations);
    }

    // Load activations[] from incoming_activations[]
    for (int n = 0; n < g->num_nodes; n++) {
      Node *node = &g->nodes[n];
      if (node->in_use) {
        if (incoming_activations[n] != UNWRITTEN) {
          if (activations[n] == UNWRITTEN) {
            activations[n] = 0.0;
          }
          double a;
          switch (node->input_acc) {
            case SUM_INCOMING:
              a = activations[n];
              break;
            case MULT_INCOMING:
            case MIN_INCOMING:
              a = 0.0;
              break;
          }
          double x = spreading_rate *
                     incoming_activations[n] *
                     pow(decay, timestep - 1);
          switch (node->activation_type) {
            case CLAMP_ONLY:
              activations[n] = clamp(a + x);
              break;
            case SIGMOID:
              {
                double slope = 6.0 * node->control;
                //if (verbose > 2) {
                  //printf("activations[%d] = %.16lf\n", n, activations[n]);
                  //printf("steep_sigmoid(%.16lf, 0.0, %.16lf) = %.16lf\n",
                    //x, slope, steep_sigmoid(x, 0.0, slope));
                //}
                activations[n] = clamp(
                  a + steep_sigmoid(x, 0.0, slope));
                  //activations[n] + steep_sigmoid(x, node->control, slope));
                  //activations[n] + steep_sigmoid(x, node->control, node->control));
              }
          }
        }
      }
    }
    //if (verbose >= 9) {
    if (outf) {
      fprintf(outf, "timestep: %d\n", timestep);
      fputs("final: ", outf);
      print_all_activations(outf, g, activations);
    }
  }

  // Update each Node's 'final_activation' and 'final_output' from activations[]
  for (int n = 0; n < g->num_nodes; n++) {
    Node *node = &g->nodes[n];
    if (node->in_use) {
      node->final_activation = activations[n];
      node->final_output = node_output(node, node->final_activation);
    }
  }

  if (verbose)
    print_phenotype(g);
}

// ------------------------------------------------------------------------

void sanity_check(World *w);

void set_phenotypes_and_fitnesses(World *w) {
  for (int n = 0; n < w->num_organisms; n++) {
    Organism *o = w->organisms[n];
    //sa(o->genotype, w->sa_timesteps, w->decay, w->spreading_rate);
    sa(w, o->genotype, verbose ? stdout : NULL);
    if (debug)
      sanity_check(w);
//    if (dot)
//      print_dot(w, o, stdout);
    o->fitness = w->genotype_fitness_func(w, o->genotype);
  }
}

void init_random_population(World *w) {
  free_organisms(w->organisms, w->num_organisms);
  w->organisms = calloc(w->num_organisms, sizeof(Organism *));

  for (int n = 0; n < w->num_organisms; n++) {
    w->organisms[n] = create_random_organism(w);
  }
  set_phenotypes_and_fitnesses(w);
}

Organism *mutate(World *w, int, Organism *);

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
Organism *crossover(World *, int, Organism *, int, Organism *);

void sanity_check_organism(World *w, Organism *o) {
  Genotype *g = o->genotype;
  assert(g);
  // check in/out nodes in use
  for (int i = 0; i < g->num_in + g->num_out; i++)
    assert(g->nodes[i].in_use);
  for (int i = 0; i < g->num_in; i++)
    assert(g->nodes[i].initial_activation != UNWRITTEN);
  for (int i = g->num_in; i < g->num_in + g->num_out; i++)
    assert(g->nodes[i].initial_activation == UNWRITTEN);
  for (int i = 0; i < g->num_nodes; i++) {
    Node *node = &g->nodes[i];
    if (is_affected_by_control(node))
      assert(node->control != UNWRITTEN);
  }
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

void sanity_check(World *w) {
  for (int p = 0; p < w->num_organisms; p++) {
    sanity_check_organism(w, w->organisms[p]);
  }
}

// -- ancestor logging -------------------------------------------------------

void verify(bool cond) {
  if (!cond) {
    fprintf(stderr, "verify failed!\n");
    exit(-1);
  }
}

void verify_msg(bool cond, char *msg) {
  if (!cond) {
    fprintf(stderr, "%s\n", msg);
    exit(-1);
  }
}

typedef struct {
  size_t world_size;
  size_t genotype_size;
  size_t organism_size;
  size_t node_size;
  size_t edge_size;
  size_t mutation_record_size;
} STRUCT_SIZES;

void log_preamble(World *w) {
  if (w->log) {
    STRUCT_SIZES sizes = {
      sizeof(World),
      sizeof(Genotype),
      sizeof(Organism),
      sizeof(Node),
      sizeof(Edge),
      sizeof(MUTATION_RECORD),
    };
    size_t sizes_size = sizeof(STRUCT_SIZES);
    verify(1 == fwrite(&sizes_size, sizeof(size_t), 1, w->log));
    verify(1 == fwrite(&sizes, sizeof(sizes), 1, w->log));
    verify(1 == fwrite(w, sizeof(World), 1, w->log));
  }
}

void update_dependent_fitness_variables(World *);

World *load_preamble(FILE *f) {
  size_t sizes_size;
  STRUCT_SIZES sizes;
  verify(1 == fread(&sizes_size, sizeof(size_t), 1, f));
  verify_msg(sizes_size == sizeof(STRUCT_SIZES), "log created with incompatible preamble!");
  verify(1 == fread(&sizes, sizeof(sizes), 1, f));
  verify_msg(sizeof(World) == sizes.world_size, "log created with incompatible world!");
  verify_msg(sizeof(Genotype) == sizes.genotype_size, "log created with incompatible genotype!");
  verify_msg(sizeof(Organism) == sizes.organism_size, "log created with incompatible organism!");
  verify_msg(sizeof(Node) == sizes.node_size, "log created with incompatible node!");
  verify_msg(sizeof(Edge) == sizes.edge_size, "log created with incompatible edge!");
  verify_msg(sizeof(MUTATION_RECORD) == sizes.mutation_record_size, "log created with incompatible mutation record!");
  World *w = create_world();
  verify(1 == fread(w, sizeof(World), 1, f));
  w->organisms = NULL;
  w->genotype_fitness_func = genotype_fitness; // FIXME: can't load from file
  w->epoch_fitness_deltas = NULL;
  w->log = NULL;
  update_dependent_fitness_variables(w);
  return w;
}

void log_generation(World *w) {
  if (w->log) {
    verify(1 == fwrite(&w->epoch, sizeof(w->epoch), 1, w->log));
    verify(1 == fwrite(&w->c1, sizeof(w->c1), 1, w->log)); // ugly
    verify(1 == fwrite(&w->generation, sizeof(w->epoch), 1, w->log));
    // log organsims
    for (int i = 0; i < w->num_organisms; i++) {
      Organism *o = w->organisms[i];
      verify(1 == fwrite(o, sizeof(Organism), 1, w->log));
    }
    // log genotypes
    for (int i = 0; i < w->num_organisms; i++) {
      Genotype *g = (w->organisms[i])->genotype;
      verify(1 == fwrite(g, sizeof(Genotype), 1, w->log));
    }
    // log nodes/edges
    for (int i = 0; i < w->num_organisms; i++) {
      Genotype *g = (w->organisms[i])->genotype;
      verify(g->num_nodes == fwrite(g->nodes, sizeof(Node), g->num_nodes, w->log));
      verify(g->num_edges == fwrite(g->edges, sizeof(Edge), g->num_edges, w->log));
    }
    // log mutation info
    for (int i = 0; i < w->num_organisms; i++) {
      Organism *o = w->organisms[i];
      if (o->birth_info.type == MUTATION) {
        int n = o->birth_info.mutation_info.num_mutations;
        verify(n == fwrite(o->birth_info.mutation_info.mutations, sizeof(MUTATION_RECORD), n, w->log));
      }
    }
    fflush(stdout);
  }
}

typedef struct {
  int num_organisms;
  Organism **organisms;
} GENERATION;

GENERATION *create_generation(World *w) {
  GENERATION *gen = calloc(1, sizeof(GENERATION));
  gen->num_organisms = w->num_organisms;
  gen->organisms = calloc(gen->num_organisms, sizeof(Organism *));
  return gen;
}

typedef struct {
  double c1; // OAOO violation - c1 in World is authoritative
  int num_generations;
  GENERATION **generations;
} EPOCH;

GENERATION *load_generation(World *w, int epoch, EPOCH *pepoch, int generation, FILE *f) {
  GENERATION *gen = create_generation(w);
  int loaded_epoch;
  verify(1 == fread(&loaded_epoch, sizeof(loaded_epoch), 1, f));
  verify_msg(epoch == loaded_epoch, "bad epoch number in ancestor file");
  double loaded_c1;
  verify(1 == fread(&loaded_c1, sizeof(loaded_c1), 1, f));
  if (pepoch->c1 == UNWRITTEN)
      pepoch->c1 = loaded_c1;
  int loaded_generation;
  verify(1 == fread(&loaded_generation, sizeof(loaded_generation), 1, f));
  verify_msg(generation == loaded_generation, "bad generation number in ancestor file");
  // load organisms
  Organism *loaded_o = calloc(gen->num_organisms, sizeof(Organism));
  verify(gen->num_organisms == fread(loaded_o, sizeof(Organism), gen->num_organisms, f));
  // load genotypes
  Genotype *loaded_g = calloc(gen->num_organisms, sizeof(Genotype));
  verify(gen->num_organisms == fread(loaded_g, sizeof(Genotype), gen->num_organisms, f));
  // load nodes/edges
  for (int i = 0; i < gen->num_organisms; i++) {
    gen->organisms[i] = &loaded_o[i];
    Organism *o = gen->organisms[i];
    o->genotype = &loaded_g[i];
    Genotype *g = o->genotype;
    g->nodes = calloc(g->num_nodes, sizeof(Node));
    verify(g->num_nodes == fread(g->nodes, sizeof(Node), g->num_nodes, f));
    g->edges = calloc(g->num_edges, sizeof(Edge));
    verify(g->num_edges == fread(g->edges, sizeof(Edge), g->num_edges, f));
  }
  // load mutation info
  for (int i = 0; i < gen->num_organisms; i++) {
    Organism *o = gen->organisms[i];
    if (o->birth_info.type == MUTATION) {
      int n = o->birth_info.mutation_info.num_mutations;
      o->birth_info.mutation_info.mutations = calloc(n, sizeof(MUTATION_RECORD));
      verify(n == fread(o->birth_info.mutation_info.mutations, sizeof(MUTATION_RECORD), n, f));
    }
  }
  return gen;
}

EPOCH *create_epoch(World *w) {
  EPOCH *epoch = calloc(1, sizeof(EPOCH));
  epoch->c1 = UNWRITTEN;
  epoch->num_generations = w->generations_per_epoch;
  epoch->generations = calloc(epoch->num_generations, sizeof(GENERATION *));
  return epoch;
}

typedef struct {
  World *w;
  int num_epochs;
  EPOCH **epochs;
} RUN;

RUN *create_run(World *w) {
  RUN *run = calloc(1, sizeof(RUN));
  run->w = w;
  run->num_epochs = w->num_epochs;
  run->epochs = calloc(run->num_epochs, sizeof(EPOCH *));
  return run;
}

#define SENTINEL_LOG_VALUE 0xdeadbeef

RUN *load_ancestor_file(char *path) {
  verify(path);
  FILE *f = fopen(path, "rb");
  verify_msg(f, "can't open ancestor file");
  World *w = load_preamble(f);
  RUN *run = create_run(w);
  for (int e = 1; e <= run->w->num_epochs; e++) {
    run->epochs[e-1] = create_epoch(w);
    EPOCH *epoch = run->epochs[e-1];
    for (int g = 1; g <= epoch->num_generations; g++) {
      epoch->generations[g-1] = load_generation(w, e, epoch, g, f);
    }
  }
  unsigned sentinel = 0;
  verify(1 == fread(&sentinel, sizeof(unsigned), 1, f));
  verify_msg(SENTINEL_LOG_VALUE == sentinel, "bad sentinel value at end of ancestor file");
  fclose(f);
  return run;
}

void close_ancestor_log(World *w) {
  if (w->log) {
    unsigned sentinel = SENTINEL_LOG_VALUE;
    verify(1 == fwrite(&sentinel, sizeof(unsigned), 1, w->log));
    fclose(w->log);
  }
}

int prev_generation(World *w) {
  return w->generation > 1 ? (w->generation - 1) : (w->generations_per_epoch - 1);
}

int maybe_prev_epoch(World *w) {
  return w->generation > 1 ? w->epoch : (w->epoch - 1);
}

int log_mut_index; // to avoid threading state through all the mut_* functions

void set_mutation_start(World *w, Organism *baby, int parent_index, int num_mutations) {
  baby->birth_info.type = MUTATION;
  ORGANISM_ID id = { maybe_prev_epoch(w), prev_generation(w), parent_index };
  baby->birth_info.mutation_info.parent = id;
  baby->birth_info.mutation_info.num_mutations = num_mutations;
  baby->birth_info.mutation_info.mutations = calloc(num_mutations,
      sizeof(MUTATION_RECORD));
  log_mut_index = 0;
}

void set_mutation(World *w, Organism *baby, MUTATION_RECORD *record) {
  assert(log_mut_index < baby->birth_info.mutation_info.num_mutations);
  baby->birth_info.mutation_info.mutations[log_mut_index++] = *record;
}

void set_crossover_info(World *w, Organism *baby, int mom_index, int dad_index, double crossover_frac) {
  baby->birth_info.type = CROSSOVER;
  baby->birth_info.crossover_info.crossover_frac = crossover_frac;
  ORGANISM_ID id = { maybe_prev_epoch(w), prev_generation(w), mom_index };
  baby->birth_info.crossover_info.mom = id;
  id.org_index = dad_index;
  baby->birth_info.crossover_info.dad = id;
}

void dump_fitness_nbhd(World *w);

void run_generation(World *w) {
  Organism **old_population = w->organisms;
  Organism **new_population = calloc(w->num_organisms, sizeof(Organism *));
  for (int p = 0; p < w->num_organisms; p++) {
    if (rand_float() <= w->crossover_freq) {
      int mommy = tournament_select(w);
      int daddy = tournament_select(w);
      new_population[p] = crossover(w, mommy, old_population[mommy], daddy, old_population[daddy]);
    } else {
      int selected_organism = tournament_select(w);
      new_population[p] = mutate(w, selected_organism, old_population[selected_organism]);
    }
  }
  free_organisms(w->organisms, w->num_organisms);
  w->organisms = new_population;
  if (debug)
    sanity_check(w);
  set_phenotypes_and_fitnesses(w);
  if (debug)
    sanity_check(w);
  if (!quiet)
    log_generation(w);
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

double find_best_fitness(World *w) {
  int best_organism_index = find_best_organism(w);
  return w->organisms[best_organism_index]->fitness;
}

double measure_coverage(World *w, Genotype *test_g);

void print_best_fitness(World *w) {
  int best_organism_index = find_best_organism(w);
  Organism *o = w->organisms[best_organism_index];
  double max_fitness = o->fitness;
  Genotype *g = o->genotype;
  printf("    best fitness=%20.16lf  index=%2d  nodes=%2d  edges=%2d  g-vector=[% lf % lf] phenotype=[% .16lf % .16lf] cov=%.2lf\n",
    max_fitness, best_organism_index,
    g->num_nodes_in_use,
    g->num_edges,
    g->nodes[0].initial_activation,
    g->nodes[1].initial_activation,
//    g->nodes[2].final_activation,
//    g->nodes[3].final_activation);
    g->nodes[2].final_output,
    g->nodes[3].final_output,
    0.0); // too slow
    //measure_coverage(w, o->genotype));
  //print_dot(w, o, stdout);
}

void print_generation_results(World *w) {
  if (!dot) {
    printf("  generation %d\n", w->generation);
    print_best_fitness(w);
  }
}

void print_genotype_c(Genotype *g, FILE *outf) {
  fprintf(outf, "Node nodes[%d] = {\n", g->num_nodes);
  for (int i = 0; i < g->num_nodes; i++) {
    Node *n = &g->nodes[i];
    fprintf(outf, "  { %s, %.16lf, %.16lf, %.16lf, %.16lf, %d, %d, %d, %d, %d },\n",
        boolstr(n->in_use),
        n->initial_activation, n->final_output, n->final_activation, n->control,
        n->input_acc, n->activation_type, n->output_type, n->control_update,
        n->initial_activation_type);
  }
  fprintf(outf, "};\nEdges edges[%d] = {\n", g->num_edges);
  for (int i = 0; i < g->num_edges; i++) {
    Edge *e = &g->edges[i];
    fprintf(outf, "  { %d, %d, %.16lf },\n",
        e->src, e->dst, e->weight);
  }
  fprintf(outf, "};\nGenotype g = { nodes, edges, %d, %d, %d, %d, %d };\n",
      g->num_nodes, g->num_nodes_in_use, g->num_edges, g->num_in, g->num_out);
}

void dump_organism_fitness_nbhd(World *w, Organism *original) {
  double m = 2;
  Organism *o = copy_organism(original);
  Genotype *g = o->genotype;
  printf("neighborhood:\n");
  for (double dx = -m * w->knob_constant; dx <= m * w->knob_constant; dx += w->knob_constant) {
    for (double dy = -m * w->knob_constant; dy <= m * w->knob_constant; dy += w->knob_constant) {
      g->nodes[0].initial_activation = original->genotype->nodes[0].initial_activation + dx;
      g->nodes[1].initial_activation = original->genotype->nodes[1].initial_activation + dy;
      //sa(o->genotype, w->sa_timesteps, w->decay, w->spreading_rate);
      sa(w, o->genotype, verbose ? stdout : NULL);
      o->fitness = w->genotype_fitness_func(w, o->genotype);
      printf("  % lf % lf % .16lf % .16lf % lf\n",
        dx,
        dy,
        g->nodes[2].final_output,
        g->nodes[3].final_output,
        o->fitness);
    }
  }
  free_organism(o);
}

void dump_fitness_nbhd(World *w) {
  int best_organism_index = find_best_organism(w);
  Organism *original = w->organisms[best_organism_index];
  dump_organism_fitness_nbhd(w, original);
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

//HILL_CLIMBING_RESULT climb_hill(World *w, Organism *o) {
//  sa(w, o->genotype, verbose ? stdout : NULL);
//  o->fitness = w->genotype_fitness_func(w, o->genotype);
//
//  const double starting_fitness = o->fitness;
//  double last_fitness = starting_fitness;
//  double nudge_amount = w->knob_constant;
//
//  int num_candidates = w->num_in << 1;
//  Organism *candidates[num_candidates];
//
//  int num_neutral_steps = 0;
//
//  //printf("starting_fitness = %lf\n", starting_fitness);
//  for ( ; ; ) {
//    // create organisms that take a step in each possible direction
//    for (int i = 0; i < num_candidates; i++) {
//      candidates[i] = copy_organism(o);
//      Organism *candidate = candidates[i];
//      nudge_candidate(candidate, i >> 1, (i & 1) ? nudge_amount : -nudge_amount);
//      sa(w, candidate->genotype, verbose ? stdout : NULL);
//      candidate->fitness = w->genotype_fitness_func(w, candidate->genotype);
//      //printf("candidate %d fitness = %lf\n", i, candidate->fitness);
//      //print_phenotype(candidate->genotype);
//    }
//    // take the best positive step, or bail out
//    int best_candidate_index = find_best_organism_in(candidates, num_candidates);
//    if (candidates[best_candidate_index]->fitness > last_fitness) {
//      free_organism(o);
//      o = copy_organism(candidates[best_candidate_index]);
//      last_fitness = o->fitness;
//      //printf("fitness, last fitness-> %lf, %lf\n", o->fitness, last_fitness);
//      for (int i = 0; i < num_candidates; i++)
//        free_organism(candidates[i]);
//      num_neutral_steps = 0;
//    } else if (candidates[best_candidate_index]->fitness == last_fitness) {
//      // neutral plateau
//      if (++num_neutral_steps >= 100) {
//        //printf("crazy plateau\n");
//        break;
//      }
//      free_organism(o);
//      o = copy_organism(candidates[rand_int(0, num_candidates - 1)]);
//      for (int i = 0; i < num_candidates; i++)
//        free_organism(candidates[i]);
//    } else {
//      // reached a peak
//      break;
//    }
//  }
//  HILL_CLIMBING_RESULT result =
//      { last_fitness - starting_fitness,
//        last_fitness };
//  free_organism(o);
//  return result;
//}

// Function that takes array of two doubles and returns fitness
typedef double (*FITNESS_FUNC)(World *, const double *);

HILL_CLIMBING_RESULT climb_hill2(
    World *w, FITNESS_FUNC ffunc, const double *startxy)
{
  double xy[2];
  xy[0] = startxy[0];
  xy[1] = startxy[1];

  struct Candidate {
    const double deltas[2];
    double xy[2];
    double fitness;
  } candidates[4] = {
    { {-w->knob_constant, 0.0}, {0.0, 0.0}, 0.0},
    { {+w->knob_constant, 0.0}, {0.0, 0.0}, 0.0},
    { {0.0, -w->knob_constant}, {0.0, 0.0}, 0.0},
    { {0.0, +w->knob_constant}, {0.0, 0.0}, 0.0}
  };

  const double starting_fitness = (*ffunc)(w, startxy);
  double last_fitness = starting_fitness;

  int num_neutral_steps = 0;
  struct Candidate *best_candidate;
  for ( ; ; ) {
    for (int i = 0; i < array_len(candidates); i++) {
      struct Candidate *c = &candidates[i];
      c->xy[0] = clamp(xy[0] + c->deltas[0]);
      c->xy[1] = clamp(xy[1] + c->deltas[1]);
      c->fitness = (*ffunc)(w, c->xy);

      if (i == 0)
        best_candidate = c;
      else {
        if (c->fitness > best_candidate->fitness)
          best_candidate = c;
        else if (c->fitness == best_candidate->fitness)
          best_candidate = coin_flip() ? best_candidate : c;
      }
    }

    if (best_candidate->fitness > last_fitness) {
      memmove(xy, best_candidate->xy, sizeof(xy));
      last_fitness = best_candidate->fitness;
      num_neutral_steps = 0;
    } else if (best_candidate->fitness == last_fitness) {
      if (num_neutral_steps >= 100) {
        break;
      } else {
        memmove(xy, best_candidate->xy, sizeof(xy));
        num_neutral_steps++;
      }
    } else {
      // xy is a peak
      break;
    }
  }

  //printf("last_fitness - starting_fitness = %lf\n", last_fitness - starting_fitness); //DEBUG
  HILL_CLIMBING_RESULT result =
      { last_fitness - starting_fitness,
        last_fitness };
  return result;
}

HILL_CLIMBING_RESULT measure_acclivity2(World *w, FITNESS_FUNC ffunc) {
  bool save_reward_coverage = w->reward_coverage;
  w->reward_coverage = false;

  HILL_CLIMBING_RESULT total = { 0.0, 0.0 };

  srand(0); //SHOULD push rng state

  for (int i = 0; i < w->num_hill_climbers; i++) {
    double xy[2] = { rand_initial_activation(w), rand_initial_activation(w) };
    HILL_CLIMBING_RESULT result = climb_hill2(w, ffunc, xy);
    total.fitness_delta += result.fitness_delta;
    total.ending_fitness += result.ending_fitness;
  }
  
  total.fitness_delta /= w->num_hill_climbers;
  total.ending_fitness /= w->num_hill_climbers;

  w->reward_coverage = save_reward_coverage;
  return total;
}

//HILL_CLIMBING_RESULT measure_acclivity(World *w, Organism *test_o) {
//  Organism *o;
//  HILL_CLIMBING_RESULT total = { 0.0, 0.0 };
//
//  srand(0);
//  for (int i = 0; i < w->num_hill_climbers; i++) {
//    o = copy_organism(test_o);
//    set_random_gvector(o);
//    //printf("gvector-> %lf, %lf\n", o->genotype->nodes[0].initial_activation, o->genotype->nodes[1].initial_activation);
//    HILL_CLIMBING_RESULT result = climb_hill(w, o);
//    total.fitness_delta += result.fitness_delta;
//    total.ending_fitness += result.ending_fitness;
//    //printf("->delta %lf, abs %lf\n", result.fitness_delta, result.ending_fitness);
//    //free_organism(o);
//  }
//
//  total.fitness_delta /= w->num_hill_climbers;
//  total.ending_fitness /= w->num_hill_climbers;
//  return total;
//}

//HILL_CLIMBING_RESULT phenotype_acclivity(World *w) {
//  Node nodes[] = {
//    { true, 0.0, 0.0, 0.0, 0.0, MULT_INCOMING, CLAMP_ONLY,
//      PASS_THROUGH, CONSTANT, HAS_INITIAL_ACTIVATION },
//    { true, 0.0, 0.0, 0.0, 0.0, MULT_INCOMING, CLAMP_ONLY,
//      PASS_THROUGH, CONSTANT, HAS_INITIAL_ACTIVATION },
//    { true, 0.0, 0.0, 0.0, 0.0, MULT_INCOMING, CLAMP_ONLY,
//      PASS_THROUGH, CONSTANT, HAS_INITIAL_ACTIVATION },
//    { true, 0.0, 0.0, 0.0, 0.0, MULT_INCOMING, CLAMP_ONLY,
//      PASS_THROUGH, CONSTANT, HAS_INITIAL_ACTIVATION }
//  };
//  Edge edges[] = {
//    { 0, 2, 1.0 },
//    { 1, 3, 1.0 }
//  };
//  Genotype genotype = { nodes, edges, 4, 4, 2, 2, 2 };
//  Organism null_organism = { &genotype, 0.0 };
//  
//  return measure_acclivity(w, &null_organism);
//}

HILL_CLIMBING_RESULT phenotype_acclivity(World *w) {
  return measure_acclivity2(w, phenotype_fitness);
}

HILL_CLIMBING_RESULT genotype_acclivity(World *w, Genotype *g) {
  double gfitness(World *w, const double *gvector) {
    Genotype *gg = copy_genotype(g);
    set_gvector(gg, gvector[0], gvector[1]);
    double fitness = genotype_fitness(w, g);
    //printf("gvector=(%lf, %lf) fitness=%lf\n", gvector[0], gvector[1], fitness);
    free(gg);
    return fitness;
  }
  return measure_acclivity2(w, gfitness);
}

int x2i(World *w, double x) {
  return (x - -1.0) / (2 * w->knob_constant);
}

void set_ridge_coords(World *w) {
  double phenotype[2];
  //TODO How to properly free old w->ridge_coords?
  w->ridge_coords = NULL;

  for (double x = -1.0; x <= 1.0; x += 2 * w->knob_constant) {
    phenotype[0] = x;
    for (double y = -1.0; y <= 1.0; y += 2 * w->knob_constant) {
      phenotype[1] = y;
      double fitness = phenotype_fitness(w, phenotype);
      if (fitness > 1.0) {
        int ix = x2i(w, x);
        int iy = x2i(w, y);
        //printf("x=%lf ix=%d  y=%lf iy=%d  fitness=%lf\n",
        //    x, ix, y, iy, fitness);
        add_coord(&w->ridge_coords, ix, iy);
      }
    }
  }
}

double measure_coverage(World *w, Genotype *test_g) {
  Genotype *g = copy_genotype(test_g);
  clear_coverage(w->ridge_coords);
  for (double x = -1.0; x <= 1.0; x += 2 * w->knob_constant) {
    for (double y = -1.0; y <= 1.0; y += 2 * w->knob_constant) {
      set_gvector(g, x, y);
      //sa(g, w->sa_timesteps, w->decay, w->spreading_rate);
      sa(w, g, verbose ? stdout : NULL);
      int ix = x2i(w, g->nodes[w->num_in].final_output);
      int iy = x2i(w, g->nodes[w->num_in + 1].final_output);
      cover_coord(w->ridge_coords, ix, iy);
    }
  }
  free_genotype(g);
  return calculate_coverage(w->ridge_coords);
}

// ----------------------------------------------------------------------

void free_world(World *w) {
  free_organisms(w->organisms, w->num_organisms);
  free_data(w->epoch_fitness_deltas);
  free(w);
}

double distance(double x1, double y1, double x2, double y2);

void update_dependent_fitness_variables(World *w) {
  switch (w->ridge_type) {
  case LINE:
    // dependent variables
    w->peak_x = w->c1;
    w->peak_y = clamp2(w->c2 * w->c1 + w->c3, w->c1_lb, w->c1_ub);
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
    if (verbose)
      printf("peak=(%lf,%lf) max_dist=%lf\n", w->peak_x, w->peak_y, w->max_dist);
    break;
  case CIRCLE:
    // dependent variables
    w->peak_x = 0.5 * cos(w->c1);
    w->peak_y = 0.5 * sin(w->c1);
    w->max_dist = 1.0;  // Furthest place on circle gets 0.0 fitness
    break;
  }
}

void change_fitness_constants(World *w) {
  //const double sqrt2 = sqrt(2.0);
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
    break;
  case CIRCLE:
    switch (w->peak_movement) {
    case JUMPY_PEAK_MOVEMENT:
      w->c1 = rand_double(0.0, 2 * M_PI);
      break;
    case GRADUAL_PEAK_MOVEMENT:
      {
        double delta = sample_normal(2 * M_PI / 10.0);
        //double delta = rand_double(2 * M_PI / 10, 2 * M_PI / 5.0);
        printf("c1 delta = %lf\n", delta);
        w->c1 += delta;
      }
      break;
    default:
      assert(false);
    }
    break;
  }
  update_dependent_fitness_variables(w);
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
    printf("\nepoch %d (c1=%lf, c2=%lf, c3=%lf) peak=(%.3lf, %.3lf)\n",
      e, w->c1, w->c2, w->c3, w->peak_x, w->peak_y);
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

/*void dump_virtual_fitness_func(World *w, FILE *outf) {
  int best_organism_index = find_best_organism(w);
  Organism *o = copy_organism(w->organisms[best_organism_index]);
  Genotype *g = o->genotype;
  double cov_reward = coverage_reward(w, g);
  bool save_reward_coverage = w->reward_coverage;
  w->reward_coverage = false;
  double delta = 0.02;
  fputs("BEGIN VFUNC", outf);
  for (double g1 = -1.0; g1 <= 1.0; g1 += delta) {
    for (double g2 = -1.0; g2 <= 1.0; g2 += delta) {
      set_gvector(g, g1, g2);
//      o->genotype->nodes[0].initial_activation = g1;
//      o->genotype->nodes[1].initial_activation = g2;
      //sa(o->genotype, w->sa_timesteps, w->decay, w->spreading_rate);
      sa(w, o->genotype, verbose ? stdout : NULL);
      o->fitness = w->genotype_fitness_func(w, o->genotype) + cov_reward;
      fprintf(outf, "%lf %lf %lf %lf %lf\n",
        g1,
        g2,
//        o->genotype->nodes[2].final_activation,
//        o->genotype->nodes[3].final_activation,
        o->genotype->nodes[2].final_output,
        o->genotype->nodes[3].final_output,
        o->fitness);
    }
  }
  fputs("END VFUNC", outf);
  w->reward_coverage = save_reward_coverage;
  free_organism(o);
  fflush(outf);
}*/

void dump_organism_virtual_fitness_func(World *w, Organism *org, bool delims, FILE *outf) {
  Organism *o = copy_organism(org);
  Genotype *g = o->genotype;
  double cov_reward = coverage_reward(w, g);
  bool save_reward_coverage = w->reward_coverage;
  w->reward_coverage = false;
  double delta = 0.02;
  if (delims)
    fputs("BEGIN VFUNC\n", outf);
  for (double g1 = -1.0; g1 <= 1.0; g1 += delta) {
    for (double g2 = -1.0; g2 <= 1.0; g2 += delta) {
      set_gvector(g, g1, g2);
//      o->genotype->nodes[0].initial_activation = g1;
//      o->genotype->nodes[1].initial_activation = g2;
      //sa(o->genotype, w->sa_timesteps, w->decay, w->spreading_rate);
      sa(w, o->genotype, verbose ? stdout : NULL);
      o->fitness = w->genotype_fitness_func(w, o->genotype) + cov_reward;
      fprintf(outf, "%lf %lf %lf %lf %lf\n",
        g1,
        g2,
//        o->genotype->nodes[2].final_activation,
//        o->genotype->nodes[3].final_activation,
        o->genotype->nodes[2].final_output,
        o->genotype->nodes[3].final_output,
        o->fitness);
    }
  }
  if (delims)
    fputs("END VFUNC\n", outf);
  w->reward_coverage = save_reward_coverage;
  free_organism(o);
  fflush(outf);
}

void dump_virtual_fitness_func(World *w, bool delims, FILE *outf) {
  int best_organism_index = find_best_organism(w);
  dump_organism_virtual_fitness_func(w, w->organisms[best_organism_index], delims, outf);
}

void dump_phenotype_fitness_func(World *w, bool delims, FILE *outf) {
  double delta = 0.02;
  //Genotype g;
  //init_random_genotype(w, &g, 0, 4, 2, 2);
  //Genotype *g = create_random_genotype(w);
  double phenotype[w->num_out];
  if (delims)
    fputs("BEGIN PHFUNC\n", outf);
  for (double p1 = -1.0; p1 <= 1.0; p1 += delta) {
    phenotype[0] = p1;
    for (double p2 = -1.0; p2 <= 1.0; p2 += delta) {
//      g->nodes[2].final_activation = p1;
//      g->nodes[3].final_activation = p2;
//      g->nodes[2].final_output = p1;
//      g->nodes[3].final_output = p2;
      phenotype[1] = p2;
      //double fitness = w->genotype_fitness_func(w, g);
      double fitness = phenotype_fitness(w, phenotype);
      fprintf(outf, "%lf %lf %lf\n", p1, p2, fitness);
    }
  }
  if (delims)
    fputs("END PHFUNC\n", outf);
  fflush(outf);
}

void print_world_params(World *w, FILE *outf) {
  fprintf(outf, "w->seed=%u;\n", w->seed);
  fprintf(outf, "w->run=%d;\n", w->run);
  fprintf(outf, "w->param_set=%d;\n", w->param_set);
  fputc('\n', outf);
  fputs("// Fitness\n", outf);
  fprintf(outf, "w->ridge_type=");
  switch (w->ridge_type) {
    case LINE:
      fprintf(outf, "LINE;\n");
      break;
    case CIRCLE:
      fprintf(outf, "CIRCLE;\n");
      break;
    default:
      assert(false);
      break;
  }
  fprintf(outf, "w->bumps=%s;\n", w->bumps ? "true" : "false");
  fprintf(outf, "w->ridge_radius=%lf;\n", w->ridge_radius);
  fprintf(outf, "w->c2=%lf; w->c3=%lf;\n", w->c2, w->c3);
  fprintf(outf, "w->c1_lb=%lf; w->c1_ub=%lf;\n", w->c1_lb, w->c1_ub);
  fprintf(outf, "w->distance_weight=%lf;\n", w->distance_weight);
  fprintf(outf, "w->peak_movement=");
  switch (w->peak_movement) {
  case JUMPY_PEAK_MOVEMENT:
    fprintf(outf, "JUMPY_PEAK_MOVEMENT;\n");
    break;
  case GRADUAL_PEAK_MOVEMENT:
    fprintf(outf, "GRADUAL_PEAK_MOVEMENT;\n");
    break;
  default:
    assert(false);
  }
  fprintf(outf, "w->reward_coverage=%s;\n", w->reward_coverage ? "true" : "false");
  fputc('\n', outf);
  fputs("// gvector\n", outf);
  fprintf(outf, "w->knob_type=%d;\n", w->knob_type);
  fprintf(outf, "w->knob_constant=%lf;\n", w->knob_constant);
  fputc('\n', outf);
  fputs("// Spreading activation\n", outf);
  fprintf(outf, "w->sa_timesteps=%d;\n", w->sa_timesteps);
  fprintf(outf, "w->spreading_rate=%lf;\n", w->spreading_rate);
  fprintf(outf, "w->decay=%lf;\n", w->decay);
  fprintf(outf, "w->initial_activation_type=%s;\n",
      initial_activation_type_string(w->initial_activation_type));
  fprintf(outf, "w->input_accs=%s;\n", input_accs_string(w->input_accs));
  fprintf(outf, "w->activation_types=%s;\n",
      activation_types_string(w->activation_types));
  fprintf(outf, "w->output_types=");
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      fputs("ONLY_PASS_THROUGH;", outf);
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      fputs("PASS_THROUGH_AND_STEEP_SIGMOID;", outf);
      break;
    case PASS_THROUGH_AND_TWO_STEP:
      fputs("PASS_THROUGH_AND_TWO_STEP;", outf);
      break;
    default:
      assert(false);
  }
  fprintf(outf, "w->control_update=%s;\n", control_update_string(w->control_update));
  fprintf(outf, "w->control_increment=%lf;\n", w->control_increment);
  fputc('\n', outf);
  fputs("// Mutation and crossover\n", outf);
  fprintf(outf, "w->mutation_type_ub=%d;\n", w->mutation_type_ub);
  fprintf(outf, "w->extra_mutation_rate=%lf;\n", w->extra_mutation_rate);
  fprintf(outf, "w->crossover_freq=%lf;\n", w->crossover_freq);
  fputc('\n', outf);
  fputs("// Edge inheritance\n", outf);
  fprintf(outf, "w->edge_inheritance=");
  switch (w->edge_inheritance) {
    case NO_EDGES_ACROSS_PARENTS:
      fprintf(outf, "NO_EDGES_ACROSS_PARENTS");
      break;
    case INHERIT_SRC_EDGES_FROM_MOMMY:
      fprintf(outf, "INHERIT_SRC_EDGES_FROM_MOMMY");
      break;
    case INHERIT_SRC_EDGES_FROM_BOTH_PARENTS:
      fprintf(outf, "INHERIT_SRC_EDGES_FROM_BOTH_PARENTS");
      break;
    case INHERIT_HALF_OF_CROSSOVER_EDGES:
      fprintf(outf, "INHERIT_HALF_OF_CROSSOVER_EDGES");
      break;
    case INHERIT_HALF_OF_ALL_EDGES:
      fprintf(outf, "INHERIT_HALF_OF_ALL_EDGES");
      break;
    case INHERIT_ALL_EDGES:
      fprintf(outf, "INHERIT_ALL_EDGES");
      break;
  }
  fputs(";\n", outf);
  fprintf(outf, "w->edge_weights=%s;\n", edge_weights_string(w->edge_weights));
  fprintf(outf, "w->multi_edges=%s;\n", w->multi_edges ? "true" : "false");
  fprintf(outf, "w->allow_move_edge=%s;\n", w->allow_move_edge ? "true" : "false");
  fputc('\n', outf);
  fputs("// Overall run parameters\n", outf);
  fprintf(outf, "w->num_epochs=%d;\n", w->num_epochs);
  fprintf(outf, "w->generations_per_epoch=%d;\n", w->generations_per_epoch);
  fprintf(outf, "w->num_organisms=%d;\n", w->num_organisms);
  fprintf(outf, "w->num_candidates=%d;\n", w->num_candidates);
  fprintf(outf, "w->num_nodes=%d;\n", w->num_nodes);
  fprintf(outf, "w->num_edges=%d;\n", w->num_edges);
  fputc('\n', outf);
  fprintf(outf, "w->num_hill_climbers=%d;\n", w->num_hill_climbers);
  fflush(outf);
}

void print_acclivity_measures_of_best(World *w) {
  int best_organism_index = find_best_organism(w);
  //print_phenotype(w->organisms[best_organism_index]->genotype);
  Organism *o = w->organisms[best_organism_index];
  HILL_CLIMBING_RESULT gvector_result =
    genotype_acclivity(w, o->genotype); //measure_acclivity(w, o);
  printf("gvector acclivity: fitness_delta = %lf, absolute_fitness = %lf, coverage = %lf\n",
    gvector_result.fitness_delta, gvector_result.ending_fitness,
    measure_coverage(w, o->genotype));
  HILL_CLIMBING_RESULT phenotype_result = phenotype_acclivity(w);
  printf("phenotype acclivity: fitness_delta = %lf, absolute_fitness = %lf\n",
    phenotype_result.fitness_delta, phenotype_result.ending_fitness);
}

void print_knob_fitness_numbers(World *w) {
  printf("pos_knob_turns = %lf\n", w->num_fitness_increases_from_knob_turn /
      (double) w->num_generations_measured);
}

void run_world(World *w) {
  printf("--------------------------------------------------------------------------------\n");
  w->c1 = 0.5;
  change_fitness_constants(w);
  set_ridge_coords(w);
  print_world_params(w, stdout);

  log_preamble(w);

  srand(w->seed);
  init_random_population(w);
  for (int e = 1; e <= w->num_epochs; e++) {
    run_epoch(w, e);
  }
  if (!quiet) {
    dump_virtual_fitness_func(w, true, stdout);
    dump_phenotype_fitness_func(w, true, stdout);
  }
  printf("epoch fitness deltas: ");
  print_stats(w->epoch_fitness_deltas);
  //print_data(w->epoch_fitness_deltas);
  //print_acclivity_measures_of_best(w);
  print_knob_fitness_numbers(w);

  close_ancestor_log(w);
}

// -- fitness ----------------------------------------------------------------

double many_small_hills(const double *phenotype) { // length is 2
  return cos(phenotype[0] * 30.0) * sin(phenotype[1] * 30.0);
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double invv(World *w, double target, double radius, double x) {
  double dist = fabs(target - x);
  if (dist >= radius) {
    return 0.0;
  } else {
    if (w->invu) {
        return 1 - pow(dist / radius, 2);
    } else {
        return (radius - dist) / radius;
    }
  }
}

double along_ridge(World *w, double x, double y) {
  //printf("x = %lf; y = %lf; fabs(%lf) = %lf\n", x, y, y - (w->c2 * x + w->c3), fabs(y - (w->c2 * x + w->c3)));
  switch (w->ridge_type) {
  case LINE:
    return invv(w, 0.0, w->ridge_radius, fabs(y - (w->c2 * x + w->c3)));
  case CIRCLE:
    return invv(w, 0.0, w->ridge_radius, fabs((x*x + y*y) - (0.5*0.5)));
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


// phenotype -> array of g->num_out doubles
void fill_phenotype_from_genotype(Genotype *g, double *phenotype) {
  for (int i = 0; i < g->num_out; i++)
    phenotype[i] = g->nodes[g->num_in + i].final_output;
}

// phenotype -> array of 2 doubles
double phenotype_fitness(World *w, const double *phenotype) {
  if (verbose >= 2) {
    printf("require_valid_region(w, %lf, %lf) = %lf\n", phenotype[0], phenotype[1], require_valid_region(w, phenotype[0], phenotype[1]));
    printf("along_ridge(%lf, %lf) = %lf\n", phenotype[0], phenotype[1], along_ridge(w, phenotype[0], phenotype[1]));
  }
  double fitness = 0.0;
  if (phenotype[0] != UNWRITTEN && phenotype[1] != UNWRITTEN) {
    double dist = distance(w->peak_x, w->peak_y, phenotype[0], phenotype[1]);
    double scaled_dist = (w->max_dist - dist) / w->max_dist;
    fitness = require_valid_region(w, phenotype[0], phenotype[1]) *
              along_ridge(w, phenotype[0], phenotype[1]) *
              (w->distance_weight * pow(scaled_dist, w->dist_exponent));
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

double coverage_reward(World *w, Genotype *g) {
  if (w->reward_coverage)
    return 2.0 * measure_coverage(w, g);
  else
    return 0.0;
}

double genotype_fitness(World *w, Genotype *g) {
  double phenotype[g->num_out];
  fill_phenotype_from_genotype(g, phenotype);
  return phenotype_fitness(w, phenotype); // DCB + coverage_reward(w, g);
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

void set_edge_weight(Genotype *g, int edge_index, double weight) {
  if (edge_index >= g->num_edges)
    return;
  g->edges[edge_index].weight = weight;
}

void add_edge(World *w, Genotype *g, int src, int dst, double weight) {
  if (!has_node(g, src) || src < 0 || !has_node(g, dst) || dst < 0)
    return;
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
  int src = select_in_use_node(g);
  int dst = select_in_use_node(g);
  double weight = rand_edge_weight(w);
  assert(has_node(g, src));
  assert(src >= 0);
  assert(has_node(g, dst));
  assert(dst >= 0);
  add_edge(w, g, src, dst, weight);
  MUTATION_RECORD rec;
  EDGE_MUTATION_RECORD edge_rec = { g->num_edges - 1, src, dst, weight };
  rec.type = MUT_ADD_EDGE;
  rec.edge_mutation = edge_rec;
  set_mutation(w, o, &rec);
}

void remove_edge(Genotype *g, int e) {
  if (e < g->num_edges - 1)
    memmove(&g->edges[e], &g->edges[e + 1],
        sizeof(Edge) * (g->num_edges - e - 1));
  g->num_edges--;
}

void mut_remove_edge(World *w, Organism *o) {
  Genotype *g = o->genotype;
  EDGE_MUTATION_RECORD edge_rec = { -1, -1, -1, 0.0 }; // default to "null" mutation
  if (g->num_edges > 0) {
    int selected_edge = rand() % g->num_edges;
    Edge *e = &g->edges[selected_edge];

    edge_rec.index = selected_edge;
    edge_rec.src = e->src;
    edge_rec.dst = e->dst;
    edge_rec.weight = e->weight;

    remove_edge(g, selected_edge);
  }
  MUTATION_RECORD rec;
  rec.type = MUT_REMOVE_EDGE;
  rec.edge_mutation = edge_rec;
  set_mutation(w, o, &rec);
}

void mut_move_edge(World *w, Organism *o) {
  Genotype *g = o->genotype;
  if (g->num_edges == 0)
    return mut_add_edge(w, o);
  int selected_edge = rand() % g->num_edges;
  Edge *e = &g->edges[selected_edge];
  MOVE_EDGE_MUTATION_RECORD edge_rec = { selected_edge, e->src, -1, e->dst, -1, e->weight };
  if (rand() & 1)
    e->src = select_in_use_node(g);
  else
    e->dst = select_in_use_node(g);
  edge_rec.src = e->src;
  edge_rec.dst = e->dst;
  MUTATION_RECORD rec;
  rec.type = MUT_MOVE_EDGE;
  rec.move_edge_mutation = edge_rec;
  set_mutation(w, o, &rec);
}

int add_node(World *w, Genotype *g) {
  int add_index = take_first_unused_node(g);

  g->num_nodes_in_use++;
  if (add_index == -1) {
    g->num_nodes++;
    g->nodes = realloc(g->nodes, sizeof(Node) * g->num_nodes);
    add_index = g->num_nodes - 1;
  }
  Node *node = &g->nodes[add_index];
  init_random_node(w, node, add_index);
  return add_index;
}

void mut_add_node(World *w, Organism *o) {
  Genotype *g = o->genotype;
  int add_index = add_node(w, g);
  Node *node = &g->nodes[add_index];
  NODE_MUTATION_RECORD node_rec = {
    add_index,
    node->initial_activation,
    node->control,
    node->input_acc,
    node->activation_type,
    node->output_type,
    node->control_update,
    node->initial_activation_type
  };
  MUTATION_RECORD rec;
  rec.type = MUT_ADD_NODE;
  rec.node_mutation = node_rec;
  set_mutation(w, o, &rec);

  if (debug)
    sanity_check_organism(w, o);
}

void remove_node(Genotype *g, int selected_node) {
  int unremovable = g->num_in + g->num_out;
  if (selected_node < unremovable || selected_node >= g->num_nodes_in_use)
    return;

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

void mut_remove_node(World *w, Organism *o) {
  Genotype *g = o->genotype;
  int unremovable = g->num_in + g->num_out;
  if (g->num_nodes_in_use <= unremovable)
    return;
  int selected_node = select_in_use_removable_node(g);
  if (selected_node == -1)
    return;
  assert(selected_node >= g->num_in + g->num_out);

  Node *n = &g->nodes[selected_node];
  NODE_MUTATION_RECORD node_rec = {
    selected_node,
    n->initial_activation,
    n->control,
    n->input_acc,
    n->activation_type,
    n->output_type,
    n->control_update,
    n->initial_activation_type
  };
  remove_node(g, selected_node);
  MUTATION_RECORD rec;
  rec.type = MUT_REMOVE_NODE;
  rec.node_mutation = node_rec;
  set_mutation(w, o, &rec);
}

void set_knob(Genotype *g, int genotype_index, double value) {
  if (genotype_index >= g->num_in)
    return;
  Node *node_to_change = &g->nodes[genotype_index];
  node_to_change->initial_activation = clamp(value);
}

void mut_turn_knob(World *w, Organism *o) {
  int genotype_index = rand() % o->genotype->num_in;
  Node *node_to_change = &o->genotype->nodes[genotype_index];
  double nudge;
  KNOB_TURN_RECORD knob_rec = {
    genotype_index,
    node_to_change->initial_activation,
    0.0,
    0.0
  };
  switch (w->knob_type) {
    case KNOB_DISCRETE:
      nudge = coin_flip() ? w->knob_constant : -w->knob_constant;
      break;
    case KNOB_NORMAL:
      nudge = sample_normal(w->knob_constant);
      break;
  }
  knob_rec.nudge = nudge;
  node_to_change->initial_activation =
      clamp(node_to_change->initial_activation + nudge);
  o->from_turned_knob = true; // FIXME: switch to use mutation info instead
  knob_rec.value = node_to_change->initial_activation;
  MUTATION_RECORD rec;
  rec.type = MUT_TURN_KNOB;
  rec.knob_turn = knob_rec;
  set_mutation(w, o, &rec);
}

void set_activation_type(Genotype *g, int node_index, ACTIVATION_TYPE type) {
  if (node_index >= g->num_nodes)
    return;
  Node *n = &g->nodes[node_index];
  n->activation_type = type;
}

void mut_alter_activation_type(World *w, Organism *o) {
  switch (w->activation_types) {
    case 0x01: // SIGMOID only
    case 0x02: // CLAMP_ONLY only
      mut_turn_knob(w, o);
      break;
    case 0x03: // SIGMOID or CLAMP_ONLY
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_ACTIVATION_TYPE_RECORD alter_rec = {
          n,
          node_to_change->activation_type,
          -1
        };
        switch (node_to_change->activation_type) {
          case SIGMOID:
            node_to_change->activation_type = CLAMP_ONLY;
            break;
          case CLAMP_ONLY:
            node_to_change->activation_type = SIGMOID;
            break;
          default:
            assert(false);
        }
        alter_rec.type = node_to_change->activation_type;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_ACTIVATION_TYPE;
        rec.alter_activation_type = alter_rec;
        set_mutation(w, o, &rec);
      }
      break;
    default:
      assert(false);
  }
}

void set_input_acc(Genotype *g, int node_index, INPUT_ACC type) {
  if (node_index >= g->num_nodes)
    return;
  Node *n = &g->nodes[node_index];
  n->input_acc = type;
}

void mut_alter_input_acc(World *w, Organism *o) {
  switch (w->input_accs) {
    case 0x01: // SUM_INCOMING only
    case 0x02: // MULT_INCOMING only
    case 0x04: // MIN_INCOMING only
      mut_turn_knob(w, o);
      break;
    case 0x03: // SUM_INCOMING or MULT_INCOMING
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_INPUT_ACC_RECORD alter_rec = {
          n,
          node_to_change->input_acc,
          -1
        };
        switch (node_to_change->input_acc) {
          case SUM_INCOMING:
            node_to_change->input_acc = MULT_INCOMING;
            break;
          case MULT_INCOMING:
            node_to_change->input_acc = SUM_INCOMING;
            break;
          default:
            assert(false);
        }
        alter_rec.acc = node_to_change->input_acc;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_INPUT_ACC;
        rec.alter_input_acc = alter_rec;
        set_mutation(w, o, &rec);
      }
      break;
    case 0x05: // SUM_INCOMING or MIN_INCOMING
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_INPUT_ACC_RECORD alter_rec = {
          n,
          node_to_change->input_acc,
          -1
        };
        switch (node_to_change->input_acc) {
          case SUM_INCOMING:
            node_to_change->input_acc = MIN_INCOMING;
            break;
          case MIN_INCOMING:
            node_to_change->input_acc = SUM_INCOMING;
            break;
          default:
            assert(false);
        }
        alter_rec.acc = node_to_change->input_acc;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_INPUT_ACC;
        rec.alter_input_acc = alter_rec;
        set_mutation(w, o, &rec);
      }
      break;
    case 0x07: // SUM_INCOMING, MULT_INCOMING, or MIN_INCOMING
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_INPUT_ACC_RECORD alter_rec = {
          n,
          node_to_change->input_acc,
          -1
        };
        switch (node_to_change->input_acc) {
          case SUM_INCOMING:
            node_to_change->input_acc =
                coin_flip() ? MULT_INCOMING : MIN_INCOMING;
            break;
          case MULT_INCOMING:
            node_to_change->input_acc =
                coin_flip() ? SUM_INCOMING : MIN_INCOMING;
            break;
          case MIN_INCOMING:
            node_to_change->input_acc =
                coin_flip() ? SUM_INCOMING : MULT_INCOMING;
            break;
          default:
            assert(false);
        }
        alter_rec.acc = node_to_change->input_acc;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_INPUT_ACC;
        rec.alter_input_acc = alter_rec;
        set_mutation(w, o, &rec);
      }
      break;
    default:
      assert(false);
  }
}

void set_output_type(Genotype *g, int node_index, OUTPUT_TYPE type) {
  if (node_index >= g->num_nodes)
    return;
  Node *n = &g->nodes[node_index];
  n->output_type = type;
}

void mut_alter_output_type(World *w, Organism *o) {
  switch (w->output_types) {
    case ONLY_PASS_THROUGH:
      break;
    case PASS_THROUGH_AND_STEEP_SIGMOID:
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_OUTPUT_TYPE_RECORD alter_rec = {
          n,
          node_to_change->output_type,
          -1
        };
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
        alter_rec.type = node_to_change->output_type;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_OUTPUT_TYPE;
        rec.alter_output_type = alter_rec;
        set_mutation(w, o, &rec);
      }
      break;
    case PASS_THROUGH_AND_TWO_STEP:
      {
        int n = select_in_use_node(o->genotype);
        Node *node_to_change = &o->genotype->nodes[n];
        ALTER_OUTPUT_TYPE_RECORD alter_rec = {
          n,
          node_to_change->output_type,
          -1
        };
        switch (node_to_change->output_type) {
          case PASS_THROUGH:
            node_to_change->output_type = TWO_STEP;
            break;
          case TWO_STEP:
            node_to_change->output_type = PASS_THROUGH;
            break;
          default:
            assert(false);
        }
        alter_rec.type = node_to_change->output_type;
        MUTATION_RECORD rec;
        rec.type = MUT_ALTER_OUTPUT_TYPE;
        rec.alter_output_type = alter_rec;
        set_mutation(w, o, &rec);
      }
  }
}

void set_control(Genotype *g, int selected, double value) {
  if (selected >= g->num_nodes_in_use)
    return;
  Node *node = &g->nodes[selected];
  node->control = clamp(value);
}

void mut_turn_control(World *w, Organism *o) {
  Genotype *g = o->genotype;
  int selected = select_in_use_node(g);
  Node *node = &g->nodes[selected];
  KNOB_TURN_RECORD knob_rec = {
    selected,
    node->control,
    0.0,
    0.0
  };
  double nudge = sample_normal(w->control_increment);
  knob_rec.nudge = nudge;
  node->control = clamp(node->control + nudge);
  knob_rec.value = node->control;
  MUTATION_RECORD rec;
  rec.type = MUT_TURN_CONTROL;
  rec.knob_turn = knob_rec;
  set_mutation(w, o, &rec);
}

Organism *mutate(World *w, int parent, Organism *old_o) {
  Organism *o = copy_organism(old_o);
  int num_mutations =
      1 + (int)(w->extra_mutation_rate *
                rand_float() *
                (o->genotype->num_nodes + o->genotype->num_edges));
  set_mutation_start(w, o, parent, num_mutations);
  for (int i = 0; i < num_mutations; i++) {
    int mutation_type = rand_int(0, w->mutation_type_ub);
    switch (mutation_type) {
    case 0:
      mut_add_node(w, o);
      break;
    case 1:
      mut_remove_node(w, o);
      break;
    case 2:
      mut_add_edge(w, o);
      break;
    case 3:
      mut_remove_edge(w, o);
      break;
    case 4:
      mut_alter_activation_type(w, o);
      break;
    case 5:
      mut_alter_input_acc(w, o);
      break;
    case 6:
      mut_alter_output_type(w, o);
      break;
    case 7:
      mut_turn_control(w, o);
      break;
    case 8:
      if (w->allow_move_edge) {
        mut_move_edge(w, o);
        break;
      }
      // else fall through
    default:
      mut_turn_knob(w, o);
      break;
    }
  }
  if (debug)
    sanity_check_organism(w, o);
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

Organism *crossover(World *w, int mom_idx, Organism *m, int dad_idx, Organism *d) {
  Genotype *mommy = m->genotype;
  Genotype *daddy = d->genotype;
  Organism *baby_o = calloc(1, sizeof(Organism));
  Genotype *baby = calloc(1, sizeof(Genotype));
  baby_o->genotype = baby;
  baby_o->fitness = 0.0;
  baby_o->from_turned_knob = false;

  assert(mommy->num_in == daddy->num_in);
  assert(mommy->num_out == daddy->num_out);
  baby->num_in = mommy->num_in;
  baby->num_out = mommy->num_out;

  double crossover_frac = rand_float();
  int mommy_crossover_point = mommy->num_nodes * crossover_frac;
  int daddy_crossover_point = daddy->num_nodes * crossover_frac;

  set_crossover_info(w, baby_o, mom_idx, dad_idx, crossover_frac);

  int num_from_mommy = mommy_crossover_point;
  int num_from_daddy = daddy->num_nodes - daddy_crossover_point;
  baby->num_nodes = num_from_mommy + num_from_daddy;
  baby->nodes = malloc(sizeof(Node) * baby->num_nodes);
  assert(baby->num_nodes >= mommy->num_in + mommy->num_out);
  int n = 0;
  for (int m = 0; m < mommy_crossover_point; m++) {
    baby->nodes[n++] = mommy->nodes[m];
  }
  for (int d = daddy_crossover_point; d < daddy->num_nodes; d++) {
    baby->nodes[n] = daddy->nodes[d];
    if (n < baby->num_in && baby->nodes[n].initial_activation == UNWRITTEN) {
      baby->nodes[n].initial_activation = mommy->nodes[n].initial_activation;
    }
    n++;
  }
  assert(n == baby->num_nodes);

  int in_use = 0;
  for (int n = 0; n < baby->num_nodes; n++) {
    if (n < baby->num_in + baby->num_out)
      baby->nodes[n].in_use = true;
    if (n >= baby->num_in && n < baby->num_in + baby->num_out)
      baby->nodes[n].initial_activation = UNWRITTEN;
    if (baby->nodes[n].in_use)
      in_use++;
  }
  baby->num_nodes_in_use = in_use;

  baby->edges = NULL;
  baby->num_edges = 0;
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
    }
  }
  for (int d = 0; d < daddy->num_edges; d++) {
    Edge *edge = &daddy->edges[d];
    int bsrc = (edge->src - daddy_crossover_point) + mommy_crossover_point;
    int bdst = (edge->dst - daddy_crossover_point) + mommy_crossover_point;
    switch (w->edge_inheritance) {
      case NO_EDGES_ACROSS_PARENTS:
        if (edge->src >= daddy_crossover_point && edge->dst >= daddy_crossover_point) {
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
    }
  }
  if (debug)
    sanity_check_organism(w, baby_o);
  return baby_o;
}

// ---------------------------------------------------------------------------

void sa_test() {
  Node nodes[6] = {
    { true, 0.2, 0.0, 0.0, 0.0 },
    { true, 0.4, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 }
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
  //Organism o = { &genotype, 0.0 };
  
  verbose = 9;
  World *w = create_world();
  w->sa_timesteps = 13;
  w->decay = 1.0;
  w->spreading_rate = 1.0;
  //sa(&genotype, 13, 1.0, 1.0);
  sa(w, &genotype, stdout);
}

void sa_test2() {
  Node nodes[] = {
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, 0.0, 0.0, 0.0, 0.0 },
    { true, -1.0, 0.0, 0.0, 0.0 }
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
  //Organism o = { &genotype, 0.0 };
  
  verbose = 9;
  World *w = create_world();
  w->sa_timesteps = 20;
  w->decay = 1.0;
  w->spreading_rate = 0.01;
  //sa(&genotype, 20, 1.0, 0.01);
  sa(w, &genotype, stdout);
}

void dump_phenotype_fitness() {
  World *w = create_world();
  w->num_organisms = 40;
  w->ridge_radius=0.200000;
  w->c2=1.000000; w->c3=0.000000;
  print_world_params(w, stdout);
  dump_phenotype_fitness_func(w, true, stdout);
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
  w->activation_types = 0x01;
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
  w->activation_types=0x02; // CLAMP_ONLY only

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
  w->activation_types=0x02; // CLAMP_ONLY only

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
  w->activation_types = 0x02; // CLAMP_ONLY only
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
  run_world(w);
  int best_organism_index = find_best_organism(w);
  print_phenotype(w->organisms[best_organism_index]->genotype);
  HILL_CLIMBING_RESULT gvector_result = genotype_acclivity(w, w->organisms[best_organism_index]->genotype); //measure_acclivity(w, w->organisms[best_organism_index]);
  printf("gvector: average delta = %lf, average fitness = %lf\n", gvector_result.fitness_delta, gvector_result.ending_fitness);
  HILL_CLIMBING_RESULT phenotype_result = phenotype_acclivity(w);
  printf("phenotype: average delta = %lf, average fitness = %lf\n", phenotype_result.fitness_delta, phenotype_result.ending_fitness);
}

void run_from_command_line_options(int argc, char **argv) {
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
    { "dist_exponent", required_argument, 0, 0 },
    { "input_accs", required_argument, 0, 0 },
    { "control_increment", required_argument, 0, 0 },
    { "run", required_argument, 0, 0 },
    { "param_set", required_argument, 0, 0 },
    { "log", required_argument, 0, 0 },
    { "reward_coverage", required_argument, 0, 0 },
    { "invu", required_argument, 0, 0 },
    { "noquiet", no_argument, 0, 0 },
    { NULL, 0, 0, 0 },
  };
  int c;

  for ( ; ; ) {
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
      case 34:
        w->dist_exponent = atof(optarg);
        break;
      case 35:
        w->input_accs = atoi(optarg);
        break;
      case 36:
        w->control_increment = atof(optarg);
        break;
      case 37:
        w->run = atoi(optarg);
        break;
      case 38:
        w->param_set = atoi(optarg);
        break;
      case 39:
        w->log = fopen(optarg, "wb");
        break;
      case 40:
        w->reward_coverage = atoi(optarg);
        break;
      case 41:
        w->invu = atoi(optarg);
        break;
      case 42:
        quiet = false;
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

void mtest() {
  measure_coverage(create_world(), NULL);
}

#ifndef WITH_SWIG
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Wunused-variable\"")
int main(int argc, char **argv) {
  int seed = get_seed(argv, argc);
  run_from_command_line_options(argc, argv);
  //mtest();
  //sa_test();
  //sa_test2();
  //quick_test(seed);
  //dot_test(seed);
  //long_test(seed);
  //long_test_start_small(seed);  //(677953487); // the main test
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
#endif
