#include "coordset.h"

#include <search.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


COORD *make_coord(int x, int y) {
  COORD *result = (COORD *)malloc(sizeof(COORD));
  result->x = x;
  result->y = y;
  result->is_covered = false;
  return result;
}

void print_coord(const COORD *c) {
  printf("(%d, %d, %s)\n", c->x, c->y, c->is_covered ? "true" : "false");
}

static int compare_coords(const void *pc1_, const void *pc2_) {
  const COORD *pc1 = pc1_, *pc2 = pc2_;

  if (pc1->x < pc2->x)
    return -1;
  if (pc1->x > pc2->x)
    return 1;
  if (pc1->y < pc2->y)
    return -1;
  if (pc1->y > pc2->y)
    return 1;
  return 0;
}

// Does nothing if coord is not in cset.
void cover_coord(COORD_SET cset, int x, int y) {
  COORD c = { x, y, false };
  void *v = tfind(&c, &cset, compare_coords);
  if (v) {
    COORD *pc = *(COORD **)v;
    pc->is_covered = true;
  }
}

void add_coord(COORD_SET *csetp, int x, int y) {
  COORD *pc = make_coord(x, y);
  void *v = tsearch(pc, csetp, compare_coords);
  if (v == NULL) {
    perror("add_coord");
    exit(-1);
  } else if (*(COORD **)v != pc) {
    // cset already has these coordinates
    free(pc);
  }
}

static void clear_coverage_action(const void *nodep, VISIT which, int depth) {
  COORD *c = *(COORD **)nodep;

  switch (which) {
    case preorder:
    case endorder:
      break;
    case postorder:
    case leaf:
      c->is_covered = false;
      break;
  }
}

void clear_coverage(COORD_SET cset) {
  twalk(cset, clear_coverage_action);
}

bool has_coord(COORD_SET cset, int x, int y) {
  COORD c = { x, y, false };
  return tfind(&c, &cset, compare_coords);
}

bool is_covered(COORD_SET cset, int x, int y) {
  COORD c = { x, y, false };
  void *v = tfind(&c, &cset, compare_coords);
  if (v) {
    COORD *pc = *(COORD **)v;
    return pc->is_covered;
  } else {
    return false;
  }
}

int oldmain(int argc, char **argv) {
  COORD *pc;
  void *v;
  void *root = NULL;

  pc = make_coord(2, 3);
  v = tsearch(pc, &root, compare_coords);
  if (v == NULL) {
    perror("tsearch");
    exit(-1);
  } else if (*(COORD **)v != pc) {
    puts("already have it");
    free(pc);
  }

  pc = make_coord(2, 4);
  if (tfind(pc, &root, compare_coords))
    puts("found");
  else
    puts("not found");
  free(pc);

  //tdestroy(root, free);
  return 0;
}

const char *boolstr(bool b) {
  return b ? "true" : "false";
}

int testmain(int argc, char **argv) {
  COORD_SET cset = NULL;

  add_coord(&cset, 2, 3);
  printf("ic(2, 3)=%s\n", boolstr(is_covered(cset, 2, 3)));
  cover_coord(cset, 2, 3);
  printf("ic(2, 3)=%s\n", boolstr(is_covered(cset, 2, 3)));
  clear_coverage(cset);
  printf("ic(2, 3)=%s\n", boolstr(is_covered(cset, 2, 3)));
  return 0;
}
