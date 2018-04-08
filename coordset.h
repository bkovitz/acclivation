#ifndef COORDSET_INCLUDED
#define COORDSET_INCLUDED

#include <stdbool.h>

typedef struct {
  int x, y;  // rounded coordinates
  bool is_covered;
} COORD;

typedef void *COORD_SET;

COORD *make_coord(int x, int y);
void add_coord(COORD_SET *csetp, int x, int y);
void cover_coord(COORD_SET cset, int x, int y);
void clear_coverage(COORD_SET cset);
bool has_coord(COORD_SET cset, int x, int y);
bool is_covered(COORD_SET cset, int x, int y);

#endif // COORDSET_INCLUDED
