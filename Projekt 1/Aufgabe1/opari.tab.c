#include "pomp_lib.h"


extern struct ompregdescr omp_rd_10;

int POMP_MAX_ID = 11;

struct ompregdescr* pomp_rd_table[11] = {
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  &omp_rd_10,
};
