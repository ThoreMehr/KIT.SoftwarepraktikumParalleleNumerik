#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallelfor", "", 0, "mandelbrot.c", 46, 46, 60, 60
};

#define POMP_DLIST_00001 shared(omp_rd_1)

