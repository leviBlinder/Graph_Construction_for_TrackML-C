#include <assert.h>
#include "segment.h"
#include "hit_list.h"
#include <stdlib.h>

using namespace std;
typedef struct segment* segment_T;

segment_T create_segment(size_t size){
  segment_T new_seg = (segment_T) malloc(sizeof(segment));
  assert(new_seg != NULL);
  new_seg->hits1 = create_hit_list(size);
  new_seg->hits2 = create_hit_list(size);
  new_seg->dr = (double*) malloc(size * sizeof(double));
  new_seg->dphi = (double*) malloc(size * sizeof(double));
  new_seg->dz = (double*) malloc(size * sizeof(double));
  new_seg->dR = (double*) malloc(size * sizeof(double));
  new_seg->size = size;
  assert(new_seg->dr != NULL);
  assert(new_seg->dphi != NULL);
  assert(new_seg->dz != NULL);
  assert(new_seg->dR != NULL);
  return(new_seg);
}

void free_segment(segment_T segment){
  free_hit_list(segment->hits1);
  free_hit_list(segment->hits2);
  free(segment->dr);
  free(segment->dphi);
  free(segment->dz);
  free(segment->dR);
  free(segment);
}
