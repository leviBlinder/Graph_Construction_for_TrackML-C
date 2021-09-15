#include "hit_list.h"


#ifndef SEGMENT_H
#define SEGMENT_H

struct segment {
  size_t size;
  hit_list_T hits1;
  hit_list_T hits2;
  double* dr;
  double* dphi;
  double* dz;
  double* dR;
};


typedef struct segment* segment_T;
// create new segment, creating 2 hit_lists and initializing array attributes to empty size-long arrays.
segment_T create_segment(size_t size);

// free segment, hit_lists, and all array attributes.
void free_segment(segment_T segment);

#endif
