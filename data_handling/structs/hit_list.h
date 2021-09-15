#include <stdlib.h>

#ifndef HIT_LIST
#define HIT_LIST

struct hit_list {
  size_t* evtid; // maybe not necessary
  double* r;
  double* phi;
  double* z;
  int* layer;
  // pretty sure we don't need index for c++ version
  size_t* index; // what actually is index -- both type-wise and practically speaking
  size_t* pos;
  size_t size;
  size_t* pid;
  int* volume_id; // only used before data transformed to add layers
  int* layer_id; // ''

};

typedef struct hit_list* hit_list_T;

// create new hit list, initializing *most* of the array attributes to empty size-long arrays.
hit_list_T create_hit_list(size_t size);

// free hit list and all array attributes.
void free_hit_list(hit_list_T hit_list);



#endif
