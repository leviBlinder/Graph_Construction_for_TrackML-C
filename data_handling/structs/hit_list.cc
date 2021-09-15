#include "hit_list.h"
#include <stdlib.h> // pulls in declaration of malloc, free, ?size_t?
#include <assert.h> // declare assert

using namespace std;
typedef struct hit_list* hit_list_T;

hit_list_T create_hit_list(size_t size){
  hit_list_T new_hits = (hit_list_T) malloc(sizeof(hit_list));
  assert(new_hits != NULL);
  new_hits->evtid = (size_t*) malloc(size * sizeof(size_t));
  new_hits->r = (double*) malloc(size * sizeof(double));
  new_hits->phi = (double*) malloc(size * sizeof(double));
  new_hits->z = (double*) malloc(size * sizeof(double));
  new_hits->layer = (int*) malloc(size * sizeof(int));
  new_hits->pos = (size_t*) malloc(size * sizeof(size_t));
  //Could move this out to segment creation since:
  // 1) not all hit_lists use it and it doesn't need to be mantained as with other vars
  // 2) Its a bit in the weeds for loops, might rather they be in the actual code rather than in struct definitions
  for(size_t i = 0; i < size; i ++){
    new_hits->pos[i] = i;
  }
  new_hits->index = (size_t*) malloc(size * sizeof(size_t));
  new_hits->size = size;
  new_hits->pid = (size_t*) malloc(size*sizeof(size_t*));
  new_hits->layer_id = (int*) malloc(size*sizeof(int*));
  new_hits->volume_id = (int*) malloc(size*sizeof(int*));
  assert(new_hits->evtid != NULL);
  assert(new_hits->r != NULL);
  assert(new_hits->phi != NULL);
  assert(new_hits->z != NULL);
  assert(new_hits->layer != NULL);
  assert(new_hits->index != NULL);
  assert(new_hits->pos != NULL);
  assert(new_hits->pid != NULL);
  assert(new_hits->layer_id != NULL);
  assert(new_hits->volume_id != NULL);
  return(new_hits);
}
void free_hit_list(hit_list_T hit_list){
  free(hit_list->evtid);
  free(hit_list->r);
  free(hit_list->phi);
  free(hit_list->z);
  free(hit_list->layer);
  free(hit_list->index);
  free(hit_list->pos);
  free(hit_list->pid);
  free(hit_list->layer_id);
  free(hit_list->volume_id);
  free(hit_list);
}
