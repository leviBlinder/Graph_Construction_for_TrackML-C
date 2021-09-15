/*
Data preparation script for GNN tracking.
This script processes the TrackML dataset and produces graph data on disk.
*/

#include "segment.h"
#include <utility>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <chrono> // for clock
#include <unordered_map>
#include <unordered_set>
#include <fstream>

#include "../data_handling/hep_data.h"
#include "../data_handling/structs/hit_list.h"


using namespace std;

#define PI 3.141592653589793238462643383279502884

struct args {
  char* config = NULL;
  int n_workers = 1;
  int task = 0;
  int n_tasks = 1;
  bool verbose = false;
  bool show_config = false;
  bool interactive = false;
  int start_evtid = 1000;
  int end_evtid = 2700;
};
typedef struct args* args_T;

template<typename T>
inline void freeContainer(T& p_container)
{
    T empty;
    using std::swap;
    swap(p_container, empty);
}

int max_layer;
int** vlids = (int**) malloc(2 * sizeof(void*));
int** layer_pairs = (int**) malloc(2 * sizeof(void*));
int n_layer_pairs;

size_t num_nodes = 0;
size_t num_edges = 0;
size_t num_true = 0;
double purity = 0;
size_t max_edges_by_hits = 0;
size_t max_edges_by_construction = 0;
double efficiency_by_hits = 0;
double efficiency_by_construction = 0;

// Need to figure out what includes to use for c++

/*
// System
#include os;
#include sys
#include time;
#include argparse;
#include logging;
#include multiprocessing as mp;
#include functools::partial;
sys.path.append("../");

// Externals
#include yaml;
#include pickle;
#include numpy as np;
#include pandas as pd;
#include trackml.dataset;
#include time;


// Locals
#include models::graph::Graph;
#include models::graph::save_graphs;
*/

args_T create_args(){
  args_T newArgs = (args_T) malloc(sizeof(args));
  return newArgs;
}

void parse_arg_bool(int argc, const char** argv, string keyword, bool* var){
  for(int i = 1; i < argc; i++){
    string argvi_string(argv[i]);
    if(argvi_string.compare(keyword) == 0){
      *var = true;
      return;
    }
  }
  return;
}

void parse_arg_int(int argc, const char** argv, string keyword, int* var){
  int i;
  for(i = 1; i < argc; i++){
    string argvi_string(argv[i]);
    if(argvi_string.compare(keyword) == 0){
      assert(argv[i+1] != NULL);
      *var = atoi(argv[i+1]); // might need to add "nullptr" argument
      return;
    }
  }
  return;
}

void parse_args(int argc, const char** argv, args_T args){
  /* Parse command line arguments. */
  assert(argv != NULL);
  //could just assign 1st arg to config, but not sure about order
  //parse_arg_int(argc, argv, "--n-workers", &(args->config));
  parse_arg_int(argc, argv, "--n-workers", &(args->n_workers));
  parse_arg_int(argc, argv, "--task", &(args->task));
  parse_arg_int(argc, argv, "--n-tasks", &(args->n_tasks));
  parse_arg_bool(argc, argv, "-v", &(args->verbose));
  parse_arg_bool(argc, argv, "--show-config", &(args->show_config));
  parse_arg_bool(argc, argv, "--interactive", &(args->interactive));
  parse_arg_int(argc, argv, "--start-evtid", &(args->start_evtid));
  parse_arg_int(argc, argv, "--end-evtid", &(args->end_evtid));
  return;
}


double* calc_eta(double* r, double* z, size_t length){
  double* eta = (double*) malloc(length * sizeof(double));
  for(size_t i = 0; i < length; i++){
    eta[i] = atan2(r[i], z[i]);
    eta[i] = -1. * log(tan(eta[i]/2));
  }
  return eta;
}

segment_T create_paired_segment(hit_list_T hits1, hit_list_T hits2){
  //cout << "start create seg" << endl;
  /*
  Construct all possible edges from the pairings
  between hits1 and hits2, which each contain all hits from one layer from one evtid.
  */
  size_t i = 0;
  segment_T newSeg = create_segment(hits1->size * hits2->size);

  for(size_t i1 = 0; i1 < hits1->size; i1++){
    for(size_t i2 = 0; i2 < hits2->size; i2++){
      newSeg->hits1->r[i] = hits1->r[i1];
      newSeg->hits1->phi[i] = hits1->phi[i1];
      newSeg->hits1->z[i] = hits1->z[i1];
      newSeg->hits1->layer[i] = hits1->layer[i1];
      newSeg->hits1->pos[i] = hits1->pos[i1];
      newSeg->hits1->pid[i] = hits1->pid[i1];
      newSeg->hits2->r[i] = hits2->r[i2];
      newSeg->hits2->phi[i] = hits2->phi[i2];
      newSeg->hits2->z[i] = hits2->z[i2];
      newSeg->hits2->layer[i] = hits2->layer[i2];
      newSeg->hits2->pos[i] = hits2->pos[i2];
      newSeg->hits2->pid[i] = hits2->pid[i2];
      newSeg->dz[i] = hits2->z[i2] - hits1->z[i1];
      newSeg->dr[i] = hits2->r[i2] - hits1->r[i1];
      i++;
    }
  }
  free_hit_list(hits1);
  free_hit_list(hits2);
  //cout << "mid create seg" << endl;

  for(i = 0; i < newSeg -> size; i++){
    double eta1 = atan2(newSeg->hits1->r[i], newSeg->hits1->z[i]);
    eta1 = -1. * log(tan(eta1/2));
    double eta2 = atan2(newSeg->hits2->r[i], newSeg->hits2->z[i]);
    eta2 = -1. * log(tan(eta2/2));
    double deta = eta2 - eta1;
    double dphi = newSeg->hits2->phi[i] - newSeg->hits1->phi[i];
    if(dphi > PI) dphi -= 2 * PI;
    if(dphi < -PI) dphi += 2 * PI;
    newSeg->dphi[i] = dphi;
    newSeg->dR[i] = sqrt((deta *  deta) + (dphi * dphi)); // depency on prev line, not prev loop iterations tho so no problem.
  }

  //cout << "end create seg" << endl;
  return newSeg;
}

segment_T select_segment(segment_T seg, double phi_slope_max, double z0_max, bool remove_intersecting_edges = false){
  //cout << "start select seg" << endl;
  /*
  Construct a list of selected segments from the pairings
  between hits1 and hits2, filtered with the specified
  phi slope and z0 criteria.
  Returns: pd DataFrame of (index_1, index_2), corresponding to the
  DataFrame hit label-indices in hits1 and hits2, respectively.
  */
  int layer1 = seg->hits1->layer[0];
  int layer2 = seg->hits2->layer[0];
  // Access all possible pairs of hits by double-indexing(note that there should only be 1
  // value for evtid for any input dataset)

  // Compute line through the points
  const size_t dvLength = seg->size * sizeof(double);
  // calc_* functions might be better in the main loop below for efficiency in unrolling, not sure.
 // WTF?????

  double* phi_slope = (double*) malloc(dvLength);
  double* z0 = (double*) malloc(dvLength);
  // This seems like it may not unroll well since it needs access to initial data input (from hits1/2) for so long
  // Consider optimizing based on OpenCL best practices manual.
  size_t index = 0;
  for(size_t i = 0; i < seg->size; i++){
    //cout << "z: " << (seg->hits1->z)[i] << "r: " << (seg->hits1->r)[i] << " dz: " << (seg->dz[i]) << " dr: " << (seg->dr[i]) << endl;
    z0[i] = (seg->hits1->z)[i] - ((seg->hits1->r)[i] * (seg->dz[i])/(seg->dr[i]));
    phi_slope[i] = (seg->dphi[i])/(seg->dr[i]);
    index++;
  }


  // Could add below vector operation to above loop as well. Haven't considered efficiency.

  // print("max(z0) = ", np.max(abs(z0)))
  // Apply the intersecting line cut
  bool* intersected_layer = (bool*) malloc(seg->size * sizeof(bool));
  for(size_t i = 0; i < seg->size; i++){
    // What is .abs()? If it were absolute value, why not just assign false. Can't find any documentation.
    intersected_layer[i] = abs(seg->dr[i]) < -1;
    if(intersected_layer[i] != false) cout << "WTF" << endl;
  }

  if(remove_intersecting_edges){
    // Innermost barrel layer --> innermost L,R endcap layers
    if( (layer1 == 0) && (layer2 == 11 or layer2 == 4) ){
      double* z_coord = (double*) malloc(dvLength);
      for(size_t i = 0; i < seg->size; i++){
        z_coord[i] = (71.56298065185547 * (seg->dz[i])/(seg->dr[i])) + z0[i];
        intersected_layer[i] = ((z_coord[i] > -490.975) && (z_coord[i] < 490.975));
      }
      free(z_coord);
    }
    if((layer1 == 1) && (layer2 == 11 or layer2 == 4)){
      double* z_coord = (double*) malloc(dvLength);
      for(size_t i = 0; i < seg->size; i++){
        z_coord[i] = (115.37811279296875 * (seg->dz[i])/(seg->dr[i])) + z0[i];
        intersected_layer[i] = ((z_coord[i] > -490.975) && (z_coord[i] < 490.975));
      }
      free(z_coord);
    }
  }

  // Filter segments according to criteria

  size_t newLength = 0;

  // Find number of cases after filter
  size_t phiCount = 0;
  size_t zCount = 0;
  size_t intersectCount = 0;
  for(size_t i = 0; i < seg->size; i++){
    if((abs(phi_slope[i]) < phi_slope_max) && (abs(z0[i]) < z0_max) && (intersected_layer[i] == false)){
      newLength++;
    }else{
      if(!(abs(phi_slope[i]) < phi_slope_max)){
         phiCount++;
      }else{
        if(!(abs(z0[i]) < z0_max)){
           zCount++;
         }else{
           if(intersected_layer[i] != false) intersectCount++;
         }
       }
    }
  }
  //cout << "layer1: " << layer1 << " layer2: "<< layer2 << endl;
  cout << "seg_size: " << seg->size << " newLength: " << newLength << endl;
  //cout << " phi_cut: " << phiCount << " z_cut: " << zCount << " intersect_cut: " << intersectCount << endl;
  //cout << "mid select seg" << endl;
  // Create segment to hold case data after filter, copy over all data
  segment_T newSeg = create_segment(newLength);
  size_t j = 0;
  //cout << "seg size: " << seg->size << endl;
  for(size_t i = 0; i < seg->size; i++){
    if((abs(phi_slope[i]) < phi_slope_max) && (abs(z0[i]) < z0_max) && (intersected_layer[i] == false)){
      newSeg->hits1->r[j] = seg->hits1->r[i];
      newSeg->hits1->phi[j] = seg->hits1->phi[i];
      newSeg->hits1->z[j] = seg->hits1->z[i];
      newSeg->hits1->layer[j] = seg->hits1->layer[i];
      newSeg->hits1->pos[j] = seg->hits1->pos[i];
      newSeg->hits1->pid[j] = seg->hits1->pid[i];
      newSeg->hits2->r[j] = seg->hits2->r[i];
      newSeg->hits2->phi[j] = seg->hits2->phi[i];
      newSeg->hits2->z[j] = seg->hits2->z[i];
      newSeg->hits2->layer[j] = seg->hits2->layer[i];
      newSeg->hits2->pos[j] = seg->hits2->pos[i];
      newSeg->hits2->pid[j] = seg->hits2->pid[i];
      newSeg->dr[j] = seg->dr[i];
      newSeg->dphi[j] = seg->dphi[i];
      newSeg->dz[j] = seg->dz[i];
      newSeg->dR[j] = seg->dR[i];
      j++;
    }
  }
  // Free old segment data
  // as well as arrays used for filter
  free_segment(seg);

  free(intersected_layer);
  free(z0);
  free(phi_slope);
  //cout << "end select seg" << endl;
  return newSeg;
}

void** construct_graph(hit_list_T* phits, double phi_slope_max, double z0_max, double* feature_scale, size_t* pn_edges, bool remove_intersecting_edges = false){
/* Construct one graph (e.g. from one event) */
  hit_list_T hits = *phits;
  chrono::time_point<std::chrono::system_clock> t0 = chrono::system_clock::now();

  // Loop over layer pairs and construct segments

  // Equivalent to hits.group_by(layers) -----------------------------------------------------
  size_t n_hits = hits->size;
  size_t** lgroup_indices = (size_t**) malloc(sizeof(int*) * max_layer);
  for(size_t i = 0; i < max_layer; i++){
    lgroup_indices[i] = (size_t*) calloc(n_hits, sizeof(int));
  }
  int* lgroup_sizes = (int*) calloc(max_layer,sizeof(int));
  for(size_t i = 0; i < n_hits; i++){
    int lindex = hits->layer[i];
    if(lindex >= max_layer){
      cout << "lindex: " << lindex << endl;
    }
    lgroup_indices[lindex][lgroup_sizes[lindex]] = i;
    lgroup_sizes[lindex]++;
  }
  // ------------------------------------------------------------------------------------------

  segment_T* seperate_segments = (segment_T*) calloc(n_layer_pairs,sizeof(segment_T));
  size_t totalSegmentSize = 0;
  for(size_t i = 0; i < n_layer_pairs; i++) {
    // If an event has no hits on a layer, just skip to the next layer pair

    int layer1 = layer_pairs[0][i];
    int layer2 = layer_pairs[1][i];
    if(lgroup_sizes[layer1] == 0){
      /*printf("skipping empty layer: %d\n", layer1);*/
      continue;
    }
    if(lgroup_sizes[layer2] == 0){
      /*printf("skipping empty layer: %d\n", layer2);*/
      continue;
    }

    // Find and join all hit pairs
    hit_list_T hits1 = create_hit_list(lgroup_sizes[layer1]);
    size_t* layer_group_1 = lgroup_indices[layer1];
    for(size_t i = 0; i < lgroup_sizes[layer1]; i++){
      size_t layer_index = layer_group_1[i];
      // hits1->evtid[i] = hits->evtid[layer_index];
      hits1->r[i] = hits->r[layer_index];
      hits1->phi[i] = hits->phi[layer_index];
      hits1->z[i] = hits->z[layer_index];
      hits1->layer[i] = hits->layer[layer_index]; //should just be layer1
      hits1->pid[i] = hits->pid[layer_index];
    }
    hit_list_T hits2 = create_hit_list(lgroup_sizes[layer2]);
    size_t* layer_group_2 = lgroup_indices[layer2];
    for(size_t i = 0; i < lgroup_sizes[layer2]; i++){
      size_t layer_index = layer_group_2[i];
      // hits2->evtid[i] = hits->evtid[layer_index];
      hits2->r[i] = hits->r[layer_index];
      hits2->phi[i] = hits->phi[layer_index];
      hits2->z[i] = hits->z[layer_index];
      hits2->layer[i] = hits->layer[layer_index]; //should just be layer2
      hits2->pid[i] = hits->pid[layer_index];
    }

    // Construct the segments
    segment_T allEdgeSegment = create_paired_segment(hits1, hits2);
    segment_T selectedEdgeSegment = select_segment(allEdgeSegment, phi_slope_max, z0_max, remove_intersecting_edges);



    // Combine segments from all layer pairs
    seperate_segments[i] = selectedEdgeSegment;
    totalSegmentSize += selectedEdgeSegment->size;
  }
  segment_T segments = create_segment(totalSegmentSize);
  size_t count = 0;
  for(size_t i = 0; i < n_layer_pairs; i++){
    if(seperate_segments[i] == NULL) continue;
    for(size_t j = 0; j < seperate_segments[i]->size; j++){
      segments->hits1->r[count] = seperate_segments[i]->hits1->r[j];
      segments->hits1->phi[count] = seperate_segments[i]->hits1->phi[j];
      segments->hits1->z[count] = seperate_segments[i]->hits1->z[j];
      segments->hits1->layer[count] = seperate_segments[i]->hits1->layer[j];
      segments->hits1->pos[count] = seperate_segments[i]->hits1->pos[j];
      segments->hits1->pid[count] = seperate_segments[i]->hits1->pid[j];
      segments->hits2->r[count] = seperate_segments[i]->hits2->r[j];
      segments->hits2->phi[count] = seperate_segments[i]->hits2->phi[j];
      segments->hits2->z[count] = seperate_segments[i]->hits2->z[j];
      segments->hits2->layer[count] = seperate_segments[i]->hits2->layer[j];
      segments->hits2->pos[count] = seperate_segments[i]->hits2->pos[j];
      segments->hits2->pid[count] = seperate_segments[i]->hits2->pid[j];
      segments->dr[count] = seperate_segments[i]->dr[j];
      segments->dphi[count] = seperate_segments[i]->dphi[j];
      segments->dz[count] = seperate_segments[i]->dz[j];
      segments->dR[count] = seperate_segments[i]->dR[j];
      count++;
    }
  }
  for(size_t i = 0; i < n_layer_pairs; i++){
    if(seperate_segments[i] == NULL) continue;
    free_segment(seperate_segments[i]);
  }
  free(seperate_segments);

  // Prepare the graph matrices
  // size_t n_hits = hits->size; already defined
  size_t n_edges = segments->size;
  cout << "n_edges: " << n_edges << endl;
  *pn_edges = n_edges;

  // Not using the feature_names = ['r','phi','z'] for copying into X since it
  // creates double loop, and hard to get hits->r using the value of a variable to get string literal 'r'.
  size_t n_features = 3;
  double** x = (double**) malloc(n_features * sizeof(double*));
  double** Ra = (double**) malloc(4 * sizeof(double*));
  assert(x != NULL);
  assert(Ra != NULL);
  for(size_t i = 0; i < n_features; i++){
    x[i] = (double*) malloc(n_hits * sizeof(double));
    assert(x[i] != NULL);
  }
  for(size_t i = 0; i < 4; i++){
    Ra[i] = (double*) malloc(n_edges * sizeof(double));
    assert(Ra[i] != NULL);
  }
  for(size_t i = 0; i < n_hits; i++){
    x[0][i] = hits->r[i]/feature_scale[0];
    x[1][i] = hits->phi[i]/feature_scale[1];
    x[2][i] = hits->z[i]/feature_scale[2];
  }

  for(size_t i = 0; i < n_edges; i++){
    Ra[0][i] = segments->dr[i]/feature_scale[0];
    Ra[1][i] = segments->dphi[i]/feature_scale[1];
    Ra[2][i] = segments->dz[i]/feature_scale[2];
    Ra[3][i] = segments->dR[i];
  }


  /* Ra = np.stack((seg_dr/feature_scale[0],
                 seg_dphi/feature_scale[1],
                 seg_dz/feature_scale[2],
                 seg_dR))
                 */
  // Ra = np.zeros(n_edges)
  // I'm confused by ^. Ra seems to be a 4xn_hits matrix based on np.stack, not n_edges x 1.

  int** edge_index = (int**) malloc(2 * sizeof(int*)); // could use char instead of int?
  assert(edge_index != NULL);


  edge_index[0] = (int*) calloc(n_edges, sizeof(int));
  edge_index[1] = (int*) calloc(n_edges, sizeof(int));
  assert(edge_index[0] != NULL);
  assert(edge_index[1] != NULL);
  bool* y = (bool*) calloc(n_edges, sizeof(double));
  assert(y != NULL);

  cout << "n_hits: " << n_hits << endl;

  // We have the segments' hits given by dataframe label,
  // so we need to translate into positional indices.

  // Now we can fill the association matrices.

  // instead of having a pos attribute for hit_list structs, could just use
  // a counter that resets when hits1->layer[i] or hits2->layer[i] changes, since
  // the segments are in order.

  for(size_t i = 0; i < n_edges; i++){
    size_t seg_start = lgroup_indices[segments->hits1->layer[i]][segments->hits1->pos[i]];
    size_t seg_end = lgroup_indices[segments->hits2->layer[i]][segments->hits2->pos[i]];
    edge_index[0][i] = seg_start;
    edge_index[1][i] = seg_end;
  }


  // Fill the segment, particle labels

  //What is the purpose of this section??
  /*
  pid = hits.particle_id
  unique_pid_map = {pid_old: pid_new
                    for pid_new, pid_old in enumerate(np.unique(pid.values))}
  pid_mapped = pid.map(unique_pid_map)
  print(pid_mapped)
  */


  //double* pid1 = (double*) malloc(n_edges * sizeof(double));
  //double* pid2 = (double*) malloc(n_edges * sizeof(double));
  num_true = 0;
  for(size_t i = 0; i < n_edges; i++){

    //pid1[i] = segments->hits1->pid[i];
    //pid2[i] = segments->hits2->pid[i];

    y[i] = (segments->hits1->pid[i] == segments->hits2->pid[i]);
    if(y[i] == true){
      num_true++;
      //cout << "i: " << i << " y[i]: " << y[i] << " pid1: " << segments->hits1->pid[i] << " pid2: " << segments->hits2->pid[i] << endl;
    }
  }
  purity = ((double) num_true)/(n_edges);
  num_nodes = n_hits;
  num_edges = n_edges;
  free_segment(segments);

  // Correct for multiple true barrel-endcap segments
  // Does very little, ignoring for now
  /*
  layer1 = hits.layer.loc[segments.index_1].values
  layer2 = hits.layer.loc[segments.index_2].values
  true_layer1 = layer1[y>0.5]
  true_layer2 = layer2[y>0.5]
  true_pid1 = pid1[y>0.5]
  true_pid2 = pid2[y>0.5]
  true_z1 = hits.z.loc[segments.index_1].values[y>0.5]
  true_z2 = hits.z.loc[segments.index_2].values[y>0.5]
  for(pid in np.unique(true_pid1)):
      temp = [0, 0, -1, -1]
      l1, l2 = true_layer1[true_pid1==pid], true_layer2[true_pid2==pid]
      z1, z2 = true_z1[true_pid1==pid], true_z2[true_pid2==pid]
      for l in range(len(l1)):
          barrel_to_EC = (l1[l] in [0,1,2,3] and l2[l] in [4,11])
          if (barrel_to_EC):
              temp[0] += 1
              if abs(temp[1]) < abs(z1[l]): // why????
                  if (temp[0] > 1):
                      y[(pid1==pid2) & (pid1==pid) &
                        (layer1==temp[2]) & (layer2==temp[3])] = 0

                  temp[1] = abs(z1[l])
                  temp[2] = l1[l]
                  temp[3] = l2[l]
  */
  chrono::time_point<std::chrono::system_clock> t1 = chrono::system_clock::now();
  cout << "took " << (t1-t0).count()/1000000 << " seconds" << endl;

  /*
  print("X.shape", X.shape)
  print("Ri.shape", Ri.shape)
  print("Ro.shape", Ro.shape)
  print("y.shape", y.shape)
  print("Ra.shape", Ra.shape)
  print("pid.shape", pid_mapped.shape)
  */
  void** graph = (void**) malloc(sizeof(void*) * 5);
  graph[0] = x;
  graph[1] = Ra;
  graph[2] = edge_index;
  graph[3] = y;
  graph[4] = (size_t*) malloc(sizeof(size_t)* n_hits);
  for(int i = 0; i < n_hits; i++){
    ((size_t*)graph[4])[i] = hits->pid[i];
  }
  free_hit_list(hits);
  return graph;
}





// now also simplifies pid list
// instead of particles as argument, use:
// double* px, double* py, double* pid
// and other columns as well?
int vlid_hash(int vlid_volume, int vlid_layer){
  return vlid_layer*10 + vlid_volume; //10 is greater that max value of vlid_volume
}

void set_vlids(bool endcaps){
  static int vlidsEnd_volume[18] = {8,8,8,8,7,7,7,7,7,7,7,9,9,9,9,9,9,9};
  static int vlidsEnd_layer[18] = {2,4,6,8,14,12,10,8,6,4,2,2,4,6,8,10,12,14};

  static int vlidsBase_volume[18] = {8,8,8,8};
  static int vlidsBase_layer[18] = {2,4,6,8};

  if(endcaps){
    max_layer = 18;
    vlids[0] = vlidsEnd_volume;
    vlids[1] = vlidsEnd_layer;
  }else{
    max_layer = 4;
    vlids[0] = vlidsBase_volume;
    vlids[1] = vlidsBase_layer;
  }
}

// Define adjacent layers
void set_layer_pairs(bool endcaps){
  static int layers1[3] = {0,1,2};
  static int layers2[3] = {1,2,3};
  static int layers1END[23] = {0,1,2,4,5,6,7,8,9,11,12,13,14,15,16,0,1,2,3,0,1,2,3};
  static int layers2END[23] = {1,2,3,5,6,7,8,9,10,12,13,14,15,16,17,4,4,4,4,11,11,11,11};
  if(!endcaps){
    n_layer_pairs = 3;
    layer_pairs[0] = layers1;
    layer_pairs[1] = layers2;
  }else{
    n_layer_pairs = 23;
    layer_pairs[0] = layers1END;
    layer_pairs[1] = layers2END;
  }
}



hit_list_T select_hits(hit_list_T hits, size_t* pid, double* px, double* py, size_t psize, double pt_min=0){

  size_t n_hits = hits->size;
  cout << "hit list select 0 size: " << (hits)-> size << endl;

  // Select barrel layers and assign convenient layer number [0-9]


  // Using an arbitrary hash function,
  // assign each vlid pair a single identifying value and
  // map that value to a convinient layer number [0-9]
  unordered_map<int, int> vlid_hashes;
  int* vlids_volume = vlids[0];
  int* vlids_layer = vlids[1];
  for(int i = 0; i < max_layer; i++){
    int hashed_value = vlid_hash(vlids_volume[i], vlids_layer[i]);
    vlid_hashes.insert(make_pair<int,int>(hashed_value,i));
  }

  // Calculate particle transverse momentum
  double* pt = (double*) malloc(psize * sizeof(double));
  for(size_t i = 0; i < psize; i++){
    pt[i] = sqrt(pow(px[i],2) + pow(py[i],2));
  }
  // create pt cut filter for true particle selection and
  // pid mapping to simplify removing duplicate hits.

  // Note that we don't use the particle data again, so no need to filter these
  // arrays (px, py, pt, pid, etc.)


  unordered_map<size_t, int> pid_to_mapped;
  pid_to_mapped.reserve(psize);
  size_t mapping = 1;
  for(size_t i = 0; i < psize; i++){
    if(pt[i] > pt_min){
      pid_to_mapped.insert(make_pair<size_t, int>(pid[i],mapping));
      mapping++;
    }
  }
  //cout << "mapping max: " << mapping << endl;


  // Changes particle id for each hit to its simplified version.
  // Applies pt cut, removes all noise hits.
  // also, remove all duplicate hits (i.e., those with same layer and p_id)
  // and remove all hits with out-of-bounds layer_id and volume_id such that the layer wasn't assigned
  // additionally, calculate derived hits variables
  // and calculate number of remaining hits
  size_t new_index = 0;
  size_t* old_indexes = (size_t*) calloc(hits->size, sizeof(size_t));
  size_t countNoPID_or_lowPT = 0;
  size_t countBadLayer = 0;
  size_t countSeen = 0;
  unordered_set<size_t> pidlayer_seen;
  for(size_t i = 0; i < hits->size; i++){
    if(pid_to_mapped.count(hits->pid[i]) != 0 ){ // Get rid of any hits with invalid pid
      size_t hashed_layer = vlid_hash(hits->volume_id[i], hits->layer_id[i]);
      if(vlid_hashes.count(hashed_layer) != 0){ // Get rid of any hits with layer_id/volume_id out of barrel
        // find convinient layer number [0-9] corresponding to layer_id/volume_id pair
        hits->layer[i] = vlid_hashes.find(hashed_layer)->second;

        // find convinient pid value corresponding to original raw pid value
        size_t pid_mapping = pid_to_mapped.find(hits->pid[i])->second;

        // arbitrary hash new pid and layer values (to avoid duplicate hits)
        size_t hashed_pidlayer = pid_mapping * max_layer + hits->layer[i];

        if(pidlayer_seen.count(hashed_pidlayer) == 0){ // get rid of duplicates
          old_indexes[new_index] = i; // keep track of which hits are valid so they can be accessed later
          hits->pid[i] = pid_mapping; // update pid to simplified value in [0 to n_hits]
          pidlayer_seen.insert(hashed_pidlayer);
          new_index++;
        }else{
          countSeen++;
          //cout << "hash: " << hashed_pidlayer << " pmap: " << pid_mapping << " i: " << i << " pid: " << hits->pid[i] << " layer: " << hits->layer[i] << " lid: " << hits->layer_id[i] << " vid: " << hits->volume_id[i] << endl;
        }
      }else{
        countBadLayer++;
      }
    }else{
      /*if(hits->pid[i] != 0){
        cout << "what!? i: " << i << " pid: " << hits->pid[i];
      }*/
      countNoPID_or_lowPT++;
    }
  }
  cout << "countNoPID_or_lowPT: " << countNoPID_or_lowPT <<
  " countBadLayer: " << countBadLayer <<
  " countSeen: " << countSeen << endl;

  // In order to delete invalid hits,
  // just copy all valid hits into a new hit_list and free old hit_list.
  size_t newLength = new_index;
  hit_list_T newHits = create_hit_list(newLength);
  for(size_t i = 0; i < newLength; i++){
    size_t old_i = old_indexes[i];
    newHits->pid[i] = hits->pid[old_i];
    newHits->r[i] = hits->r[old_i];
    newHits->phi[i] = hits->phi[old_i];
    newHits->z[i] = hits->z[old_i];
    newHits->layer[i] = hits->layer[old_i];
  }
  cout << "hit list select 1 size: " << (newHits)-> size << endl;
  free(old_indexes);
  freeContainer(pidlayer_seen);
  freeContainer(vlid_hashes);
  freeContainer(pid_to_mapped);
  free_hit_list(hits);
  return newHits;
}

hit_list_T* split_detector_sections(hit_list_T hits, double* phi_edges, double* eta_edges, size_t n_phi_edges, size_t n_eta_edges){
    cout << "start split_detector_sections" << endl;
    /* Split hits according to provided phi and eta boundaries. */

    size_t n_hits = hits->size;
    size_t n_sections = n_phi_edges * n_eta_edges;
    size_t** section_indices = (size_t**) malloc(sizeof(void*) * n_sections);
    for(size_t i = 0; i < n_sections; i++){
      section_indices[i] = (size_t*) calloc(n_hits, sizeof(size_t));
    }
    size_t* section_sizes = (size_t*) calloc(n_sections, sizeof(size_t));

    double* temp_eta = calc_eta(hits->r, hits->z, hits->size);
    size_t out_phi = 0;
    size_t out_eta = 0;
    for(size_t i = 0; i < n_hits; i++){
      double phi_cur = hits->phi[i];
      double eta_cur = temp_eta[i];
      // Find phi section of hit
      double phi_max;
      size_t iPhi;
      for(iPhi = 0; iPhi < n_phi_edges+1; iPhi++){
        phi_max = phi_edges[iPhi];
        if(i == 0){
          cout << "iPhi: " << iPhi << " phi_max: " << phi_max << " phi_cur: " << phi_cur << endl;
        }
        if(phi_cur <= phi_max) break;
      }

      double phi_min = phi_edges[iPhi - 1];
      if(i == 0){
        cout << "phi_min: " << phi_min << " phi_max: " << phi_max << endl;
      }
      // Find eta section of hit
      double eta_max;
      size_t iEta;
      for(iEta = 0; iEta < n_eta_edges+1; iEta++){
        eta_max = eta_edges[iEta];
        if(eta_cur <= eta_max) break;
        if(i == 0){
          cout << "iEta: " << iEta << " eta_max: " << eta_max << " eta_cur: " << eta_cur << endl;
        }
      }
      if(i == 0){
        cout << "iEta: " << iEta << endl;
      }

      if(iPhi == 0 || iPhi == n_phi_edges+1){ //get rid of out of bounds phi values
        out_phi++;
        continue;
      }
      if(iEta == 0 || iEta == n_eta_edges+1){ //get rid of out of bounds eta values
        out_eta++;
        continue;
      }

      // Center hit on phi=0
      hits->phi[i] -= (phi_min + phi_max) / 2;

      // Update section size and indices for section corresponding to hit
      size_t sectionIndex = ((iPhi-1) * (n_eta_edges-1)) + (iEta-1);
      section_indices[sectionIndex][section_sizes[sectionIndex]] = i;
      section_sizes[sectionIndex]++;
    }
    free(temp_eta);
    cout << "phi out of bounds num: " << out_phi << endl;
    cout << "eta out of bounds num: " << out_eta << endl;

    hit_list_T* hit_sections = (hit_list_T*) malloc(n_sections * sizeof(hit_list_T));
    for(size_t sid = 0; sid < n_sections; sid++){
      hit_list_T hit_section = create_hit_list(section_sizes[sid]);
      for(size_t i = 0; i < section_sizes[sid]; i++){
        hit_section->pid[i] = hits->pid[section_indices[sid][i]];
        hit_section->r[i] = hits->r[section_indices[sid][i]];
        hit_section->phi[i] =hits->phi[section_indices[sid][i]];
        hit_section->layer[i] = hits->layer[section_indices[sid][i]];
        hit_section->z[i] = hits->z[section_indices[sid][i]];
      }
      hit_sections[sid] = hit_section;
    }

    for(size_t i = 0; i < n_sections; i++){
      free(section_indices[i]);
    }
    free(section_indices);
    free(section_sizes);

    cout << "num_sections: " << n_sections << endl;
    cout << "end split_detector_sections" << endl;
    return hit_sections;
  }

double* linspace(pair<double,double> range, size_t n_sections){
  double interval = (range.second - range.first)/n_sections;
  double cur_place = range.first;
  double* return_edges = (double*) malloc((n_sections + 1) * sizeof(double));
  for(int i = 0; i < n_sections + 1; i++){
    cout << "i: " << i << " cur_place: " << cur_place << endl;
    return_edges[i] = cur_place;
    cur_place += interval;
  }
  return return_edges;
}
// main can be run cpu/gpu-side?
void process_event(char* prefix, char* output_dir, double pt_min, size_t n_eta_sections, size_t n_phi_sections, pair<double,double> eta_range, pair<double,double> phi_range, double phi_slope_max, double z0_max, bool phi_reflect, bool endcaps, bool remove_intersecting_edges){
  // Load the data
  size_t evtid = 1000;
  // logging.info('Event %i, loading data' % evtid)
  hit_list_T hits;
  void** particles;
  read_limited(evtid, &hits, &particles);
  cout << "hit list out of read size: " << (hits)-> size << endl;

  size_t* pid = (size_t*) particles[0];
  double* px = (double*) particles[1];
  double* py = (double*) particles[2];
  size_t* ppsize = (size_t*) particles[3];
  size_t* p_hits = (size_t*) particles[4];
  size_t psize = *ppsize;

  max_edges_by_hits = 0;
  for(size_t i = 0; i < psize; i++){
    double pt = sqrt(pow(px[i],2) + pow(py[i],2));
    if(pt >= pt_min){
      max_edges_by_hits += p_hits[i];
      if(p_hits[i] != 0) max_edges_by_hits--;
    }
  }
  cout << "max_edges_by_hits: " << max_edges_by_hits << endl;
  // Apply hit selection
  // logging.info('Event %i, selecting hits' % evtid)
  set_vlids(endcaps);
  hits = select_hits(hits, pid, px, py, psize, pt_min);

  cout << "hit list after select size: " << (hits)-> size << endl;

  // Divide detector into sections

  // phi_range = (-np.pi, np.pi)
  double* phi_edges = linspace(phi_range, n_phi_sections);
  double* eta_edges = linspace(eta_range, n_eta_sections);
  hit_list_T* hits_sections = split_detector_sections(hits, phi_edges, eta_edges, n_phi_sections, n_eta_sections);
  size_t n_sections = n_phi_sections * n_eta_sections;
  cout << "n_graphs: " <<  n_sections << endl;
  // Graph features and scale
  // feature_names = ['r', 'phi', 'z']
  double feature_scale[3] = {1000., PI / n_phi_sections, 1000.};
  if(phi_reflect) feature_scale[1] *= -1;

  // define adjacent layer pairs
  set_layer_pairs(endcaps);


  // Construct the graph
  // logging.info('Event %i, constructing graphs' % evtid)
  void*** graphs = (void***) malloc(n_sections * sizeof(void*));
  size_t* n_edges_array = (size_t*) malloc(n_sections * sizeof(size_t));
  size_t* n_hits_array = (size_t*) malloc(n_sections * sizeof(size_t));
  for(size_t i = 0; i < n_sections; i++){
    graphs[i] = (void**) construct_graph(hits_sections+i, phi_slope_max, z0_max, feature_scale, n_edges_array+i, remove_intersecting_edges);
    n_hits_array[i] = hits_sections[i]->size;
  }


  // Write these graphs to the output directory
  char** filenames = (char**) malloc(sizeof(void*) * n_sections);
  size_t buf = 100;
  for(size_t i = 0; i < n_sections; i++){
    filenames[i] = (char*) calloc(buf, sizeof(char));
    (void) sprintf(filenames[i], "%s%03lu.csv", output_dir, i);
  }

  //logging.info('Event %i, writing graphs', evtid)
  save_graphs(graphs, filenames, n_hits_array, n_edges_array, n_sections);
  free(n_edges_array);
}



/*hit_list_T soft_select_hits(hit_list_T hits, size_t* pid, size_t psize){
  size_t n_hits = hits->size;
  unordered_map<int, int> vlid_hashes;
  int* vlids_volume = vlids[0];
  int* vlids_layer = vlids[1];
  for(int i = 0; i < max_layer; i++){
    int hashed_value = vlid_hash(vlids_volume[i], vlids_layer[i]);
    vlid_hashes.insert(make_pair<int,int>(hashed_value,i));
  }

  // Remove all duplicate hits (i.e., those with same layer and p_id)
  // and remove all hits with out-of-bounds layer_id and volume_id such that the layer wasn't assigned
  // also calculate number of remaining hits
  size_t new_index = 0;
  size_t* old_indexes = (size_t*) calloc(hits->size, sizeof(size_t));
  unordered_set<size_t> pidlayer_seen;
  for(size_t i = 0; i < hits->size; i++){
    if(hits->pid[i] != 0){ // Get rid of any hits with pid == 0
      size_t hashed_layer = vlid_hash(hits->volume_id[i], hits->layer_id[i]);
      if(vlid_hashes.count(hashed_layer) != 0){ // Get rid of any hits with layer_id/volume_id out of barrel
        hits->layer[i] = vlid_hashes.find(hashed_layer)->second; // find convinient layer number [0-9] corresponding to layer_id/volume_id pair
        if(hits -> layer[i] > max_layer){
          cout << "bad i: " << i << " layer: " << hits->layer[i] << endl;
        }
        size_t hashed_pidlayer = hits->pid[i] * max_layer + hits->layer[i]; // arbitrary hash pid and layer values (to avoid duplicate hits)
        if(pidlayer_seen.count(hashed_pidlayer) == 0){ // get rid of duplicates
          old_indexes[new_index] = i; // keep track of which hits are valid so they can be accessed later
          pidlayer_seen.insert(hashed_pidlayer);
          new_index++;
        }
      }
    }
  }

  size_t newLength = new_index;
  hit_list_T newHits = create_hit_list(newLength);
  for(size_t i = 0; i < newLength; i++){
    size_t old_i = old_indexes[i];
    newHits->pid[i] = hits->pid[old_i];
    newHits->r[i] = hits->r[old_i];
    newHits->phi[i] = hits->phi[old_i];
    newHits->z[i] = hits->z[old_i];
    newHits->layer[i] = hits->layer[old_i];
    if(newHits -> layer[i] > max_layer){
      cout << "bad new i: " << i << " layer: " << newHits->layer[i] << endl;
    }
  }
  free(old_indexes);
  freeContainer(pidlayer_seen);
  freeContainer(vlid_hashes);
  free_hit_list(hits);
  return newHits;
}*/

void soft_construct_graph(hit_list_T* phits, double phi_slope_max, double z0_max, double* feature_scale, size_t* pn_edges, bool remove_intersecting_edges = false){
/* Construct one graph (e.g. from one event) */
  hit_list_T hits = *phits;

  // Loop over layer pairs and construct segments

  // Equivalent to hits.group_by(layers) -----------------------------------------------------
  size_t n_hits = hits->size;
  size_t** lgroup_indices = (size_t**) malloc(sizeof(int*) * max_layer);
  for(size_t i = 0; i < max_layer; i++){
    lgroup_indices[i] = (size_t*) calloc(n_hits, sizeof(int));
  }
  int* lgroup_sizes = (int*) calloc(max_layer,sizeof(int));
  for(size_t i = 0; i < n_hits; i++){
    int lindex = hits->layer[i];
    lgroup_indices[lindex][lgroup_sizes[lindex]] = i;
    lgroup_sizes[lindex]++;
  }
  // ------------------------------------------------------------------------------------------

  segment_T* seperate_segments = (segment_T*) calloc(n_layer_pairs,sizeof(segment_T));
  size_t totalSegmentSize = 0;
  for(size_t i = 0; i < n_layer_pairs; i++) {
    // If an event has no hits on a layer, just skip to the next layer pair

    int layer1 = layer_pairs[0][i];
    int layer2 = layer_pairs[1][i];
    if(lgroup_sizes[layer1] == 0){
      continue;
    }
    if(lgroup_sizes[layer2] == 0){
      continue;
    }

    // Find and join all hit pairs
    hit_list_T hits1 = create_hit_list(lgroup_sizes[layer1]);
    size_t* layer_group_1 = lgroup_indices[layer1];
    for(size_t i = 0; i < lgroup_sizes[layer1]; i++){
      size_t layer_index = layer_group_1[i];
      hits1->r[i] = hits->r[layer_index];
      hits1->phi[i] = hits->phi[layer_index];
      hits1->z[i] = hits->z[layer_index];
      hits1->layer[i] = hits->layer[layer_index]; //should just be layer1
      hits1->pid[i] = hits->pid[layer_index];
    }
    hit_list_T hits2 = create_hit_list(lgroup_sizes[layer2]);
    size_t* layer_group_2 = lgroup_indices[layer2];
    for(size_t i = 0; i < lgroup_sizes[layer2]; i++){
      size_t layer_index = layer_group_2[i];
      hits2->r[i] = hits->r[layer_index];
      hits2->phi[i] = hits->phi[layer_index];
      hits2->z[i] = hits->z[layer_index];
      hits2->layer[i] = hits->layer[layer_index]; //should just be layer2
      hits2->pid[i] = hits->pid[layer_index];
    }

    // Construct the segments
    segment_T allEdgeSegment = create_paired_segment(hits1, hits2);
    segment_T selectedEdgeSegment = select_segment(allEdgeSegment, phi_slope_max, z0_max, remove_intersecting_edges);



    // Combine segments from all layer pairs
    seperate_segments[i] = selectedEdgeSegment;
    totalSegmentSize += selectedEdgeSegment->size;
  }
  segment_T segments = create_segment(totalSegmentSize);
  size_t count = 0;
  for(size_t i = 0; i < n_layer_pairs; i++){
    if(seperate_segments[i] == NULL) continue;
    for(size_t j = 0; j < seperate_segments[i]->size; j++){
      segments->hits1->r[count] = seperate_segments[i]->hits1->r[j];
      segments->hits1->phi[count] = seperate_segments[i]->hits1->phi[j];
      segments->hits1->z[count] = seperate_segments[i]->hits1->z[j];
      segments->hits1->layer[count] = seperate_segments[i]->hits1->layer[j];
      segments->hits1->pos[count] = seperate_segments[i]->hits1->pos[j];
      segments->hits1->pid[count] = seperate_segments[i]->hits1->pid[j];
      segments->hits2->r[count] = seperate_segments[i]->hits2->r[j];
      segments->hits2->phi[count] = seperate_segments[i]->hits2->phi[j];
      segments->hits2->z[count] = seperate_segments[i]->hits2->z[j];
      segments->hits2->layer[count] = seperate_segments[i]->hits2->layer[j];
      segments->hits2->pos[count] = seperate_segments[i]->hits2->pos[j];
      segments->hits2->pid[count] = seperate_segments[i]->hits2->pid[j];
      segments->dr[count] = seperate_segments[i]->dr[j];
      segments->dphi[count] = seperate_segments[i]->dphi[j];
      segments->dz[count] = seperate_segments[i]->dz[j];
      segments->dR[count] = seperate_segments[i]->dR[j];
      count++;
    }
  }
  for(size_t i = 0; i < n_layer_pairs; i++){
    if(seperate_segments[i] == NULL) continue;
    free_segment(seperate_segments[i]);
  }
  free(seperate_segments);

  // Prepare the graph matrices
  // size_t n_hits = hits->size; already defined
  size_t n_edges = segments->size;
  *pn_edges = n_edges;


  // for soft graph construct, we actually only need num true edges, not full graph
  bool* y = (bool*) calloc(n_edges, sizeof(double));
  assert(y != NULL);

  max_edges_by_construction = 0;
  for(size_t i = 0; i < n_edges; i++){
    y[i] = (segments->hits1->pid[i] == segments->hits2->pid[i]);
    if(y[i] == true){
      max_edges_by_construction++;
    }
  }
  free_segment(segments);

}

void process_max_event(bool endcaps, bool phi_reflect, double pt_min){
  size_t evtid = 1000;
  hit_list_T hits;
  void** particles;
  read_limited(evtid, &hits, &particles);
  size_t* pid = (size_t*) particles[0];
  double* px = (double*) particles[1];
  double* py = (double*) particles[2];
  size_t* ppsize = (size_t*) particles[3];
  size_t psize = *ppsize;
  set_vlids(endcaps);
  hits = select_hits(hits, pid, px, py, psize, pt_min);
  double feature_scale[3] = {1000., PI , 1000.};
  if(phi_reflect) feature_scale[1] *= -1;
  size_t* pn_edges =  (size_t*) malloc(sizeof(size_t));
  soft_construct_graph(&hits,2 * PI, 1000000, feature_scale, pn_edges, true);
}



int main(int argc, const char** argv){
  /* Main function */

  // Parse the command line
  args_T args = create_args();
  parse_args(argc, argv, args);
  char prefix[] = {'t','e','s','t','G','r','a','p','h'};
  char output_dir[] = {'.','.','/','g','r','a','p','h','s','/','t','e','s','t','G','r','a','p','h'};
  ofstream exportStatsFile("../statsCompatible.csv", ios::app);
  exportStatsFile << "pt,purity,purity_err,efficiency,efficiency_err,n_edges,n_edges_err,n_nodes,n_nodes_err\n";
  double pt_mins[15] = {2,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1,.9,.8,.7,.6};
  double phi_slopes[15] = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.0062, 0.0064, 0.0066, 0.0068, 0.007, 0.0072, 0.008, 0.0085, 0.0085};
  double z0_maxes[15] = {275, 275, 275, 275, 275, 275, 290, 305, 320, 335, 350, 375, 425, 500, 500};
  for(int graph_index = 0; graph_index < 15; graph_index++){
    double pt_min = pt_mins[graph_index];
    size_t n_eta_sections = 1;
    size_t n_phi_sections = 1;
    pair<double,double> eta_range (-5,5);
    pair<double,double> phi_range (-PI, PI);
    double phi_slope_max = phi_slopes[graph_index];
    double z0_max = z0_maxes[graph_index];
    bool phi_reflect = false;
    bool endcaps = true;
    bool remove_intersecting_edges = true;
    process_event(prefix, output_dir, pt_min, n_eta_sections, n_phi_sections, eta_range, phi_range, phi_slope_max, z0_max, phi_reflect, endcaps, remove_intersecting_edges);
    process_max_event(endcaps, phi_reflect, pt_min);
    // stats
    efficiency_by_hits = ((double) num_true) / max_edges_by_hits;
    efficiency_by_construction = ((double) num_true) / max_edges_by_construction;
    cout << "pt_min: " << pt_min << " phi_slope_max: " << phi_slope_max << " z0_max: " << z0_max << endl;
    cout << "num_nodes: " << num_nodes << " num_edges: " << num_edges << " num_true: " << num_true << " purity: " << purity << endl;
    cout << "max_edges_by_hits: " << max_edges_by_hits << " max_edges_by_construction: " << max_edges_by_construction << " efficiency_by_hits: " << efficiency_by_hits << " efficiency_by_construction: " << efficiency_by_construction << endl;
    //ofstream statsFile("../stats.csv", ios::app);
    //statsFile << phi_slope_max << "," << z0_max << "," << pt_min << "," << num_nodes << "," << num_edges << "," << num_true << "," << purity << "," << max_edges_by_construction << "," << efficiency_by_construction << "," << max_edges_by_hits << "," << efficiency_by_hits << "\n";
    double purity_err = 0;
    double efficiency_err = 0;
    double n_edges_err = 0;
    double n_nodes_err = 0;
    exportStatsFile << pt_min << "," << purity << "," << purity_err << "," << efficiency_by_construction << "," << efficiency_err << "," << num_edges << "," << n_edges_err << "," << num_nodes << "," << n_nodes_err << "\n";
  }
  exportStatsFile.close();
  return 0;
}
