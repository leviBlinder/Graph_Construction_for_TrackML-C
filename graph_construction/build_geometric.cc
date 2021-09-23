/*
Data preparation script for GNN tracking.
This script processes the TrackML dataset and produces graph data on disk.
*/


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

// #include "../data_handling/hep_data.h"
// #include "../data_handling/structs/hit_list.h"
// #include "../data_handling/structs/segment.h"
#include "../data_handling/structs/args.h"
#include "geometric_functions.h"


using namespace std;

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


// main can be run cpu/gpu-side?
void process_event(char* prefix, char* output_dir, geometric_functions* G){
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

  G->max_edges_by_hits = 0;
  for(size_t i = 0; i < psize; i++){
    double pt = sqrt(pow(px[i],2) + pow(py[i],2));
    if(pt >= G->pt_min){
      G->max_edges_by_hits += p_hits[i];
      if(p_hits[i] != 0) G->max_edges_by_hits--;
    }
  }
  cout << "max_edges_by_hits: " << G->max_edges_by_hits << endl;
  // Apply hit selection
  // logging.info('Event %i, selecting hits' % evtid)
  G->set_vlids();
  hits = G->select_hits(hits, pid, px, py, psize);

  cout << "hit list after select size: " << (hits)-> size << endl;

  // Divide detector into sections

  // phi_range = (-np.pi, np.pi)
  double* phi_edges = G->linspace(G->phi_range, G->n_phi_sections);
  double* eta_edges = G->linspace(G->eta_range, G->n_eta_sections);
  hit_list_T* hits_sections = G->split_detector_sections(hits, phi_edges, eta_edges, G->n_phi_sections, G->n_eta_sections);
  size_t n_sections = G->n_phi_sections * G->n_eta_sections;
  cout << "n_graphs: " <<  n_sections << endl;
  // Graph features and scale
  // feature_names = ['r', 'phi', 'z']
  double feature_scale[3] = {1000., PI / G->n_phi_sections, 1000.};
  if(G->phi_reflect) feature_scale[1] *= -1;

  // define adjacent layer pairs
  G->set_layer_pairs();


  // Construct the graph
  // logging.info('Event %i, constructing graphs' % evtid)
  void*** graphs = (void***) malloc(n_sections * sizeof(void*));
  size_t* n_edges_array = (size_t*) malloc(n_sections * sizeof(size_t));
  size_t* n_hits_array = (size_t*) malloc(n_sections * sizeof(size_t));
  for(size_t i = 0; i < n_sections; i++){
    graphs[i] = (void**) G->construct_graph(hits_sections+i, feature_scale, n_edges_array+i);
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

int main(int argc, const char** argv){
  /* Main function */

  // Parse the command line
  args_T args = new struct args();
  args->parse_args(argc, argv);
  char prefix[] = {'t','e','s','t','G','r','a','p','h'};
  char output_dir[] = {'.','.','/','g','r','a','p','h','s','/','t','e','s','t','G','r','a','p','h'};
  double pt_min = 1;
  size_t n_eta_sections = 1;
  size_t n_phi_sections = 1;
  pair<double,double> eta_range (-5,5);
  pair<double,double> phi_range (-PI, PI);
  double phi_slope_max = 0.008;
  double z0_max = 500;
  bool phi_reflect = false;
  bool endcaps = true;
  bool remove_intersecting_edges = true;
  geometric_functions* GEVT = new geometric_functions(pt_min, n_eta_sections, n_phi_sections, eta_range, phi_range, phi_slope_max, z0_max, phi_reflect, endcaps, remove_intersecting_edges);
  process_event(prefix, output_dir, GEVT);


  return 0;
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
