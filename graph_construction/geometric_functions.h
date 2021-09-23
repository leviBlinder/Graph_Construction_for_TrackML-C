/*
Data preparation script for GNN tracking.
This script processes the TrackML dataset and produces graph data on disk.
*/
#ifndef BUILD_GEOMETRIC_H
#define BUILD_GEOMETRIC_H



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
#include "../data_handling/structs/segment.h"


using namespace std;

#define PI 3.141592653589793238462643383279502884

template<typename T>
inline void freeContainer(T& p_container)
{
    T empty;
    using std::swap;
    swap(p_container, empty);
}
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

class geometric_functions{
public:
  size_t max_layer;
  int** vlids;
  int** layer_pairs;
  size_t n_layer_pairs;

  size_t num_nodes;
  size_t num_edges;
  size_t num_true;
  double purity;
  size_t max_edges_by_hits;
  size_t max_edges_by_construction;
  double efficiency_by_hits;
  double efficiency_by_construction;

  double pt_min;
  size_t n_eta_sections;
  size_t n_phi_sections;
  pair<double,double> eta_range;
  pair<double,double> phi_range;
  double phi_slope_max;
  double z0_max;
  bool phi_reflect;
  bool endcaps;
  bool remove_intersecting_edges;

  geometric_functions(double _pt_min, size_t _n_eta_sections, size_t _n_phi_sections, pair<double,double> _eta_range, pair<double,double> _phi_range, double _phi_slope_max, double _z0_max, bool _phi_reflect, bool _endcaps, bool _remove_intersecting_edges);

  /*args_T create_args();

  void parse_arg_bool(int argc, const char** argv, string keyword, bool* var);

  void parse_arg_int(int argc, const char** argv, string keyword, int* var);

  void parse_args(int argc, const char** argv, args_T args); */

  double* calc_eta(double* r, double* z, size_t length);

  segment_T create_paired_segment(hit_list_T hits1, hit_list_T hits2);

  segment_T select_segment(segment_T seg, double phi_slope_cut, double z0_cut, bool intersect_cut);

  void** construct_graph(hit_list_T* phits, double* feature_scale, size_t* pn_edges);


  int vlid_hash(int vlid_volume, int vlid_layer);


  void set_vlids();

  // Define adjacent layers
  void set_layer_pairs();


  hit_list_T select_hits(hit_list_T hits, size_t* pid, double* px, double* py, size_t psize);


  hit_list_T* split_detector_sections(hit_list_T hits, double* phi_edges, double* eta_edges, size_t n_phi_edges, size_t n_eta_edges);


  double* linspace(pair<double,double> range, size_t n_sections);

  void soft_construct_graph(hit_list_T* phits, double* feature_scale, size_t* pn_edges);

};

#endif
