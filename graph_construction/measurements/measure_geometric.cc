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


//#include "../../data_handling/hep_data.h"
//#include "../../data_handling/structs/hit_list.h"
//#include "../../data_handling/structs/segment.h"
#include "../../data_handling/structs/args.h"
#include "../geometric_functions.h"


using namespace std;

void process_event(geometric_functions* G){
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


  // Construct the graph, don't save; only constructing to find n_edges, n_nodes, other features.
  // logging.info('Event %i, constructing graphs' % evtid)
  for(size_t i = 0; i < n_sections; i++){
    (void) G->construct_graph(hits_sections+i, feature_scale, NULL);
  }
}

void process_max_event(geometric_functions* G){
  size_t evtid = 1000;
  hit_list_T hits;
  void** particles;
  read_limited(evtid, &hits, &particles);
  size_t* pid = (size_t*) particles[0];
  double* px = (double*) particles[1];
  double* py = (double*) particles[2];
  size_t* ppsize = (size_t*) particles[3];
  size_t psize = *ppsize;
  G->set_vlids();
  hits = G->select_hits(hits, pid, px, py, psize);
  double feature_scale[3] = {1000., PI , 1000.};
  if(G->phi_reflect) feature_scale[1] *= -1;
  size_t* pn_edges =  (size_t*) malloc(sizeof(size_t));
  G->soft_construct_graph(&hits, feature_scale, pn_edges);
}



int main(int argc, const char** argv){
  /* Main function */

  // Parse the command line
  args_T args = new struct args();
  args->parse_args(argc, argv);
  ofstream exportStatsFile("statsCompatible.csv", ios::app);
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
    geometric_functions* GEVT = new geometric_functions(pt_min, n_eta_sections, n_phi_sections, eta_range, phi_range, phi_slope_max, z0_max, phi_reflect, endcaps, remove_intersecting_edges);
    process_event(GEVT);
    process_max_event(GEVT);
    // stats
    GEVT->efficiency_by_hits = ((double) GEVT->num_true) / GEVT->max_edges_by_hits;
    GEVT->efficiency_by_construction = ((double) GEVT->num_true) / GEVT->max_edges_by_construction;
    cout << "pt_min: " << GEVT->pt_min << " phi_slope_max: " << GEVT->phi_slope_max << " z0_max: " << GEVT->z0_max << endl;
    cout << "num_nodes: " << GEVT->num_nodes << " num_edges: " << GEVT->num_edges << " num_true: " << GEVT->num_true << " purity: " << GEVT->purity << endl;
    cout << "max_edges_by_hits: " << GEVT->max_edges_by_hits << " max_edges_by_construction: " << GEVT->max_edges_by_construction << " efficiency_by_hits: "
    << GEVT->efficiency_by_hits << " efficiency_by_construction: " << GEVT->efficiency_by_construction << endl;
    //ofstream statsFile("../stats.csv", ios::app);
    //statsFile << phi_slope_max << "," << z0_max << "," << pt_min << "," << num_nodes << "," << num_edges << "," << num_true << "," << purity << "," << max_edges_by_construction << "," << efficiency_by_construction << "," << max_edges_by_hits << "," << efficiency_by_hits << "\n";
    double purity_err = 0;
    double efficiency_err = 0;
    double n_edges_err = 0;
    double n_nodes_err = 0;
    exportStatsFile << pt_min << "," << GEVT->purity << "," << purity_err << "," << GEVT->efficiency_by_construction << "," << efficiency_err << "," << GEVT->num_edges << "," << n_edges_err << "," << GEVT->num_nodes << "," << n_nodes_err << "\n";
  }
  exportStatsFile.close();
  return 0;
}
