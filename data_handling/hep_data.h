#include<iostream>
#include<fstream>
#include<cmath>
#include"hit_list.h"
using namespace std;
#ifndef HEP_DATA
#define HEP_DATA
// There may be a better but probably slower way to read files on GormAnalysis "Reading And Writing CSV Files With C++"


// void read_hits_limited(char* filename, hit_list_T hits);
// void read_truth_limited(char* filename, hit_list_T hits);
// void read_particles_limited(size_t size, char* filename, void*** pparticles);
void read_limited(size_t evtid, hit_list_T* phits, void*** pparticles);
void write_hits(char* filename, hit_list_T hits);

// void free_graph_partial(void** graph, size_t n_hits);
// void free_sparse(void* sparse);
void** graph_to_sparse(void** graph, size_t n_hits, size_t n_edges);

// void free_sparse_partial(void* sparse);
void** sparse_to_graph(void** sparse, size_t n_hits, size_t n_edges);

void save_graphs(void*** graphs, char** filenames, size_t* n_hits_array, size_t* n_edges_array, size_t n_files);
void save_graph(void** graph, char* filename, size_t n_hits, size_t n_edges);
#endif
