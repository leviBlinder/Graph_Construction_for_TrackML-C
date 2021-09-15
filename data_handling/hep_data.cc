#include<iostream>
#include<fstream>
#include<cmath>
#include <vector>
#include"hep_data.h"
#include"hit_list.h"
#include <stdlib.h>  // malloc
#include <assert.h>
using namespace std;

// The is a better but probably slower way to read files on GormAnalysis "Reading And Writing CSV Files With C++"
void read_hits_limited(char* filename, hit_list_T hits){
  size_t i = 0;
  char line[100];
  ifstream file(filename, ifstream::in);
  if(!file.is_open()) throw runtime_error("Could not open file");
  // skip rid of colnames row
  file.getline(line,100);
  while(file.good()){
    (void) file.getline(line,100,',');
    //hits->hit_id[i] = stoul(line);
    // hit_id only used for merging with truth, which will be done by default here
    (void) file.getline(line,100,',');
    double x = stod(line);
    (void) file.getline(line,100,',');
    double y = stod(line);
    hits->r[i] = sqrt(pow(x,2) + pow(y,2));
    hits->phi[i] = atan2(y,x);
    if(hits->phi[i] > 4){
      cout << "badPhi: " << hits->phi[i] << " x: " << x << " y: " << y << endl;
    }
    (void) file.getline(line,100,',');
    hits->z[i] = stod(line);
    (void) file.getline(line,100,',');
    hits->volume_id[i] = stoi(line);
    (void) file.getline(line,100,',');
    hits->layer_id[i] = stoi(line);
    (void) file.getline(line,100);
    // double module_id = stoi(line);
    if(file.peek() == EOF) break;
    i++;
  }
  file.close();
}
void read_truth_limited(char* filename, hit_list_T hits){
  size_t i = 0;
  char line[100];
  ifstream file(filename, ifstream::in);
  if(!file.is_open()) throw runtime_error("Could not open file");
  // skip rid of colnames row
  file.getline(line,100);
  while(file.good()){
    (void) file.getline(line,100,',');
    //hits->hit_id[i] = stoul(line);
    (void) file.getline(line,100,',');
    hits->pid[i] = stoul(line);
    (void) file.getline(line,100);
    //we don't use other columns (tx,ty,tz,tpx,tpy,tpz,weight)
    if(file.peek() == EOF) break;
    i++;
  }
  file.close();
}
void read_particles_limited(size_t size, char* filename, void*** pparticles){
  size_t i = 0;
  *(pparticles) = (void**) malloc(5 * sizeof(double*));
  assert(*pparticles != NULL);
  (*pparticles)[0] = (size_t*) malloc(size * sizeof(size_t));
  (*pparticles)[1] = (double*) malloc(size * sizeof(double));
  (*pparticles)[2] = (double*) malloc(size * sizeof(double));
  (*pparticles)[3] = (size_t*) malloc(sizeof(size_t));
  (*pparticles)[4] = (size_t*) malloc(size * sizeof(double));
  size_t* particle_id = (size_t*)(*pparticles)[0];
  double* px = (double*)(*pparticles)[1];
  double* py = (double*)(*pparticles)[2];
  size_t* ppsize = (size_t*)(*pparticles)[3];
  size_t* p_hits = (size_t*)(*pparticles)[4];
  assert(particle_id != NULL);
  assert(px != NULL);
  assert(py != NULL);
  assert(ppsize != NULL);
  assert(p_hits != NULL);
  *ppsize = size;
  char line[100];
  ifstream file(filename, ifstream::in);
  if(!file.is_open()) throw runtime_error("Could not open file");
  // skip rid of colnames row
  file.getline(line,100);
  while(file.good()){
    (void) file.getline(line,100,',');
    particle_id[i] = stoul(line);

    //skip vx, vy, vz columns
    (void) file.getline(line,100,',');
    (void) file.getline(line,100,',');
    (void) file.getline(line,100,',');

    (void) file.getline(line,100,',');
    px[i] = stod(line);
    (void) file.getline(line,100,',');
    py[i] = stod(line);
    /* // skip remaining columns (pz, q, nhits)
    (void) file.getline(line,100); */
    (void) file.getline(line,100,',');
    (void) file.getline(line,100,',');
    (void) file.getline(line,100);
    //cout << "nhits line: " << line << endl
    p_hits[i] = stoi(line);
    if(file.peek() == EOF) break;
    i++;
  }
  file.close();

}
void read_limited(size_t evtid, hit_list_T* phits, void*** pparticles){
    char pfilename[100];
    char hfilename[100];
    char tfilename[100];
    char buffer[100];
    size_t psize = 0;
    size_t htsize = 0;
    (void) sprintf(pfilename, "../train_100_events/event00000%zu-particles.csv", evtid);
    ifstream pfile(pfilename, ifstream::in);
    if(!pfile.is_open()) throw runtime_error("Could not open file");
    while(pfile.good()){
      (void) pfile.getline(buffer,100);
      if(pfile.peek() == EOF) break;
      psize++;
    }
    pfile.close();

    (void) sprintf(hfilename, "../train_100_events/event00000%zu-hits.csv", evtid);
    (void) sprintf(tfilename, "../train_100_events/event00000%zu-truth.csv", evtid);
    ifstream htfile(hfilename, ifstream::in);
    if(!htfile.is_open()) throw runtime_error("Could not open file");
    while(htfile.good()){
      (void) htfile.getline(buffer,100);
      if(htfile.peek() == EOF) break;
      htsize++;
    }
    htfile.close();

    cout << "particle size: " << psize << endl;
    cout << "hit/truth size: " << htsize << endl;
    *phits = create_hit_list(htsize);
    cout << "hit list in read size: " << (*phits)-> size << endl;
    read_hits_limited(hfilename, *phits);
    read_truth_limited(tfilename, *phits);
    read_particles_limited(psize, pfilename, pparticles);
  }
void write_hits(char* filename, hit_list_T hits){
    // Make a CSV file with hits columns
      // Create an output filestream object
      std::ofstream file(filename);
      size_t size = hits->size;
      // Send column names to the stream
      file << "evtid,r,phi,z,layer,index,volume_id,layer_id\n";

      // Send data to the stream
      for(int i = 0; i < size; i++)
      {
          file << hits->evtid[i] << ',' << hits->r[i] << ','
          << hits->phi[i] << ',' << hits->z[i] << ','
          << hits->layer[i] << ',' << hits->index[i] << ','
          << hits->volume_id[i] << ',' << hits->layer_id[i] << "\n";
      }

      // Close the file
      file.close();
  }


// maybe not working properly??
void free_graph_partial(void** graph, size_t n_hits){
  for(int i = 2; i <= 3; i++){
    // Ro and Ri are not in sparse representation and have n_hits columns
    void** graphArray = (void**) graph[i];
    for(size_t j = 0; j < n_hits; j++){
      free(graphArray[j]);
    }
    free(graph[i]);
  }
  free(graph);
}
/*
void free_graph(void** graph, size_t n_hits){
  for(int i = 0; i < 6; i++){
    if(i == 0){
      // X has 3 columns
      void** graphArray = graph[i];
      for(int j = 0; j < 3; j++){
        free(graphArray[j]);
      }
    }
    if(i == 2 || i == 3){
      // Ro and Ri have n_hits columns
      void** graphArray = graph[i];
      for(size_t j = 0; j < n_hits; j++){
        free(graphArray[j]);
      }
    }
    free(graph[i]);
  }
}
*/


void free_sparse(void** sparse){
  for(int i = 0; i < 5; i++){
    if(i == 0){
      // X has 3 columns
      void** sparseArray = *((void***)sparse);
      for(int j = 0; j < 3; j++){
        free(sparseArray[j]);
      }
    }
    if(i == 1){
      // Ra has 4 columns
      void** sparseArray = *((void***)sparse+1);
      for(int j = 0; j < 4; j++){
        free(sparseArray[j]);
      }
    }
    if(i == 2){
      // edge_index has 2 columns
      void** sparseArray = *((void***)sparse+2);
      for(int j = 0; j < 2; j++){
        free(sparseArray[j]);
      }
    }
    free(((void**)(sparse))[i]);
  }
  free(sparse);
}
void** graph_to_sparse(void** graph, size_t n_hits, size_t n_edges){
  vector<size_t>* Ri_rows = new vector<size_t>;
  vector<size_t>* Ri_cols = new vector<size_t>;
  vector<size_t>* Ro_rows = new vector<size_t>;
  vector<size_t>* Ro_cols = new vector<size_t>;
  int** Ri = (int**) graph[2];
  int** Ro = (int**) graph[3];
  size_t i;
  size_t j;
  for(i = 0; i < n_hits; i++){
    for(j = 0; j < n_edges; j++){
      if((Ri[i][j]) != 0){
        Ri_rows->push_back(i);
        Ri_cols->push_back(j);
      }
      if((Ro[i][j]) != 0){
        Ro_rows->push_back(i);
        Ro_cols->push_back(j);
      }
    }
  }
  void** sparse = (void**) malloc(8 * sizeof(void*));
  ((double***) sparse)[0] = (double**) graph[0];
  ((int**) sparse)[1] = (int*) graph[1];
  ((vector<size_t>**)sparse)[2] = Ri_rows;
  ((vector<size_t>**)sparse)[3] = Ri_cols;
  ((vector<size_t>**)sparse)[4] = Ro_rows;
  ((vector<size_t>**)sparse)[5] = Ro_cols;
  ((bool**) sparse)[6] = (bool*) graph[4];
  ((int**) sparse)[7] = (int*) graph[5];
  free_graph_partial(graph, n_hits);
  return sparse;
}

void free_sparse_partial(void** sparse){

  for(int i = 2; i <= 5; i++){
    // Ri_row, Ri_col, Ro_row, Ro_col aren't in graph representation
    vector<size_t>* psparseVector = ((vector<size_t>**)sparse)[i];
    *psparseVector = vector<size_t>();
    free(psparseVector);
  }
  free(sparse);
}

void** sparse_to_graph(void** sparse, size_t n_hits, size_t n_edges){
  int** Ri = (int**) malloc(n_hits * sizeof(void*));
  int** Ro = (int**) malloc(n_hits * sizeof(void*));
  for(size_t i = 0; i < n_hits; i++){
    Ri[i] = (int*) calloc(n_edges, sizeof(int));
    Ro[i] = (int*) calloc(n_edges, sizeof(int));
  }
  vector<size_t> Ri_rows = *((vector<size_t>*)sparse + 0);
  vector<size_t> Ri_cols = *((vector<size_t>*)sparse + 1);
  vector<size_t> Ro_rows = *((vector<size_t>*)sparse + 2);
  vector<size_t> Ro_cols = *((vector<size_t>*)sparse + 3);
  for(size_t i = 0; i < Ri_rows.size(); i++){
    Ri[Ri_cols[i]][Ri_rows[i]] = 1;
  }
  for(size_t i = 0; i < Ro_rows.size(); i++){
    Ro[Ro_cols[i]][Ro_rows[i]] = 1;
  }
  void** graph = (void**) malloc(6 * sizeof(void*));
  void** p_to_arrays = (void**) sparse;
  graph[0] = p_to_arrays[0];
  graph[1] = p_to_arrays[1];
  graph[2] = Ri;
  graph[3] = Ro;
  graph[4] = p_to_arrays[6];
  graph[5] = p_to_arrays[7];
  free_sparse_partial(sparse);
  return graph;
}


void save_graphs(void*** graphs, char** filenames, size_t* n_hits_array, size_t* n_edges_array, size_t n_files){
  for(int i = 0; i < n_files; i++){
    save_graph(graphs[i], filenames[i], n_hits_array[i], n_edges_array[i]);
  }
}
void save_graph(void** sparse, char* filename, size_t n_hits, size_t n_edges){
    // Make a CSV file with columns equiv to sparse version of graph tuple
    // Create an output filestream object


    // NEW: I think sparse no longer needed. Normally the void** argument is called 'graph', then changed to sparse.
    // In interest of not changing too much, I just changed the argument name

    //void** sparse = graph_to_sparse(graph, n_hits, n_edges);


    //filename = "../graphs/TestData.csv"
    std::ofstream file(filename);

    //size_t ri_size = (((vector<size_t>**)sparse)[2])->size();
    //size_t ro_size = (((vector<size_t>**)sparse)[4])->size();

    // Send column names to the stream
    file << "X_0(r),X_1(phi),X_2(z),Ra_0(dr),Ra_1(dphi),Ra_2(dz),Ra_3(dR),edge_index_1,edge_index_2,y,pid\n";

    // Send data to the stream
    size_t maxsize = n_edges;
    for(size_t i = 0; i < maxsize; i++)
    {
        // 3 cols in X are length n_hits
        if(i < n_hits){
          file << ((double**)*(sparse+0))[0][i] << ',' << ((double**)*(sparse+0))[1][i] << ',' << ((double**)*(sparse+0))[2][i] << ',';
        }else{
          file << ',' << ',' << ',';
        }
        // 4 cols in Ra are length n_edges.
        // NEW: 2 cols in edge_index also length n_edges
        if(i < n_edges){
          file << ((double**)*(sparse+1))[0][i] << ',' << ((double**)*(sparse+1))[1][i] << ',' << ((double**)*(sparse+1))[2][i] << ',' << ((double**)*(sparse+1))[3][i]
          << ',' << ((int**)*(sparse+2))[0][i] << ',' << ((int**)*(sparse+2))[1][i] << ',';
        }else{
          file << ',' << ',' << ',' << ',' << ',' << ',';
        }



/*
        // Ri_rows, Ri_cols are length ri_size
        if(i < ri_size){
          file << (((vector<size_t>**)sparse)[2])->at(i) << ',' << (((vector<size_t>**)sparse)[3])->at(i) << ',';
        }else{
          file << ',' << ',';
        }
        // Ro_rows, Ro_cols are length ro_size
        if(i < ro_size){
          file << (((vector<size_t>**)sparse)[4])->at(i) << ',' << (((vector<size_t>**)sparse)[5])->at(i) << ',';
        }else{
          file << ',' << ',';
        }
*/

        // y col is length n_edges
        if(i < n_edges){
          file << (*((bool**)sparse+3))[i] << ',';
        }else{
          file << ',';
        }
        // pid col is length n_hits
        if(i < n_hits){
          file << (*((int**)sparse+4))[i];
        }
        file << '\n';
    }
    // Close the file
    file.close();
    free_sparse(sparse);
}
/* void** load_graph(graph, filename){


  graph = sparse_to_graph(sparse);
  return graph;
}
*/
