#include <utility>
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>

#include "args.h"

using namespace std;
typedef struct args* args_T;

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

void args::parse_args(int argc, const char** argv){
  /* Parse command line arguments. */
  assert(argv != NULL);
  //could just assign 1st arg to config, but not sure about order
  //parse_arg_int(argc, argv, "--n-workers", &(args->config));
  parse_arg_int(argc, argv, "--n-workers", &(n_workers));
  parse_arg_int(argc, argv, "--task", &(task));
  parse_arg_int(argc, argv, "--n-tasks", &(n_tasks));
  parse_arg_bool(argc, argv, "-v", &(verbose));
  parse_arg_bool(argc, argv, "--show-config", &(show_config));
  parse_arg_bool(argc, argv, "--interactive", &(interactive));
  parse_arg_int(argc, argv, "--start-evtid", &(start_evtid));
  parse_arg_int(argc, argv, "--end-evtid", &(end_evtid));
  return;
}

args::args(){
  config = NULL;
  n_workers = 1;
  task = 0;
  n_tasks = 1;
  verbose = false;
  show_config = false;
  interactive = false;
  start_evtid = 1000;
  end_evtid = 2700;
}
