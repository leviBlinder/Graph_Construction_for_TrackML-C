#ifndef ARGS_H
#define ARGS_H

struct args {
public:
  char* config;
  int n_workers;
  int task;
  int n_tasks;
  bool verbose;
  bool show_config;
  bool interactive;
  int start_evtid;
  int end_evtid;

  args();
  void parse_args(int argc, const char** argv);
};

typedef struct args* args_T;

#endif
