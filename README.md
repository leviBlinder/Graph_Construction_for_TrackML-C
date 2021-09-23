# Graph construction for charged particle tracking [C++]

## Tutorial
- (All steps assume starting in home directory)
- In case executables aren't present or up to date, run
  ```
  make
  ```
- **constructing graphs for plots**

  ***1)*** Navigate to and run graph building executable:
  ```
  cd graph_construction
  ./build_geometric
  ```
  [Output graph in .csv form will now be in graphs folder]

  ***2)*** Navigate to and run python script to convert graph data to .npz format
  ```
  cd ../data_handling
  python3 conversionToNPZ.py
  ```
  [Output graph in .npz form will now be in graphs/graphsNPZ folder]

- **making plots for graphs**

  ***1)*** Open jupyter notebooks
  ```
  jupyter notebook
  ```
  ***2)*** Open plotting folder then paper_plots.ipynb in jupyter notebooks.

- **taking measurements for graph construction**

  ***1)*** Navigate to graph construction measurement executable
  ```
  cd graph_construction/measurements
  ```
  ***2)*** Remove previous measurements (this step will likely be automatic in the future)
  ```
  rm statsCompatible.csv
  ```
  ***3)*** Run graph construction measurement executable
  ```
  ./measure_geometric
  ```
  [stats from measurement script will now be accessible in /graph_construction/measurements/statsCompatible.csv]

## Repo Organization
- **data_handling/**
  - structs/: Folder containing definitions of "segment" and "hitlist" objects used to simplify/clarify data structures used in graph construction.

  - conversionToNPZ.py: Transforms all graphs from .csv to .npz. Currently needs to be manually run with:
    ```
    python3 conversionToNPZ.py
    ```

  - hep_data.cc: implements functions to load data from kaggle event files and well as to and from .csv graphs.
- **graph_construction/** [Contains implementations and measurements of graph construction algorithms]
  - measurements/: Folder which will contain averaged measurements of various graph features (purity, efficiency, etc.) of graphs constructed for various pt cuts on 100 different training events.

  - build_geometric.cc: HEP.TrkX construction algorithm modified to build graphs in the pixel barrel+endcaps layers.
    - NOTE: Currently includes measurement-specific functions and process in main [these will be separated moved into measurement folder].
    - NOTE: Currently only works on "evtid = 1000", the first event.
    - NOTE: Configs not yet implemented.
  - build_geometric: Executable file for above algorithm.

- **graphs/**
  [Stores graphs constructed by graph_construction algorithms]
  - graphsNPZ/: Stores graphs transformed from .csv to .npz files for input into plotting.

    - Note .npz file type is not strictly necessary, currently used for compatibility paper_plots.ipynb [modified from [plots](https://github.com/GageDeZoort/interaction_network_paper/blob/pytorch_geometric/plotting/paper_plots.ipynb) in python-based Interaction Network (trackML project)]
- **plotting/**
  - paper_plots.ipynb: Jupyter notebook that loads and plots saved graphs, and loads and plots measurements.
- **train_100_events/**
  [cells, hits, particles, and truth files for 100 events from [Kaggle trackml](https://www.kaggle.com/c/trackml-particle-identification) training dataset. Each event consists of many sensor 'hits' recorded over a very small time frame]
