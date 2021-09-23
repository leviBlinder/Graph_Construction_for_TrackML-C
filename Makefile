INCLUDE_DIRECTORIES = data_handling data_handling/structs graph_construction
VPATH = $(INCLUDE_DIRECTORIES)
# Compiler:
CC = g++

# Compiler Flags:
# -g (debugging info)
# -Wall (compiler warnings)
CFLAGS = -g -Wall

# Compiler Path Flags:
CPFLAGS = $(addprefix -I , $(INCLUDE_DIRECTORIES))

# The build targets
TARGETBUILD = graph_construction/build_geometric
TARGETMEASURE = graph_construction/measurements/measure_geometric

# Local .cc files to be linked
DEPEND = data_handling/hep_data.cc data_handling/structs/segment.cc data_handling/structs/hit_list.cc data_handling/structs/args.cc graph_construction/geometric_functions.cc


all: $(TARGETBUILD) $(TARGETMEASURE)

$(TARGETMEASURE): $(TARGETMEASURE).cc $(DEPEND)
	$(CC) $(CFLAGS) $(CPFLAGS) $(DEPEND) $(TARGETMEASURE).cc -o $(TARGETMEASURE)

$(TARGETBUILD): $(TARGETBUILD).cc $(DEPEND)
	$(CC) $(CFLAGS) $(CPFLAGS) $(DEPEND) $(TARGETBUILD).cc -o $(TARGETBUILD)


clean:
	rm $(TARGETBUILD) $(TARGETMEASURE)
