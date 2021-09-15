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
TARGET = graph_construction/build_geometric

# Local .cc files to be linked
DEPEND = data_handling/hep_data.cc data_handling/structs/segment.cc data_handling/structs/hit_list.cc


all: $(TARGET)

$
$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) $(CPFLAGS) $(DEPEND) $(TARGET).cc -o $(TARGET)

clean:
	rm $(TARGET)
