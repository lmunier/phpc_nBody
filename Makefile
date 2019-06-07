CXX=g++
NVCC=nvcc
CXXFLAGS= -std=c++11 -g -Wall -O3
CUDAFLAGS= -pg -g -std=c++11 -arch=sm_35 -O3
LDFLAGS=
LIBS= -lcudadevrt

TARGET=nbody

DOC = doc
SRCDIR=src
OBJDIR=build
BINDIR=.

SOURCES  := $(wildcard $(SRCDIR)/*.cu)
INCLUDES := $(wildcard $(SRCDIR)/*.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cu=$(OBJDIR)/%.o)

RM = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(NVCC) $(OBJECTS) $(CUDAFLAGS) -o $@ $(LIBS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cu
	$(NVCC) $(CUDAFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
	@echo "Cleaned up objects."

.PHONY: cleanall
cleanall: clean
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Cleaned up binaries."

.PHONY: doc
doc:
	cd $(DOC) && doxygen Doxyfile
