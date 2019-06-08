CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3
LDFLAGS=
LDLIBS=

TARGET=nbody

DOC = doc
SRCDIR=src
OBJDIR=build
BINDIR=.
CXXFLAGS+= -I./$(SRCDIR)/*.hpp

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

RM = rm -f

CD = cd doc/
CD_ = cd..

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

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

