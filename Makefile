# Makefile template from
# http://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/

CC := g++
SRCDIR := src
BUILDDIR := build
TARGET := bin/disorder_cooling

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
WARNING := -Wall -Werror -Wextra -Wfloat-equal -ansi -pedantic
CFLAGS := -pipe -O2 -std=c++17 -march=native -mtune=native -flto -funroll-loops \
	-finline-functions -fno-stack-protector -ftree-vectorize -fopenmp
LIB := -L lib -fopenmp
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

test:
	$(CC) $(CFLAGS) test/test.cpp $(INC) $(LIB) -o bin/test

.PHONY: clean
