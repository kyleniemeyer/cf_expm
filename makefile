SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .c .o

# Compilers
CC    = gcc
LINK   = $(CC)

DOXY = doxygen

# L=0 for testing, L=4 for optimization
ifndef L
  L = 4
endif

# object file directory
ODIR = obj

# source files
_OBJ = main.o cf.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Paths
INCLUDES = -I.
LIBS = -lm -llapack -lfftw3

#flags
ifeq ("$(CC)", "gcc")
ifeq ("$(L)", "0")
 FLAGS = -O0 -g3 -fbounds-check -Wunused-variable -Wunused-parameter \
	       -Wall -pedantic -ftree-vrp -std=c99 \
			   -Wextra -fno-strict-aliasing -fwrapv -v -save-temps
else ifeq ("$(L)", "4")
  FLAGS = -O3 -std=c99
  #FLAGS += -ffast-math 
endif
endif

default: $(ODIR) all

all: cf_expm
	
$(ODIR):
	mkdir $(ODIR)

$(ODIR)/%.o : %.c
	$(CC) -c -o $@ $< $(FLAGS) $(INCLUDES)

cf_expm : $(OBJ)
		$(LINK) -o $@ $(OBJ) $(LIBS) $(FLAGS)

doc : $(OBJ)
		$(DOXY)

.PHONY : clean		
clean :
	rm -f $(OBJ) cf_expm