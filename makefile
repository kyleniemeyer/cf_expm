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

OBJ = main.o cf.o

# Paths
INCLUDES    = -I.
LIBS  = -lm -llapack -lfftw3

#flags
ifeq ("$(CC)", "gcc")
  ifeq ("$(L)", "0")
    FLAGS = -O0 -g3 -fbounds-check -Wunused-variable -Wunused-parameter \
	        -Wall -pedantic -ftree-vrp -std=c99
  else ifeq ("$(L)", "4")
    FLAGS = -O3 -ffast-math -std=c99
  endif
endif

all: cf_expm

%.o : %.c
	$(CC) -c -o $@ $< $(FLAGS) $(INCLUDES)

cf_expm : $(OBJ)
		$(LINK) -o $@ $(OBJ) $(LIBS) $(FLAGS)

doc : $(OBJ)
		$(DOXY)

.PHONY : clean		
clean :
	rm -f $(OBJ) cf_expm