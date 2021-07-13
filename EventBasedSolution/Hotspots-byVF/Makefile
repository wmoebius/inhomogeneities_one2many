# event driven model  makefile


# The source file
SOURCE = event_driven_model_hotspots.c

#The executable
EXECUTABLES = event_model

#Objects
OBJS = event_driven_model_hotspots.o


#######################################################
############### FLAGS ##################
MPI		 = YES
DEBUG		 = YES
PERFECT          = NO
#######################################

CC = gcc

ifeq ($(MPI), YES)
	CC	= mpicc
	EDEFS	+= -DMPI
endif

ifeq ($(DEBUG), YES)
        CFLAGS = -O0 -g  
	EDEFS	+= -DDEBUG
else
        CFLAGS = -O3 
endif

ifeq ($(PERFECT), YES)
       EDEFS    +=  -DPERFECT
endif

LIBS   = -lm
CFLAGS += $(EDEFS)


# Makefile rules
all: $(EXECUTABLES)

%.o: %.c $(EDEFS)
	$(CC) -c $(CFLAGS) $< -o $@

event_model: $(SOURCE) $(OBJS) 
	$(CC) $(CFLAGS) -o event_model $(OBJS) $(LIBS)

# Tell make that these are phone targets
.PHONY: all help clean

help: doc/Doxyfile $(SOURCE)
	doxygen doc/Doxyfile
	cd doc/latex && make 

clean:
	@rm -f *.o *.c~
	echo "Cleaned."

veryclean:
	@rm -f *.dat *.out *.bin data/*
