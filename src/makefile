# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

gccopt= gcc -O2 -funroll-loops
gccdeb= gcc -ggdb -Wall -pedantic
gppdeb = g++ $(CCDEB) 
#CCC = $(gccopt)
CCC = $(gccdeb)

CFLAGS=
LINUXLIBS= -lm -lpthread

LIBS= $(LINUXLIBS) 


#---------------- GUROBI -----------------------------------
kikiGRB = /home/dano/thirdparty/gurobi/gurobi201/linux32
furrybabyGRB = /opt/gurobi602/linux64
rockwaldGRB = /opt/gurobi/gurobi451/linux64
fatpuppyGRB = /opt/gurobi/gurobi563/linux64

GRB = $(furrybabyGRB)
GRBLNFLAGS = -L$(GRB)/lib -lgurobi60

#GRB = $(fatpuppyGRB)

GRBLPINCDIR = $(GRB)/include

#cplex
#------------------------- true Unix ----------------------------
LIBNAME = libcplex.a 
#------------------------------- cygwin -------------------------
#LIBNAME = cplex112.lib

CPX = $(fatbabyCPX)
#CPX = $(houndbabyCPX)
#CPX = $(rockyCPX)
#CPX = $(racepointCPX)

CPXLNFLAGS = $(CPX)/libcplex.a

CPXLPINCDIR = /opt/ibm/ILOG/CPLEX_Enterprise_Server126/CPLEX_Studio/cplex/include

#------------ LPSOLVER CHOICE MADE HERE ---------------

#CCCLNFLAGS = $(CPXLNFLAGS) $(LIBS)
#CCCLNFLAGS = $(GRBLNFLAGS) $(LIBS)
#LPINCDIR = $(GRBLPINCDIR) #none for this project

CCCLNFLAGS = $(LIBS)
LPINCDIR =

CCCFLAGS = 

../bin/%.o: %.c
	@echo compiling $*.c with $(CCC) $(CCCFLAGS)
	@$(CCC) $(CCCFLAGS) -c $< -o $@

RPOWER = ../bin/mainrpower.o ../bin/power.o ../bin/utilities.o

rpower: $(RPOWER)
	$(CCC) $(CCCFLAGS) -o ../bin/rpower $(RPOWER) $(CCCLNFLAGS)

clean:
	rm ../bin/*

seep:
	rm *~
	rm *#
