
# ===== Compiler selection =====
FC ?= gfortran

# ===== Flags per compiler =====
ifeq ($(FC),gfortran)
    FFLAGS = -O3 -march=native -ffast-math -funroll-loops -fopenmp
    LAPACK = -llapack -lblas
endif

ifeq ($(FC),ifort)
    FFLAGS = -O3 -xHost -qopenmp -ipo -fp-model fast=2
    LAPACK = -mkl
endif

ifeq ($(FC),ifx)
    FFLAGS = -O3 -xHost -qopenmp -ipo -fp-model fast=2
    LAPACK = -mkl
endif

SRC = src


# --- Core modules (NO main here) ---
CORE_OBJS = \
       lobatto.o \
       constants.o \
       structure_parameters.o \
       dynamic_parameters.o \
       exploit_parameters.o \
       math_util.o \
       util.o \
       fedvr.o \
       fedvr_topology.o \
       fedvr_conf_struct.o \
       fedvr_derivative_ops.o \
       space_time_ops.o \
       global_assembly.o \
       propagation.o \
       observables.o \
       io_module.o \
       conv_tests.o

# --- Executables ---
STRUCTURE_EXE = structure
DYNAMIC_EXE  = dynamic
EXPLOIT_EXE  = exploit

all: $(STRUCTURE_EXE)

# --- Structure build ---
$(STRUCTURE_EXE): $(CORE_OBJS) main_structure.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Future ---
$(DYNAMIC_EXE): $(CORE_OBJS) main_dynamic.o
	$(FC) $^ -o $@ $(LAPACK)

$(EXPLOIT_EXE): $(CORE_OBJS) main_exploit.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Compilation rule ---
%.o: $(SRC)/%.f
	$(FC) $(FFLAGS) -c $<

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

# --- Dependency include ---
-include *.d

.PHONY: clean
clean:
	rm -f *.o *.d *.mod fort.* $(STRUCTURE_EXE) $(DYNAMIC_EXE) $(EXPLOIT_EXE)

