TARGET = 3D

#SRC  = modules.f90 SPmain.f90 parser.f90 init.f90 allocation.f90 allocateell.f90 3D.f90 cadenas.f90 cadenasMK.f90 fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 monomers.definitions-onck.f90 chains.definitions.f90 sphere.f90 kapfromfile.f90

FF = mpif77 #${F90}
SRC = modules.f90 SPmain.f90 channel.f90 PBC.f90 parser.f90 init.f90 allocation.f90 allocatencha.f90 allocateell.f90 3D.f90  allocatecpp.f90  cadenas.f90 cadenas_b.f90 cadenas_b2.f90  fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 transform.f90 testsystem.f90 testsystemc.f90 testsystemr.f90 monomers.definitions.f90 chains.definitions.f90

HOST=$(shell hostname)
$(info HOST is ${HOST})


# some definitions
SHELL = /bin/bash
FFLAGS= -O3 #-fbacktrace -fbounds-check # -O3

ifeq ($(HOST),piluso.rosario-conicet.gov.ar)
LFLAGS = -L/home/mtagliazucchi.inquimae/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif


ifeq ($(HOST),pear)
LFLAGS=-L/home/mario/software/KINSOL/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux
endif

ifeq ($(HOST),cnode01)
LFLAGS = -L/home/mtagliazucchi/software/Kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif

ifeq ($(HOST),master) 
LFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
endif

ifeq ($(HOST),mate.bme.northwestern.edu) 
LFLAGS = -L/home/mario/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
endif

ifeq ($(HOST),quser13) 
LFLAGS=-L/home/khl4149/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4 -L/hpc/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64/ -L/lib/../lib64 -L/lib/../lib64/ -L/usr/lib/../lib64 -L/usr/lib/../lib64/ -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../ -L/lib64 -L/lib/ -L/usr/lib64 -L/usr/lib -limf -lm -lifport -lifcore -lsvml -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl
endif

ifeq ($(HOST),quser12)
LFLAGS=-L/home/khl4149/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4 -L/hpc/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64/ -L/lib/../lib64 -L/lib/../lib64/ -L/usr/lib/../lib64 -L/usr/lib/../lib64/ -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../ -L/lib64 -L/lib/ -L/usr/lib64 -L/usr/lib -limf -lm -lifport -lifcore -lsvml -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl
endif

ifeq ($(HOST),quser11)
LFLAGS=-L/home/khl4149/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4 -L/hpc/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64/ -L/lib/../lib64 -L/lib/../lib64/ -L/usr/lib/../lib64 -L/usr/lib/../lib64/ -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../ -L/lib64 -L/lib/ -L/usr/lib64 -L/usr/lib -limf -lm -lifport -lifcore -lsvml -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl
endif

ifeq ($(HOST),quser10)
LFLAGS=-L/home/khl4149/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4 -L/hpc/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64/ -L/lib/../lib64 -L/lib/../lib64/ -L/usr/lib/../lib64 -L/usr/lib/../lib64/ -L/opt/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/ipp/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/ -L/opt/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../ -L/lib64 -L/lib/ -L/usr/lib64 -L/usr/lib -limf -lm -lifport -lifcore -lsvml -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl
endif

ifeq ($(HOST),aci-service-1.chtc.wisc.edu)
LFLAGS=-L/home/khuang28/sundials-2.5.0-build/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
FF = /home/khuang28/mpich/bin/mpif90
endif

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

VER = ~/bin/3D_curvature

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(GFLAGS)
	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(GFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif
