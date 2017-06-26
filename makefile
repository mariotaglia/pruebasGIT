# put the name of the target program here
TARGET = 3D

SRC = modules.f90 SPmain.f90 channel.f90 PBC.f90 parser.f90 init.f90 allocation.f90 allocatencha.f90 allocateell.f90 3D.f90  allocatecpp.f90  cadenas.f90 cadenas_b.f90 cadenas_b2.f90  fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 ellipsoid.f90 dielectric.f90 transform.f90 testsystem.f90 testsystemc.f90 monomers.definitions.f90 chains.definitions.f90

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 #-fbounds-check ${F90FLAGS} # -O3 -traceback -check all ${F90FLAGS}

#LDFLAGS= -L/export/apps/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

#LFFLAGS=$(LDFLAGS)

LDFLAGS=-L/home/khuang28/sundials-2.5.0-build/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

LFLAGS=$(LDFLAGS)

#LFLAGS= -L/home/khuang/sundials-2.5.0-build/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

VER = ~/bin/3d_channel_g2

#FF = /state/partition1/apps/openmpi-1.6.3/bin/mpif77
#FF = /share/apps/openmpi-1.10.0_no_ib/bin/mpif90
FF = /home/khuang28/mpich/bin/mpif90
#FF = mpif77 #${F90}
#FF = /opt/openmpi/bin/mpif90


all:    $(TARGET)

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
