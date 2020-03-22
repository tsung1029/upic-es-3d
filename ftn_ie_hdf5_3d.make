#Makefile for 3D parallel PIC codes in new_pbeps32.source

CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile mpiifort compiler at LLNL

#MPIFC = mpiifort
#MPIFC77 = mpiifort

#FC90 = mpiifort
#FC77 = mpiifort
#CC = mpiicc
#FC90 = mpif90
#FC77 = mpif90
#CC = mpicc

#OPTS90 = -O3  -r8 -ipo
#OPTS77 = -O3  -r8 -ipo
#OPTS90 = -r8 -xSSE4.1 -fp-model source -assume minus0
#OPTS77 = -r8 -xSSE4.1 -fp-model source -assume minus0

#CCOPTS = -O3 -ipo
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = nullLOG.o

#HDF_DIR = /usr/local/tools/hdf4
#H5_DIR = /usr/local/tools/hdf5-intel-serial
#SZ_DIR = /usr/local/tools/szip

#INCPATH = -I$(HDF_DIR)/include -I$(H5_DIR)/include -I$(H5_DIR)/lib

#LIBS =  -L$(HDF_DIR)/lib -lmfhdf -ldf -ljpeg -lz \
#        -L$(H5_DIR)/lib -lsz -lhdf5_fortran -lhdf5 \
#        -L$(SZ_DIR)/lib

#LIBS =  -L$(HDF_DIR)/lib -L$(H5_DIR)/lib -L$(SZ_DIR)/lib \
#        -Bstatic -lmfhdf -ldf -lhdf5_fortran -lhdf5 \
#        -Bdynamic -ljpeg -lz -lsz

#################################
#Debug version
MPIFC = ftn
MPIFC77 = ftn
#MPIFC = mpiifort
#MPIFC77 = mpiifort
#MPIFC = mpif90
#MPIFC77 = mpif77

CC = cc
#CC = mpiicc
#CC = mpicc

FC90 = $(MPIFC)
FC77 = $(MPIFC77)

#OPTS90 = -O0 -g -s real64
OPTS90 = -O3 -r8
#OPTS90 = -O3 -s real64
#OPTS90 = -O3 -g -r8 -ipo -no-prec-div -axCORE-AVX-I,AVX,SSE4.2 -xSSE2 -msse2
#OPTS90 = -O3 -r8 -ipo -no-prec-div
OPTS77 = $(OPTS90)
#OPTS90 = -r8 -xSSE4.1 -fp-model source -assume minus0
#OPTS77 = -r8 -xSSE4.1 -fp-model source -assume minus0

#CCOPTS = -O0 -g
CCOPTS = -O3
#CCOPTS = -O3 -ipo -no-prec-div
#CCOPTS = -O3 -g -ipo -no-prec-div -axCORE-AVX-I,AVX,SSE4.2 -xSSE2 -msse2
MOPTS = 
MBOPTS = 
LOPTS =
LEGACY =

# Link convention
UNDERSCORE = FORTRANSINGLEUNDERSCORE

MPIOBJS = nullLOG.o

HDF_VERSION = 1.8.11
#HDF_DIR = /global/u2/b/bwinjum/hdf5-1.8.10
HDF_DIR = $(HDF5_DIR)
#HDF_DIR = /opt/cray/hdf5-parallel/$(HDF_VERSION)/cray/81
#HDF_DIR = /usr/local/tools/hdf5-intel-parallel-$(HDF_VERSION)
HDF_LIBPATH = -L$(HDF_DIR)/lib -lhdf5_fortran -lhdf5  -lhdf5_hl -lhdf5hl_fortran -lz
HDF_INCPATH = -I$(HDF_DIR)/include

#SZIP_LIBPATH = -L$(SZIP_DIR)/lib -lsz
#SZIP_LIBPATH = -L/usr/local/tools/szip/lib -lsz


INCPATH = $(HDF_INCPATH)
LIBPATH = $(HDF_LIBPATH)

################################

# Makefile Absoft compiler with MacOS X

#MPIFC = f90
#MPIFC77 = f77
#FC90 = f90
#FC77 = f77
#CC = gcc

#OPTS90 = -O3 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_
#OPTS77 = -O -N113 -f -N15
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s -N11
#LOPTS = -plainappl
#LEGACY =

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS = $(CARBON)

# Makefile Nag compiler with MacOS X

#MPIFC = f95
#MPIFC77 = f95
#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS = -framework carbon 
#LEGACY = -dusty

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS =

# Makefile IBM compiler with MacOS X

#MPIFC = xlf90
#MPIFC77 = xlf
#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS = $(CARBON)

# Makefile IBM compiler with LAM MPI and MacOS X

#MPIFC = mpif77
#MPIFC77 = mpif77
#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPxlf.o LnxMP.o

#LIBS = $(CARBON)

# Makefile g95 compiler with MacOS X

#MPIFC = g95
#MPIFC77 = g95
#FC90 = g95
#FC77 = g95
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -fstatic
#MBOPTS = -fstatic
#LOPTS =
#LEGACY =

#MPIOBJS = MacMPIf77.o MacMPI_S.o
#MPOBJS = MacMPf77.o MacMP.o

#LIBS = $(CARBON) -lSystemStubs

# Makefile Intel compiler with Linux

#MPIFC = ifc
#MPIFC77 = ifc
#FC90 = ifc
#FC77 = ifc
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPIOBJS = MacMPIf77.o LnxMPI_S.o
#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

# Makefile Intel compiler with MPI and Linux

#MPIFC = mpiifort
#MPIFC77 = mpiifort
#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

# Makefile Pathscale compiler with AMD Opteron and Linux

#MPIFC = mpipathf90
#MPIFC77 = mpipathf90
#FC90 = pathf90
#FC77 = pathf90
#CC = gcc

#OPTS90 = -O3 -r8 
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS = -lpthread
#LEGACY =

#MPIOBJS = nullLOG.o
#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

# Makefile PGI compiler with AMD Opteron and Linux

#MPIFC = mpipgf90
#MPIFC77 = mpipgf90
#FC90 = pgf90
#FC77 = pgf90
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -Msave
#MBOPTS = -Msave
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

#

#GESOBJS = pespush32mod.o psimul32mod.o

SYSOBJS = os-sys-multi-c.o os-sys-multi.o

GESOBJSIE = psimul32mod_ie.o pespush32mod_ie.o

MGESOBJS = mpespush32mod.o mpsimul32mod.o

ESOBJS = globals.o pinit32mod.o \
pfft32mod.o pfield32mod.o pdiag32mod.o p0mod.o mp0mod.o \
pinit32lib.o pfft32lib.o pfield32lib.o \
pdiag32lib.o p0lib.o 

ESOBJSJF = pinit32mod_jf.o \
ext_driver32_jf.o p32mod_jf.o pinit32lib_jf.o
#ppush2mod_jf.o ppush2lib_jf.o pfield2mod_jf.o

ESOBJSIE = pinit32mod_ie.o sqnmod_ie.o p32mod_ie.o p32lib_ie.o prbpush32mod_ie.o \
prbpush32lib_ie.o ppush32mod_ie.o ppush32lib_ie.o \
string.o os_util_hdf5.o data_utils.o p_os_util_hdf5.o pdata_utils.o \
hdf5_write32_ie.o diag32_ie.o par_track_ie.o

MESOBJS = mprbpush32mod.o mppush32mod.o mpfft32mod.o mp32mod.o \
mprbpush32lib.o mppush32lib.o mpfft32lib.o mp32lib.o 

GEMOBJS = pempush32mod.o

MGEMOBJS = mpempush32mod.o

EMOBJS = pbpush32mod.o pbpush32lib.o

MEMOBJS = mpbpush32mod.o mpbpush32lib.o

DESOBJS = pdfield32mod.o pbfield32mod.o pcfield32mod.o pnpfield32mod.o \
pdfield32lib.o pbfield32lib.o pcfield32lib.o

NPOBJS = nullMP.o

# Linkage rules

all : new_pbeps32_ie.out

#all : new_pbeps32.out new_d0_pbeps32.out new_pbbeps32.out new_d0_pbbeps32.out \
#     threaded

threaded: new_mpbeps32.out new_d0_mpbeps32.out new_mpbbeps32.out new_d0_mpbbeps32.out

new_pbeps32_ie.out : new_pbeps32_ie.o $(SYSOBJS) $(GESOBJSIE) $(ESOBJS) $(ESOBJSJF) $(ESOBJSIE) $(MPIOBJS) $(NPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_pbeps32_ie.out \
	new_pbeps32_ie.o $(SYSOBJS) $(GESOBJS) $(GESOBJSIE) $(ESOBJS) $(ESOBJSJF) $(ESOBJSIE) $(MPIOBJS) $(NPOBJS) $(LIBPATH) $(INCPATH)

new_pbeps32.out : new_pbeps32.o $(GESOBJS) $(ESOBJS) $(MPIOBJS) $(NPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_pbeps32.out \
	new_pbeps32.o $(GESOBJS) $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(LIBS)

new_d0_pbeps32.out : new_d0_pbeps32.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(MPIOBJS) \
                     $(NPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbeps32.out \
	new_d0_pbeps32.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) $(LIBS)

new_mpbeps32.out : new_mpbeps32.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPIOBJS) \
                   $(MPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbeps32.out \
    new_mpbeps32.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) $(LIBS)

new_d0_mpbeps32.out : new_d0_mpbeps32.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) \
                      $(MPIOBJS) $(MPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbeps32.out \
	new_d0_mpbeps32.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) $(MPIOBJS) \
	$(MPOBJS) $(LIBS)

new_pbbeps32.out : new_pbbeps32.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS)
	$(MPIFC) $(OPTS90) $(OPTS90) $(LOPTS) -o new_pbbeps32.out \
	new_pbbeps32.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) $(LIBS)

new_d0_pbbeps32.out : new_d0_pbbeps32.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) \
                      $(MPIOBJS) $(NPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbbeps32.out \
	new_d0_pbbeps32.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) \
    $(LIBS)

new_mpbbeps32.out : new_mpbbeps32.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                    $(MEMOBJS) $(MPIOBJS) $(MPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbbeps32.out \
    new_mpbbeps32.o $(ESOBJS) $(MESOBJS) $(EMOBJS) $(MGEMOBJS) $(MEMOBJS) $(MPIOBJS) \
    $(MPOBJS) $(LIBS)

new_d0_mpbbeps32.out : new_d0_mpbbeps32.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                       $(MEMOBJS) $(DESOBJS) $(MPIOBJS) $(MPOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbbeps32.out \
    new_d0_mpbbeps32.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
	$(DESOBJS) $(MPIOBJS) $(MPOBJS) $(LIBS)

# Compilation rules

#MacMPIcf.o : MacMPIcf.c
#	$(CC) $(CCOPTS) -c MacMPIcf.c

MacMPIf77.o : MacMPIf77.c
	$(CC) $(CCOPTS) -c MacMPIf77.c

MacMPIxlf.o : MacMPIxlf.c
	$(CC) $(CCOPTS) -c MacMPIxlf.c

MacMPI_S.o : MacMPI_S.c
	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMPI_S.c

LnxMPI_S.o : LnxMPI_S.c
	$(CC) $(CCOPTS) -c LnxMPI_S.c

#MacMPcf.o : MacMPcf.c
#	$(CC) $(CCOPTS) -c MacMPcf.c

MacMPf77.o : MacMPf77.c
	$(CC) $(CCOPTS) -c MacMPf77.c

#MacMPxlf.o : MacMPxlf.c
#	$(CC) $(CCOPTS) -c MacMPxlf.c

MacMP.o : MacMP.c
	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMP.c

LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c

os-sys-multi-c.o : os-sys-multi-c.c
	$(CC) $(CCOPTS) -D$(UNDERSCORE) -c os-sys-multi-c.c

os-sys-multi.o : os-sys-multi.f90 os-sys-multi-c.c
	$(FC90) $(OPTS90) -c os-sys-multi.f90

nullLOG.o : nullLOG.f
	$(FC77) $(OPTS77) -c nullLOG.f

nullMP.o : nullMP.f
	$(FC77) $(OPTS77) -c nullMP.f

p0lib.o : p0lib.f
	$(MPIFC77) $(OPTS77) $(LEGACY) -c p0lib.f

mp32lib.o : mp32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mp32lib.f

p32lib_ie.o : p32lib_ie.f
	$(FC77) $(OPTS77) $(LEGACY) -c p32lib_ie.f

p32lib.o : p32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c p32lib.f

pinit32lib.o : pinit32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pinit32lib.f

pinit32lib_jf.o : pinit32lib_jf.f
	$(FC77) $(OPTS77) $(LEGACY) -c pinit32lib_jf.f

mppush32lib.o : mppush32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mppush32lib.f

ppush32lib.o : ppush32lib.f
	$(FC77) $(OPTS77) -c ppush32lib.f

ppush32lib_ie.o : ppush32lib_ie.f
	$(FC77) $(OPTS77) -c ppush32lib_ie.f

mpbpush32lib.o : mpbpush32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mpbpush32lib.f

pbpush32lib.o : pbpush32lib.f
	$(FC77) $(OPTS77) -c pbpush32lib.f

mprbpush32lib.o : mprbpush32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mprbpush32lib.f

prbpush32lib.o : prbpush32lib.f
	$(FC77) $(OPTS77) -c prbpush32lib.f

prbpush32lib_ie.o : prbpush32lib_ie.f
	$(FC77) $(OPTS77) -c prbpush32lib_ie.f

mpfft32lib.o : mpfft32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mpfft32lib.f

pfft32lib.o : pfft32lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pfft32lib.f

pfield32lib.o : pfield32lib.f
	$(FC77) $(OPTS77) -c pfield32lib.f

pdfield32lib.o : pdfield32lib.f
	$(FC77) $(OPTS77) -c pdfield32lib.f

pbfield32lib.o : pbfield32lib.f
	$(FC77) $(OPTS77) -c pbfield32lib.f

pcfield32lib.o : pcfield32lib.f
	$(FC77) $(OPTS77) -c pcfield32lib.f

pdiag32lib.o : pdiag32lib.f
	$(FC77) $(OPTS77) -c pdiag32lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

pinit32mod.o : pinit32mod.f globals.o
	$(FC90) $(OPTS90) -c pinit32mod.f

mppush32mod.o : mppush32mod.f ppush32mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mppush32mod.f

ppush32mod.o : ppush32mod.f p0mod.o
	$(FC90) $(OPTS90) -c ppush32mod.f

ppush32mod_ie.o : ppush32mod_ie.f p0mod.o
	$(FC90) $(OPTS90) -c -free ppush32mod_ie.f

mpbpush32mod.o : mpbpush32mod.f pbpush32mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpbpush32mod.f

pbpush32mod.o : pbpush32mod.f p0mod.o
	$(FC90) $(OPTS90) -c pbpush32mod.f

mprbpush32mod.o : mprbpush32mod.f prbpush32mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mprbpush32mod.f

prbpush32mod.o : prbpush32mod.f p0mod.o
	$(FC90) $(OPTS90) -c prbpush32mod.f

prbpush32mod_ie.o : prbpush32mod_ie.f p0mod.o
	$(FC90) $(OPTS90) -c -free prbpush32mod_ie.f

mpfft32mod.o : mpfft32mod.f pfft32mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpfft32mod.f

pfft32mod.o : pfft32mod.f p0mod.o
	$(FC90) $(OPTS90) -c pfft32mod.f

pfield32mod.o : pfield32mod.f globals.o
	$(FC90) $(OPTS90) -c pfield32mod.f

pdfield32mod.o : pdfield32mod.f globals.o
	$(FC90) $(OPTS90) -c pdfield32mod.f

pbfield32mod.o : pbfield32mod.f globals.o
	$(FC90) $(OPTS90) -c pbfield32mod.f

pcfield32mod.o : pcfield32mod.f globals.o
	$(FC90) $(OPTS90) -c pcfield32mod.f

pnpfield32mod.o : pnpfield32mod.f pfield32mod.o pdfield32mod.o pbfield32mod.o \
                 pcfield32mod.o
	$(FC90) $(OPTS90) -c pnpfield32mod.f

pdiag32mod.o : pdiag32mod.f pinit32mod.o
	$(FC90) $(OPTS90) -c pdiag32mod.f
	
par_track_ie.o : par_track_ie.f pinit32mod_jf.o
	$(FC90) $(OPTS90) $(HDF_INCPATH) -c -free par_track_ie.f

diag32_ie.o : diag32_ie.f pinit32mod.o pinit32mod_ie.o pinit32mod_jf.o hdf5_write32_ie.o p32mod_jf.o
	$(FC90) $(OPTS90) -c -free diag32_ie.f
	
sqnmod_ie.o : sqnmod_ie.f
	$(FC90) $(OPTS90) -c -free sqnmod_ie.f

hdf5_write32_ie.o : hdf5_write32_ie.f string.o os_util_hdf5.o data_utils.o p_os_util_hdf5.o pdata_utils.o pinit32mod_jf.o
	$(FC90) $(OPTS90) -c -free hdf5_write32_ie.f

string.o : string.f
	$(FC90)	$(OPTS90) -c -free string.f

os_util_hdf5.o : os_util_hdf5.f
	$(FC90) $(OPTS90) $(HDF_INCPATH) -c -free os_util_hdf5.f

data_utils.o : data_utils.f os_util_hdf5.o string.o
	$(FC90) $(OPTS90) -c -free data_utils.f

p_os_util_hdf5.o : p_os_util_hdf5.f os_util_hdf5.o
	$(FC90) $(OPTS90) $(HDF_INCPATH) -c -free p_os_util_hdf5.f

pdata_utils.o : pdata_utils.f p_os_util_hdf5.o string.o
	$(FC90) $(OPTS90) -c -free pdata_utils.f

#pdb_write32_ie.o : pdb_write32_ie.f
#	$(FC90) $(OPTS90) -c -free pdb_write32_ie.f

pinit32mod_jf.o : pinit32mod_jf.f
	$(FC90) $(OPTS90) -c -free pinit32mod_jf.f

pinit32mod_ie.o : pinit32mod_ie.f
	$(FC90) $(OPTS90) -c -free pinit32mod_ie.f
	
ext_driver32_jf.o : ext_driver32_jf.f
	$(FC90) $(OPTS90) -c -free ext_driver32_jf.f

p0mod.o : p0mod.f globals.o
	$(FC90) $(OPTS90) -c p0mod.f

mp32mod.o : mp32mod.f p32mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mp32mod.f

p32mod.o : p32mod.f p0mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c p32mod.f

p32mod_ie.o : p32mod_ie.f p0mod.o par_track_ie.o
	$(FC90) $(OPTS90) $(LEGACY) -c -free p32mod_ie.f

p32mod_jf.o : p32mod_jf.f p0mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c -free p32mod_jf.f

mpespush32mod.o : mpespush32mod.f mprbpush32mod.o mppush32mod.o mpfft32mod.o \
                  mp32mod.o
	$(FC90) $(OPTS90) -c mpespush32mod.f

pespush32mod.o : pespush32mod.f prbpush32mod.o ppush32mod.o pfft32mod.o p32mod.o
	$(FC90) $(OPTS90) -c pespush32mod.f

pespush32mod_ie.o : pespush32mod_ie.f prbpush32mod_ie.o ppush32mod_ie.o pfft32mod.o p32mod_ie.o
	$(FC90) $(OPTS90) -c pespush32mod_ie.f

mpempush32mod.o : mpempush32mod.f mprbpush32mod.o mpbpush32mod.o mppush32mod.o \
                  mpfft32mod.o mp32mod.o
	$(FC90) $(OPTS90) -c mpempush32mod.f

pempush32mod.o : pempush32mod.f prbpush32mod.o pbpush32mod.o ppush32mod.o \
                 pfft32mod.o p32mod.o
	$(FC90) $(OPTS90) -c pempush32mod.f

psimul32mod.o : psimul32mod.f pdiag32mod.o pespush32mod.o pfield32mod.o
	$(FC90) $(OPTS90) -c psimul32mod.f

psimul32mod_ie.o : psimul32mod_ie.f pdiag32mod.o pespush32mod_ie.o pfield32mod.o hdf5_write32_ie.o os-sys-multi.o
	$(FC90) $(OPTS90) $(HDF_INCPATH) -free -c psimul32mod_ie.f

mpsimul32mod.o : psimul32mod.f pdiag32mod.o mpespush32mod.o pfield32mod.o
	$(FC90) $(OPTS90) -o mpsimul32mod.o -c psimul32mod.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f

new_pbeps32.o : new_pbeps32.f psimul32mod.o pfield32mod.o pdiag32mod.o \
                mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbeps32.f

new_pbeps32_ie.o : new_pbeps32_ie.f psimul32mod_ie.o pfield32mod.o pdiag32mod.o \
                mp0mod.o $(SYSOBJS) $(ESOBJSJF) $(ESOBJSIE) pinit32mod_ie.o
	$(FC90) $(OPTS90) $(MOPTS) -c -free new_pbeps32_ie.f

new_mpbeps32.o : new_pbeps32.f mpsimul32mod.o pfield32mod.o pdiag32mod.o \
                 mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mpbeps32.o -c new_pbeps32.f

new_pbbeps32.o : new_pbbeps32.f pempush32mod.o pfield32mod.o pdiag32mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbbeps32.f

new_mpbbeps32.o : new_pbbeps32.f mpempush32mod.o pfield32mod.o pdiag32mod.o \
                  pdiag32mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mpbbeps32.o -c new_pbbeps32.f

new_d0_pbeps32.o : new_d0_pbeps32.f pespush32mod.o pnpfield32mod.o pdiag32mod.o \
                   mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_pbeps32.f

new_d0_mpbeps32.o : new_d0_pbeps32.f mpespush32mod.o pnpfield32mod.o pdiag32mod.o \
                    mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_d0_mpbeps32.o -c new_d0_pbeps32.f

new_d0_pbbeps32.o : new_d0_pbbeps32.f pempush32mod.o pnpfield32mod.o pdiag32mod.o \
                    mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_pbbeps32.f

new_d0_mpbbeps32.o : new_d0_pbbeps32.f mpempush32mod.o pnpfield32mod.o pdiag32mod.o \
                     mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -o new_d0_mpbbeps32.o -c new_d0_pbbeps32.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f *.out
