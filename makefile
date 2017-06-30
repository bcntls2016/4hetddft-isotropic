#
# Makefile para montar el programa DFT4He3d 
#  (Density Functional Theory 4He 3dimensional program)
#

COMP = ifort
CFLAGS = -c -O3 -xAVX -align array64byte -qopenmp -parallel -mkl=parallel -unroll0 -module ./modules
LD_FLAGS = -threads -parallel -qopt-matmul -I${MKLROOT}/include/fftw -mkl=parallel
PROGNAME=4hetddft-isotropic

objs=init_deriv_parallel.o	modules.o	V_impur.o	DFT4He3d.o	derden.o	dimen.o		energy.o\
		fforma.o	fft.o		initcg.o	morse.o\
		mates.o		poten.o		printoutc.o	r_cm.o	redef.o\
		readenc.o	respar.o	term_alfa.o	timer.o		titols.o\
		tstgrid.o	s13adf.o	newder.o\
		steprk.o	steppc.o	potenimp.o	potenimpini.o

.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
	$(COMP)	-o $(PROGNAME) $(objs) $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;

clean:
	rm -f *.o *.bak *.lst modules/*.mod;
distclean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
