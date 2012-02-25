LIB=minpack
FFLAGS=-O
OBJ = \
	covar.o	\
	dmchar.o	\
	dogleg.o	\
	dpmpar.o	\
	enorm.o	\
	errjac.o	\
	fdjac1.o	\
	fdjac2.o	\
	grdfcn.o	\
	hesfcn.o	\
	hybipt.o	\
	hybrd.o	\
	hybrd1.o	\
	hybrj.o	\
	hybrj1.o	\
	lhesfcn.o	\
	lmder.o	\
	lmder1.o	\
	lmdif.o	\
	lmdif1.o	\
	lmdipt.o	\
	lmpar.o	\
	lmstr.o	\
	lmstr1.o	\
	objfcn.o	\
	ocpipt.o	\
	qform.o	\
	qrfac.o	\
	qrsolv.o	\
	r1mpyq.o	\
	r1updt.o	\
	rwupdt.o	\
	ssqfcn.o	\
	ssqjac.o	\
	vecfcn.o	\
	vecjac.o

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	ln -s /netlib/netlib/minpack/lib$(LIB).a /usr/local/lib
	rm *.o

test: test.o
	f77 test.o -l$(LIB)
	time a.out
