#!/bin/bash

#
# Download the old minpack and test files from Netlib
# and compile and run all the tests.
#
# This was used to generate the comparisons for the unit tests.
#

COMPILER=gfortran

rm -rf lib
rm -rf test

mkdir lib
mkdir test

# minpack library files:
wget https://netlib.org/minpack/chkder.f
wget https://netlib.org/minpack/dogleg.f
wget https://netlib.org/minpack/dpmpar.f
wget https://netlib.org/minpack/enorm.f
wget https://netlib.org/minpack/fdjac1.f
wget https://netlib.org/minpack/fdjac2.f
wget https://netlib.org/minpack/hybrd.f
wget https://netlib.org/minpack/hybrd1.f
wget https://netlib.org/minpack/hybrj.f
wget https://netlib.org/minpack/hybrj1.f
wget https://netlib.org/minpack/lmder.f
wget https://netlib.org/minpack/lmder1.f
wget https://netlib.org/minpack/lmdif.f
wget https://netlib.org/minpack/lmdif1.f
wget https://netlib.org/minpack/lmpar.f
wget https://netlib.org/minpack/lmstr.f
wget https://netlib.org/minpack/lmstr1.f
wget https://netlib.org/minpack/qform.f
wget https://netlib.org/minpack/qrfac.f
wget https://netlib.org/minpack/qrsolv.f
wget https://netlib.org/minpack/r1mpyq.f
wget https://netlib.org/minpack/r1updt.f
wget https://netlib.org/minpack/rwupdt.f

mv *.f lib

# test programs:
wget https://netlib.org/minpack/ex/file15
wget https://netlib.org/minpack/ex/file16
wget https://netlib.org/minpack/ex/file17
wget https://netlib.org/minpack/ex/file18
wget https://netlib.org/minpack/ex/file19
wget https://netlib.org/minpack/ex/file20

mv file15 ./test/file15.f
mv file16 ./test/file16.f
mv file17 ./test/file17.f
mv file18 ./test/file18.f
mv file19 ./test/file19.f
mv file20 ./test/file20.f

# modify the tests to read the data from the files
# note: this is mac sed:
sed -i -e 's|         READ (NREAD,50) NPROB,N,NTRIES|         if (ic==0) open(unit=NREAD,file="file21",status="OLD");READ (NREAD,50) NPROB,N,NTRIES|g' ./test/file15.f
sed -i -e 's|         READ (NREAD,50) NPROB,N,NTRIES|         if (ic==0) open(unit=NREAD,file="file21",status="OLD");READ (NREAD,50) NPROB,N,NTRIES|g' ./test/file16.f
sed -i -e 's|         READ (NREAD,50) NPROB,N,M,NTRIES|         if (ic==0) open(unit=NREAD,file="file22",status="OLD");READ (NREAD,50) NPROB,N,M,NTRIES|g' ./test/file17.f
sed -i -e 's|         READ (NREAD,50) NPROB,N,M,NTRIES|         if (ic==0) open(unit=NREAD,file="file22",status="OLD");READ (NREAD,50) NPROB,N,M,NTRIES|g' ./test/file18.f
sed -i -e 's|         READ (NREAD,50) NPROB,N,M,NTRIES|         if (ic==0) open(unit=NREAD,file="file22",status="OLD");READ (NREAD,50) NPROB,N,M,NTRIES|g' ./test/file19.f
sed -i -e 's|      LDFJAC = 10|      LDFJAC = 10; open(unit=NREAD,file="file23",status="OLD")|g' ./test/file20.f

sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file15.f
sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file16.f
sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file17.f
sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file18.f
sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file19.f
sed -i -e 's|      DATA NREAD,NWRITE /5,6/|      DATA NREAD,NWRITE /500,6/|g' ./test/file20.f

# input files:
wget https://netlib.org/minpack/ex/file21
wget https://netlib.org/minpack/ex/file22
wget https://netlib.org/minpack/ex/file23

# compile tests:
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file15.f -o test_hybrd_O0_$COMPILER
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file16.f -o test_hybrj_O0_$COMPILER
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file17.f -o test_lmder_O0_$COMPILER
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file18.f -o test_lmstr_O0_$COMPILER
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file19.f -o test_lmdif_O0_$COMPILER
$COMPILER -O0 -ffixed-line-length-none ./lib/*.f ./test/file20.f -o test_chkder_O0_$COMPILER

# run tests:
./test_hybrd_O0_$COMPILER > output_test_hybrd_O0_$COMPILER.txt
./test_hybrj_O0_$COMPILER > output_test_hybrj_O0_$COMPILER.txt
./test_lmder_O0_$COMPILER > output_test_lmder_O0_$COMPILER.txt
./test_lmstr_O0_$COMPILER > output_test_lmstr_O0_$COMPILER.txt
./test_lmdif_O0_$COMPILER > output_test_lmdif_O0_$COMPILER.txt
./test_chkder_O0_$COMPILER > output_test_chkder_O0_$COMPILER.txt

# compile tests:
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file15.f -o test_hybrd_O2_$COMPILER
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file16.f -o test_hybrj_O2_$COMPILER
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file17.f -o test_lmder_O2_$COMPILER
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file18.f -o test_lmstr_O2_$COMPILER
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file19.f -o test_lmdif_O2_$COMPILER
$COMPILER -O2 -ffixed-line-length-none ./lib/*.f ./test/file20.f -o test_chkder_O2_$COMPILER

# run tests:
./test_hybrd_O2_$COMPILER > output_test_hybrd_O2_$COMPILER.txt
./test_hybrj_O2_$COMPILER > output_test_hybrj_O2_$COMPILER.txt
./test_lmder_O2_$COMPILER > output_test_lmder_O2_$COMPILER.txt
./test_lmstr_O2_$COMPILER > output_test_lmstr_O2_$COMPILER.txt
./test_lmdif_O2_$COMPILER > output_test_lmdif_O2_$COMPILER.txt
./test_chkder_O2_$COMPILER > output_test_chkder_O2_$COMPILER.txt