#-*- mode: makefile; mode: font-lock; vc-back-end: RCS -*-

SHELL = /bin/bash

TARGET=ising_gen
BIN=../bin

# compilers etc
CC=gcc
CFLAGS=-O3 -pg -g
LD=gcc
LDFLAGS=-lm -O3 -pg -g

# List of object files needed to build program
OBJECTS=isinglib2.o wangmpi.o wang1D.o ising_gen.o ising_error.o ising_temp.o ising_mag.o ising_dbg.o WL-post.o

all: wang wangmpi ising_gen

dosmanip: dos-manip.o isinglib2.o
	$(LD) -o $(BIN)/dos-manip dos-manip.o isinglib2.o $(LDFLAGS)
ising_gen: isinglib2.o ising_gen.o
	$(LD) -o $(BIN)/ising_gen ising_gen.o isinglib2.o $(LDFLAGS)

ising_dbg: isinglib2.o ising_dbg.o
	$(LD) -o $(BIN)/ising_dbg ising_dbg.o isinglib2.o $(LDFLAGS)

wang1D: isinglib2.o wang1D.o
	$(LD) -o $(BIN)/wang1D wang1D.o isinglib2.o $(LDFLAGS)

WL-post: isinglib2.o WL-post.o
	$(LD) -o $(BIN)/WL-post WL-post.o isinglib2.o $(LDFLAGS)

wang: isinglib2.o wang.o
	$(LD) -o $(BIN)/wang wang.o isinglib2.o $(LDFLAGS)

jar_sweep: isinglib2.o jar_sweep.o
	$(LD) -o $(BIN)/jar_sweep jar_sweep.o isinglib2.o $(LDFLAGS)
	
ti: isinglib2.o ti.o
	$(LD) -o $(BIN)/ti ti.o isinglib2.o $(LDFLAGS)
		
wang3d: isinglib2.o wang3d.o
	$(LD) -o $(BIN)/wang3d wang3d.o isinglib2.o $(LDFLAGS)
	
wangmpi: isinglib2.o wangmpi.o
#	module load gnu/ompi
	mpicc -o $(BIN)/wangmpi wangmpi.o isinglib2.o $(LDFLAGS)
	
replica: isinglib2.o replica.o
#	module load gnu/ompi
	mpicc -o $(BIN)/replica replica.o isinglib2.o $(LDFLAGS)
	
ising_temp: isinglib2.o ising_temp.o
	$(LD) -o $(BIN)/ising_temp isinglib2.o ising_temp.o $(LDFLAGS)
	
ising_mag: isinglib2.o ising_mag.o
	$(LD) -o $(BIN)/ising_mag isinglib2.o ising_mag.o $(LDFLAGS)

ising_mag_vert_wolff: isinglib2.o ising_mag_vert_wolff.o
	$(LD) -o $(BIN)/ising_mag_vert_wolff isinglib2.o ising_mag_vert_wolff.o $(LDFLAGS)

ising_error:  isinglib2.o ising_error.o
	$(LD) -o $(BIN)/ising_error isinglib2.o ising_error.o $(LDFLAGS)
clean:
	rm $(OBJECTS) $(TARGET) $(FILES)

wangmpi.o: wangmpi.c
#	module load gnu/ompi
	mpicc -c wangmpi.c -o wangmpi.o
	
replica.o: replica.c
#	module load gnu/ompi
	mpicc -c replica.c -o replica.o

.PRECIOUS: %.o
.PHONY:  clean

%: %.o
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
