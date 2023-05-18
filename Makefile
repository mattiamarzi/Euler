FC = gfortran
FLAGS = -fmax-errors=3 -O3
TARGET = euler

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

%.o: %.f90
	$(FC) $(FLAGS) -c $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)

common.o: precision.o
main.o: precision.o common.o 

clean:
	rm *.o *.mod *~ *.dat $(TARGET)
