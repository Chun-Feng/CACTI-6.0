LIBS = -lm
#FLAGS = -Wall -O2 -fomit-frame-pointer -msse -march=athlon
#FLAGS = -Wall -O0
CC = gcc -g -pthread
#CC = gcc -g -pthread -pg

SRCS = main.c time.c area.c io.c leakage.c technology.c basic_circuit.c def.h areadef.h leakage.h basic_circuit.h io.h time.h cacti_interface.h router.h router.c

OBJS = main.o time.o area.o io.o leakage.o technology.o basic_circuit.o router.o

all: cacti

pythonlib : time.o area.o io.o leakage.o technology.o basic_circuit.o cacti_wrap.o router.o
		gcc -shared $(FLAGS) area.o time.o leakage.o technology.o basic_circuit.o io.o cacti_wrap.o router.o -L /usr/lib/python2.4/config -lpython2.4 -o _cacti.o

cacti : main.o time.o area.o io.o leakage.o technology.o basic_circuit.o router.o
	  $(CC) $(FLAGS) $(OBJS) -o cacti $(LIBS)

main.o : main.c def.h areadef.h leakage.h basic_circuit.h io.h technology.c
	  $(CC) $(FLAGS) -c main.c -o main.o

leakage.o : leakage.h leakage.c
	  $(CC) $(FLAGS) -c leakage.c -o leakage.o

technology.o : def.h areadef.h technology.c
	  $(CC) $(FLAGS) -c technology.c -o technology.o

time.o :  time.c def.h areadef.h leakage.h basic_circuit.h cacti_interface.h
	   $(CC) $(FLAGS) -c time.c -o time.o

area.o : area.c def.h areadef.h cacti_interface.h
	   $(CC) $(FLAGS) -c area.c -o area.o 

io.o : def.h io.c areadef.h cacti_interface.h router.h io.h
	  $(CC) $(FLAGS) -c io.c -o io.o

router.o : router.c router.h cacti_interface.h def.h
	 $(CC) $(FLAGS) -c router.c -o router.o

basic_circuit.o : basic_circuit.h basic_circuit.c
		   $(CC) $(FLAGS) -c basic_circuit.c -o basic_circuit.o 

cacti_wrap.o :  cacti_wrap.c
		$(CC) -c io.c area.c time.c \
		basic_circuit.c leakage.c cacti_wrap.c \
		-I /usr/include/python2.4 \
		-I /usr/lib/python2.4/config

cacti_wrap.c: cacti.i
			swig -python cacti.i

clean:
	  rm -f *.o cacti cache_params.aux core

