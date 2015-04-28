#-finline-limit=10000 -ffast-math -mcpu=athlon -m3dnow -fomit-frame-pointer 

LIBRARIES = -lm
OBJECTS = p480.o nrutil.o physics.o
HEADERS = p480.h nrutil.h physics.h
#OPTS = -pg -O3

#clean: 
	#rm *.o

rk4debug: rk4debug.c $(OBJECTS) $(HEADERS)
	gcc -o rk4debug rk4debug.c $(OBJECTS) $(LIBRARIES)

poisson: poisson.c $(OBJECTS) $(HEADERS)
	gcc -o poisson poisson.c $(OBJECTS) $(LIBRARIES)

quicksch: quicksch.c $(OBJECTS) $(HEADERS)
	gcc -o quicksch quicksch.c $(OBJECTS) $(LIBRARIES)

schdebug: schdebug.c $(OBJECTS) $(HEADERS)
	gcc -o schdebug schdebug.c $(OBJECTS) $(LIBRARIES)

Edebug: Edebug.c $(OBJECTS) $(HEADERS)
	gcc -o Edebug Edebug.c $(OBJECTS) $(LIBRARIES)

getEdebug: getEdebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getEdebug getEdebug.c $(OBJECTS) $(LIBRARIES)

getallEsdebug: getallEsdebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getallEsdebug getallEsdebug.c $(OBJECTS) $(LIBRARIES)

getPsidebug: getPsidebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getPsidebug getPsidebug.c $(OBJECTS) $(LIBRARIES)

getallPsisdebug: getallPsisdebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getallPsisdebug getallPsisdebug.c $(OBJECTS) $(LIBRARIES)

getRhodebug: getRhodebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getRhodebug getRhodebug.c $(OBJECTS) $(LIBRARIES)

getphidebug: getphidebug.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o getphidebug getphidebug.c $(OBJECTS) $(LIBRARIES)

test_excpPZ: test_excpPZ.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o test_excpPZ test_excpPZ.c $(OBJECTS) $(LIBRARIES)

test_getxc: test_getxc.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o test_getxc test_getxc.c $(OBJECTS) $(LIBRARIES)

test_atoms: test_atoms.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o test_atoms test_atoms.c $(OBJECTS) $(LIBRARIES)

test_getg: test_getg.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o test_getg test_getg.c $(OBJECTS) $(LIBRARIES)

test_atomC: test_atomC.c $(OBJECTS) $(HEADERS)
	gcc  $(OPTS) -o test_atomC test_atomC.c $(OBJECTS) $(LIBRARIES)

test_atomSb: test_atomSb.c $(OBJECTS) $(HEADERS)
	gcc $(OPTS) -o test_atomSb test_atomSb.c $(OBJECTS) $(LIBRARIES)

test_lu: test_lu.c $(OBJECTS) $(HEADERS)
	gcc  $(OPTS) -o test_lu test_lu.c $(OBJECTS) $(LIBRARIES)

test_it: test_it.c $(OBJECTS) $(HEADERS)
	gcc  $(OPTS) -o test_it test_it.c $(OBJECTS) $(LIBRARIES)

crap: crap.c $(OBJECTS) $(HEADERS)
	gcc  -o crap crap.c $(OBJECTS) $(LIBRARIES)
