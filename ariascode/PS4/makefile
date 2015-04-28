LIBRARIES = -lm
OBJECTS = p480.o nrutil.o physics.o
HEADERS = p480.h nrutil.h physics.h

rk4debug: rk4debug.c $(OBJECTS) $(HEADERS)
	cc -o rk4debug rk4debug.c $(OBJECTS) $(LIBRARIES)

poisson: poisson.c $(OBJECTS) $(HEADERS)
	cc -o poisson poisson.c $(OBJECTS) $(LIBRARIES)

quicksch: quicksch.c $(OBJECTS) $(HEADERS)
	cc -o quicksch quicksch.c $(OBJECTS) $(LIBRARIES)

schdebug: schdebug.c $(OBJECTS) $(HEADERS)
	cc -o schdebug schdebug.c $(OBJECTS) $(LIBRARIES)

Edebug: Edebug.c $(OBJECTS) $(HEADERS)
	cc -o Edebug Edebug.c $(OBJECTS) $(LIBRARIES)

getEdebug: getEdebug.c $(OBJECTS) $(HEADERS)
	cc -o getEdebug getEdebug.c $(OBJECTS) $(LIBRARIES)

getallEsdebug: getallEsdebug.c $(OBJECTS) $(HEADERS)
	cc -o getallEsdebug getallEsdebug.c $(OBJECTS) $(LIBRARIES)

getPsidebug: getPsidebug.c $(OBJECTS) $(HEADERS)
	cc -o getPsidebug getPsidebug.c $(OBJECTS) $(LIBRARIES)

getallPsisdebug: getallPsisdebug.c $(OBJECTS) $(HEADERS)
	cc -o getallPsisdebug getallPsisdebug.c $(OBJECTS) $(LIBRARIES)

getRhodebug: getRhodebug.c $(OBJECTS) $(HEADERS)
	cc -o getRhodebug getRhodebug.c $(OBJECTS) $(LIBRARIES)

getphidebug: getphidebug.c $(OBJECTS) $(HEADERS)
	cc -o getphidebug getphidebug.c $(OBJECTS) $(LIBRARIES)
