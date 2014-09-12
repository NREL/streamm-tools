#Makefile for OpenCubMan

PROJECT= OPENCUBMAN
SOURCES=  fetocv.cpp cvtofe.cpp opencubman.cpp
OBJECTS=  fetocv.o opencubman.o -lm 
CC= g++
#CC= icpc


$(PROJECT): $(OBJECTS)
#	$(CC) -O2 -g  $(OBJECTS) -o $(PROJECT)
	$(CC) -O2 -g  cvtofe.o opencubman.o -lm -o cvtofe
	$(CC) -O2 -g  fetocv.o opencubman.o -lm -o fetocv
	@echo Compilation Complete

$(OBJECTS): $(SOURCES)
	@echo Compiling Sources
	$(CC) -O2 -fconserve-space -g -c $(SOURCES)
# this will be used once the code fulfills ansi
#        $(CC) -O2 -Wall -ansi -pedantic -fconserve-space -g -c $(SOURCES)
clean:
	@echo Deleting up $(OBJECTS) $(PROJECT)
	rm -f *.o;rm $(PROJECT)
									
