CC          = c++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   = -g -Wno-deprecated

#GL_LIB = -framework OpenGL -framework GLUT -framework foundation
GL_LIB = -lGL -lglut -lGLU

#-----------------------------------------

TARGETS = stub

#OBJECTS = kdTree.o

#-----------------------------------------

LIBS = 
INCS = -I/usr/local/include/eigen3 -I/usr/include/eigen3

CCOPTS = $(OPT) $(INCS) 
LDOPTS = $(OPT) $(INCS) $(LIBS) 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)


clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------


stub: $(OBJECTS) stub.cpp
	$(CC) $(OBJECTS) stub.cpp $(LDOPTS) -o stub

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















