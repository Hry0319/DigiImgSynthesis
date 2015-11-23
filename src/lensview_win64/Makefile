ARCH = $(shell uname)

CFLAGS = -g -Wall -O2 -funroll-loops
CC = gcc
CXX = g++
LDFLAGS = -L/usr/lib -L/usr/X11R6/lib -L../core

ifeq ($(ARCH), Darwin)
LIBS = -framework GLUT -framework OpenGL -lobjc -lpbrt 
else
LIBS = -lglut -lGLU -lXmu -lXi -lpbrt -lm 
endif

INCLUDE = -I. -I../core

TARGET = lensview
SRCS = lensview.cpp lensdata.cpp

$(TARGET) : $(SRCS)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(INCLUDE) $(SRCS) -o $@ $(LIBS)

clean:
	rm -f $(TARGET) *.o a.out *~ 


