#####
# MAIN PROGRAM
#####

CXX = nvcc
CC = nvcc
CPPFLAGS=-I/usr/include -I/usr/local/cuda/include -DLIBFREENECT_INTERFACE

OBJ = kfusion.o helpers.o test.o interface.o kinect.o

CXXFLAGS=-g -m64 -O3 -use_fast_math -L. -L/usr/lib/x86_64-linux-gnu/

LIBS = libovr.a -lGL -lGLU -lglut -lX11 -lpthread -lglfw3 -ludev -lXxf86vm -lXrandr -lXi

LDFLAGS= -g -m64 -lfreenect -lglut -lGL -L/usr/local/cuda/lib64

#####
# RULES
#####

kinect.o: kinect.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBS) -c $^ 

%.o: %.cu
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBS) -c $^ 

all: kinect test

test: kfusion.o helpers.o test.o

kinect: kfusion.o helpers.o kinect.o interface.o

clean:
	rm -f $(OBJ)
	rm -f test kinect