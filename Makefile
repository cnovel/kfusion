#####
# MAIN PROGRAM
#####

CXX = nvcc

CPPFLAGS=-I/usr/include -I/usr/local/cuda/include -DLIBFREENECT_INTERFACE

OBJ = kfusion.o helpers.o test.o interface.o kinect.o

CXXFLAGS= -g -m64 -O2 -use_fast_math -L. -L/usr/lib/x86_64-linux-gnu -L/usr/local/cuda/lib64

LIBS = -lGL -lGLU -lglut -lX11 -lpthread -lglfw3 -ludev -lXxf86vm -lXrandr -lXi -lfreenect libovr.a

#####
# RULES
#####

all: kinect test

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBS) -c $^ 

%.o: %.cu
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LIBS) -c $^ 

test: kfusion.o helpers.o test.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	chmod +x $@


kinect: kfusion.o helpers.o kinect.o interface.o track.o gridHelpers.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	chmod +x $@


clean:
	rm -f $(OBJ)
	rm -f test kinect