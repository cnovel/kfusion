CXX = nvcc
CC = nvcc
CPPFLAGS=-I/usr/include -I/usr/local/cuda/include -DLIBFREENECT_INTERFACE
## -I../include/libfreenect  -I/opt/local/include
CXXFLAGS=-g -m64 -O3 -use_fast_math
## -ptx -src-in-ptx
LDFLAGS= -g -m64 -lfreenect -lglut -lGL -L/usr/local/cuda/lib64
## -g -m64 -L/lib -lfreenect
%.o: %.cu
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $^

all: kinect test

test: kfusion.o helpers.o test.o

kinect: kfusion.o helpers.o kinect.o interface.o

clean:
	rm *.o test kinect