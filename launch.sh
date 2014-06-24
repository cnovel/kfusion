#!/bin/bash

make clean && make -j 8 && sudo ldconfig /usr/local/cuda/lib64 && ./kinect
