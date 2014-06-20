#!/bin/bash

make clean && make && sudo ldconfig /usr/local/cuda/lib64 && ./kinect
