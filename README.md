
OFusion 0.1
=============

OFusion is based on KFusion 0.4, copyright TU Graz, Gerhard Reitmayr, 2011 - 2013

OFusion uses the KFusion reconstruction and displays it into the Oculus Rift.


Requirements
------------
OFusion has been tested only with Linux Ubuntu 13.10 64bits.

OFusion depends mainly on the following libraries:

* TooN http://www.edwardrosten.com/cvd/toon.html
* GLUT
* libfreenect http://openkinect.org/
* CUDA 5 SDK http://developer.nvidia.com/cuda
* Oculus Rift SDK

Install
-----
Tweak the Makefile for your system. I tried to be as clear as possible. Then execute **launch.sh**. It will compile and launch the program. Please check that your Kinect and Rift are plugged in.


To do
----------
- ~~Deformation of the 3D image for display in the Rift~~ done
- Calibration of the Rift and the Kinect
- Add augmented reality assets

Contributors
-------------
Hartmut Seichter, Gerhard Reitmayr

Thanks to
-------------
Andrew Davison, Steven Lovegrove
