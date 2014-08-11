#ifndef GRIDHELPERS_H
#define GRIDHELPERS_H

#include <vector>
#include <algorithm>
#include "kfusion.h"

void addNeighbours(Image<bool, HostDevice> & gridWroteOn, std::vector<int> & vCubePositionsToClear, int coordX, int coordY, int coordZ, int resX, int resY, int resZ);
float medianDepth(const Image<uint16_t> & rawDepth, const int2 curPos, int radius);
float medianEstimation(std::vector<float> vLastTenDepths);
void updateLastTenDepths(std::vector<float> & vLastTenDepths, float medDepth);

Mat3 getRot(const Matrix4 & pose);
Mat3 transpose(const Mat3 & rot);
Mat3 multiply(const Mat3 & A, const Mat3 & B);
void correctView(Matrix4 & ovrPose, const Mat3 & rotHMDtoOVR);
#endif // GRIDHELPERS_H