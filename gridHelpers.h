#ifndef GRIDHELPERS_H
#define GRIDHELPERS_H

#include <vector>
#include <algorithm>
#include "kfusion.h"

void addNeighbours(Image<bool, HostDevice> & gridWroteOn, std::vector<int> & vCubePositionsToClear, int coordX, int coordY, int coordZ, int resX, int resY, int resZ);
float medianDepth(const Image<uint16_t> & rawDepth, const int2 curPos, int radius);
float medianEstimation(std::vector<float> vLastTenDepths);
void updateLastTenDepths(std::vector<float> & vLastTenDepths, float medDepth);
#endif // GRIDHELPERS_H