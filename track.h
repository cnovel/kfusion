#ifndef TRACK_H
#define TRACK_H

#include "kfusion.h"
#include "perfstats.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <GL/glut.h>
#include <GL/gl.h>
#include <cuda_gl_interop.h> // includes cuda_gl_interop.h

float3 learnColor(Image<uchar3> & texture);
float3 hsvOf(const Image<uchar3> & texture, uint2 pos);
void computeMap(const Image<uchar3> & texture, char map[50*50], const float3 hsvToTrack, const int2 pos);
void computeHeatmap(char map[50*50], int heatmap[50*50], int radius);
std::vector<uint2> findHeatPoints(int heatmap[50*50], int radius);
int2 newPosition(const std::vector<uint2>& vHeatPoints, int2 curPos);
void drawPosition(Image<uchar3> & texture, const int2 curPos);

// float medianDepth(Image<float> rawDepth, const int2 curPos, int radius);
float3 rotateVec( const Matrix4 & M, const float3 & v);
float4 multiply(const Matrix4 & M, const float4 & v);

float distanceBtwPoints(const float3 & pt1, const float3 & pt2);


#endif // TRACK_H