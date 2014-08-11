#include "gridHelpers.h"

void addNeighbours(Image<bool, HostDevice> & gridWroteOn, std::vector<int> & vCubePositionsToClear, int coordX, int coordY, int coordZ, int resX, int resY, int resZ) {
	int min = -3;
	int max = 1-min;
	for(int i = min; i < max; i++) {
		for(int j = min; j < max; j++) {
			for(int k = min; k < max; k++) {
				if(0 <= coordX + i && coordX + i < resX && 0 <= coordY + j && coordY + j < resY && 0 <= coordZ + k && coordZ + k < resZ) {
					gridWroteOn[make_uint2(coordX+i+resX*(coordY+j)+resX*resY*(coordZ+k),0)] = true;
                    vCubePositionsToClear.push_back(coordX+i+resX*(coordY+j)+resX*resY*(coordZ+k));
				}
			}
		}
	}
}

float medianDepth(const Image<uint16_t> & rawDepth, const int2 curPos, int radius) {
    std::vector<float> vDepths;

    for(int i = -radius + curPos.x; i < radius + curPos.x; i++) {
        for(int j = -radius + curPos.y; j < radius + curPos.y; j++) {
            if(i >= 0 && i < rawDepth.size.x && j >= 0 && j < rawDepth.size.y) {
                if (rawDepth[make_uint2(i,j)] != 0) {
                    vDepths.push_back(rawDepth[make_uint2(i,j)]/1000.0f);
                }
            }
        }
    }

    std::sort(vDepths.begin(), vDepths.end());
    if(!vDepths.empty()) {
        int meanInt = vDepths.size()/2;
        return vDepths[meanInt]+0.01f;
    }
    else
        return -1.0f;
}

float medianEstimation(std::vector<float> vLastTenDepths) {
	std::sort(vLastTenDepths.begin(), vLastTenDepths.end());
	return vLastTenDepths[vLastTenDepths.size()/2];
}

void updateLastTenDepths(std::vector<float> & vLastTenDepths, float medDepth) {
	for(int i=0; i < vLastTenDepths.size()-1; i++) {
		vLastTenDepths[i] = vLastTenDepths[i+1];
	}
	vLastTenDepths[vLastTenDepths.size()-1] = medDepth;
}

Mat3 getRot(const Matrix4 & pose) {
    Mat3 rot;
    rot.data[0].x = pose.data[0].x;
    rot.data[0].y = pose.data[0].y;
    rot.data[0].z = pose.data[0].z;

    rot.data[1].x = pose.data[1].x;
    rot.data[1].y = pose.data[1].y;
    rot.data[1].z = pose.data[1].z;

    rot.data[2].x = pose.data[2].x;
    rot.data[2].y = pose.data[2].y;
    rot.data[2].z = pose.data[2].z;

    return rot;
}

Mat3 transpose(const Mat3 & rot) {
    Mat3 trans;
    trans.data[0].x = rot.data[0].x;
    trans.data[0].y = rot.data[1].x;
    trans.data[0].z = rot.data[2].x;

    trans.data[1].x = rot.data[0].y;
    trans.data[1].y = rot.data[1].y;
    trans.data[1].z = rot.data[2].y;

    trans.data[2].x = rot.data[0].z;
    trans.data[2].y = rot.data[1].z;
    trans.data[2].z = rot.data[2].z;

    return trans;
}

Mat3 multiply(const Mat3 & A, const Mat3 & B) {
    Mat3 res;

    res.data[0].x = A.data[0].x*B.data[0].x + A.data[0].y*B.data[1].x + A.data[0].z*B.data[2].x;
    res.data[0].y = A.data[0].x*B.data[0].y + A.data[0].y*B.data[1].y + A.data[0].z*B.data[2].y;
    res.data[0].z = A.data[0].x*B.data[0].z + A.data[0].y*B.data[1].z + A.data[0].z*B.data[2].z;

    res.data[1].x = A.data[1].x*B.data[0].x + A.data[1].y*B.data[1].x + A.data[1].z*B.data[2].x;
    res.data[1].y = A.data[1].x*B.data[0].y + A.data[1].y*B.data[1].y + A.data[1].z*B.data[2].y;
    res.data[1].z = A.data[1].x*B.data[0].z + A.data[1].y*B.data[1].z + A.data[1].z*B.data[2].z;

    res.data[2].x = A.data[2].x*B.data[0].x + A.data[2].y*B.data[1].x + A.data[2].z*B.data[2].x;
    res.data[2].y = A.data[2].x*B.data[0].y + A.data[2].y*B.data[1].y + A.data[2].z*B.data[2].y;
    res.data[2].z = A.data[2].x*B.data[0].z + A.data[2].y*B.data[1].z + A.data[2].z*B.data[2].z;

    return res;
}

void correctView(Matrix4 & ovrPose, const Mat3 & rotHMDtoOVR) {
    Mat3 rotated;
    rotated.data[0].x = rotHMDtoOVR.data[0].x*ovrPose.data[0].x + rotHMDtoOVR.data[0].y*ovrPose.data[1].x + rotHMDtoOVR.data[0].z*ovrPose.data[2].x;
    rotated.data[0].y = rotHMDtoOVR.data[0].x*ovrPose.data[0].y + rotHMDtoOVR.data[0].y*ovrPose.data[1].y + rotHMDtoOVR.data[0].z*ovrPose.data[2].y;
    rotated.data[0].z = rotHMDtoOVR.data[0].x*ovrPose.data[0].z + rotHMDtoOVR.data[0].y*ovrPose.data[1].z + rotHMDtoOVR.data[0].z*ovrPose.data[2].z;

    rotated.data[1].x = rotHMDtoOVR.data[1].x*ovrPose.data[0].x + rotHMDtoOVR.data[1].y*ovrPose.data[1].x + rotHMDtoOVR.data[1].z*ovrPose.data[2].x;
    rotated.data[1].y = rotHMDtoOVR.data[1].x*ovrPose.data[0].y + rotHMDtoOVR.data[1].y*ovrPose.data[1].y + rotHMDtoOVR.data[1].z*ovrPose.data[2].y;
    rotated.data[1].z = rotHMDtoOVR.data[1].x*ovrPose.data[0].z + rotHMDtoOVR.data[1].y*ovrPose.data[1].z + rotHMDtoOVR.data[1].z*ovrPose.data[2].z;

    rotated.data[2].x = rotHMDtoOVR.data[2].x*ovrPose.data[0].x + rotHMDtoOVR.data[2].y*ovrPose.data[1].x + rotHMDtoOVR.data[2].z*ovrPose.data[2].x;
    rotated.data[2].y = rotHMDtoOVR.data[2].x*ovrPose.data[0].y + rotHMDtoOVR.data[2].y*ovrPose.data[1].y + rotHMDtoOVR.data[2].z*ovrPose.data[2].y;
    rotated.data[2].z = rotHMDtoOVR.data[2].x*ovrPose.data[0].z + rotHMDtoOVR.data[2].y*ovrPose.data[1].z + rotHMDtoOVR.data[2].z*ovrPose.data[2].z;

    ovrPose.data[0].x = rotated.data[0].x;
    ovrPose.data[0].y = rotated.data[0].y;
    ovrPose.data[0].z = rotated.data[0].z;

    ovrPose.data[1].x = rotated.data[1].x;
    ovrPose.data[1].y = rotated.data[1].y;
    ovrPose.data[1].z = rotated.data[1].z;

    ovrPose.data[2].x = rotated.data[2].x;
    ovrPose.data[2].y = rotated.data[2].y;
    ovrPose.data[2].z = rotated.data[2].z;
}