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