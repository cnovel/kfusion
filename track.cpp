#include "track.h"

float3 learnColor(Image<uchar3> & texture) {
	std::vector<float> vHue, vSat, vVal;
    for (int i = 310; i < 331; i++) {
        for (int j = 230; j < 251; j++) {
        	uchar3 pixColor = texture[make_uint2(i,j)];
        	float r = ((float)pixColor.x) / 255.0f;
    		float g = ((float)pixColor.y) / 255.0f;
    		float b = ((float)pixColor.z) / 255.0f;

    		float hue, sat, val;
    		float max = (r > g) ? (r > b) ? r : b : (g > b) ? g : b;
    		float min = (r < g) ? (r < b) ? r : b : (g < b) ? g : b;

    		if (max == min) {
    		    hue = 0;
    		} else if (max == r) {
    		    hue = (60.0f*(g-b)/(max-min) + 360.0f);
    		    if (hue > 360.0f) {
    		        hue -= 360.0f;
    		    }
    		} else if (max == g){
    		    hue = (60.0f*(b-r)/(max-min) + 120.0f);
    		} else {
    		    hue = (60.0f*(r-g)/(max-min) + 240.0f);
    		}

    		hue = hue*3.14159/360.0f;

    		if (max == 0)
    		    sat = 0;
    		else
    		    sat = 1.0f - min/max;

    		val = max;

    		vHue.push_back(hue);
    		vSat.push_back(sat);
    		vVal.push_back(val);
        }
    }
    float sumSin, sumCos, meanHue, meanSat, meanVal;
    sumCos = sumSin = meanVal = meanSat = meanHue = 0.0f;
    for (int i = 0; i < vHue.size(); i++) {
    	sumCos += cos(vHue[i]);
    	sumSin += sin(vHue[i]);
    	meanSat += vSat[i];
    	meanVal += vVal[i];
    }
    meanVal /= vHue.size();
    meanSat /= vHue.size();
    sumCos /= vHue.size();
    sumSin /= vHue.size();

    meanHue = atan2(sumSin, sumCos);

    float3 hsv = make_float3(meanHue, meanSat, meanVal);
    return hsv;
}

float3 hsvOf(const Image<uchar3> & texture, uint2 pos) {
	uchar3 pixColor = texture[pos];
	float r = ((float)pixColor.x) / 255.0f;
	float g = ((float)pixColor.y) / 255.0f;
	float b = ((float)pixColor.z) / 255.0f;

	float epsiSat = 0.2f;
	float epsiHue = 0.2f;

	float hue, sat, val;
	float max = (r > g) ? (r > b) ? r : b : (g > b) ? g : b;
	float min = (r < g) ? (r < b) ? r : b : (g < b) ? g : b;

	if (max == min) {
		hue = 0;
	} else if (max == r) {
		hue = (60.0f*(g-b)/(max-min) + 360.0f);
		if (hue > 360.0f) {
			hue -= 360.0f;
		}
	} else if (max == g){
		hue = (60.0f*(b-r)/(max-min) + 120.0f);
	} else {
		hue = (60.0f*(r-g)/(max-min) + 240.0f);
	}

	hue = hue*3.14159/360.0f;

	if (max == 0)
		sat = 0;
	else
		sat = 1.0f - min/max;

	val = max;

	return make_float3(hue, sat, val);
}

void computeMap(const Image<uchar3> & texture, char map[50*50], const float3 hsvToTrack, const int2 pos) {
	float epsiHue = .2f;

	for(int i = 0; i < 50; i++) {
		for (int j = 0; j < 50; j++) {
			map[i*50+j] = 0;
			float diffSat, diffHue;
			int2 curPos = make_int2(pos.x + i - 25, pos.y + j - 25);
			if (curPos.x >= 0 && curPos.y >= 0 && curPos.x < 640 && curPos.y < 480) {
				float3 hsv = hsvOf(texture, make_uint2(curPos.x, curPos.y));

				diffSat = (hsv.y > hsvToTrack.y) ? hsv.y - hsvToTrack.y : hsvToTrack.y - hsv.y;
				diffHue = (hsv.x + 3.14159 - hsvToTrack.x);
				if (diffHue > 3.14159)
					diffHue -= 3.14159;

				if (diffHue < epsiHue && diffSat < hsvToTrack.y/4.0f) {
					map[i*50+j]++;
				}
			}
		}
	}
}

void computeHeatmap(char map[50*50], int heatmap[50*50], int radius) {
	for(int i = 0; i < 50; i++) {
		for(int j = 0; j < 50; j++) {
			heatmap[i*50+j] = 0;
			for(int k = i-radius; k < i+radius+1; k++) {
				for(int l = j-radius; l < j+radius+1; l++) {
					if (k>=0 && k<50 && l>=0 && l<50)
						heatmap[i*50+j] += map[k*50+l];
				}
			}
		}
	}
}

std::vector<uint2> findHeatPoints(int heatmap[50*50], int radius) {
	std::vector<uint2> vHeatPoints;
	for(int i = 0; i < 50; i++) {
		for(int j = 0; j < 50; j++) {
			bool heat = true;
			if (heatmap[i*50+j] < 10)
				heat = false;
			for(int k = i-radius; k < i+radius+1; k++) {
				for(int l = j-radius; l < j+radius+1; l++) {
					if (heat && k>=0 && k<50 && l>=0 && l<50) {
						if (heatmap[i*50+j] < heatmap[k*50+l])
							heat = false;
					}
				}
			}
			if (heat) {
				uint2 pos = make_uint2(i,j);
				vHeatPoints.push_back(pos);
			}
		}
	}
	return vHeatPoints;
}

int2 newPosition(const std::vector<uint2>& vHeatPoints, int2 curPos) {
	float dist = 50.0f*50.0f;
	int bestHP = -1;
	for(int i = 0; i < vHeatPoints.size(); i++) {
		float curDist = (25.0 - vHeatPoints[i].x)*(25.0 - vHeatPoints[i].x) + (25.0 - vHeatPoints[i].y)*(25.0 - vHeatPoints[i].y);
		if (curDist < dist) {
			dist = curDist;
			bestHP = i;
		}
	}
	int2 newPos;
	if(bestHP != -1) {
		newPos = make_int2(curPos.x - 25 + vHeatPoints[bestHP].x, curPos.y - 25 + vHeatPoints[bestHP].y);
	} else {
		newPos = make_int2(-1,-1);
		//newPos = curPos;
	}
	return newPos;
}

void drawPosition(Image<uchar3> & texture, const int2 curPos) {
	for(int i = 0; i < 10; i++) {
		for(int j = 0; j < 10; j++) {
			int2 drawOn = make_int2(curPos.x - 5 + i, curPos.y -5 + j);
			if (drawOn.x >= 0 && drawOn.y >= 0 && drawOn.x < 640 && drawOn.y < 480) {
				texture[make_uint2(drawOn.x, drawOn.y)] = make_uchar3(255,0,255);
			}
		}
	}
}

// float medianDepth(const Image<float> & rawDepth, const int2 curPos, int radius) {
// 	std::vector<float> vDepths;
// 	for(int i = -radius + curPos.x; i < radius + curPos.x; i++) {
// 		for(int j = -radius + curPos.y; j < radius + curPos.y; j++) {
// 			if(i >= 0 && i < 640 && j >= 0 && j < 480) {
// 				if (rawDepth[make_uint2(i,j)] != 0)
// 					vDepths.push_back(rawDepth[make_uint2(i,j)]);
// 			}
// 		}
// 	}

// 	std::sort(vDepths.begin(), vDepths.end());
// 	if(!vDepths.empty())
// 		return vDepths[vDepths.size()/2];
// 	else
// 		return -1.0f;
// }

float3 rotateVec( const Matrix4 & M, const float3 & v){
    return make_float3(
        dot(make_float3(M.data[0]), v),
        dot(make_float3(M.data[1]), v),
        dot(make_float3(M.data[2]), v));
}

float4 multiply(const Matrix4 & M, const float4 & v) {
	float x, y, z, w;
	x = M.data[0].x * v.x + M.data[0].y * v.y + M.data[0].z * v.z + M.data[0].w * v.w;
	y = M.data[1].x * v.x + M.data[1].y * v.y + M.data[1].z * v.z + M.data[1].w * v.w;
	z = M.data[2].x * v.x + M.data[2].y * v.y + M.data[2].z * v.z + M.data[2].w * v.w;
	w = M.data[3].x * v.x + M.data[3].y * v.y + M.data[3].z * v.z + M.data[3].w * v.w;
	return make_float4(x,y,z,w);
}

float distanceBtwPoints(const float3 & pt1, const float3 & pt2) {
	return (pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y) + (pt1.z - pt2.z)*(pt1.z - pt2.z);
}