/*
Copyright (c) 2011-2013 Gerhard Reitmayr, TU Graz

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "kfusion.h"
#include "perfstats.h"
#include <iostream>
#include <GL/glut.h>
#include <GL/gl.h>

using namespace std;

PerfStats Stats;

__global__ void setSphere( Volume volume, const float3 center, const float radius, const float val ){
    uint3 pos = make_uint3(thr2pos2());
    for(pos.z = 0; pos.z < volume.size.z; ++pos.z) {
        const float d = length(volume.pos(pos) - center);
        if(d < radius)
            volume.set(pos, make_float2(val, 0.0f));
    }
}

__global__ void setBox( Volume volume, const float3 min_corner, const float3 max_corner, const float val ){
    uint3 pos = make_uint3(thr2pos2());
    for(pos.z = 0; pos.z < volume.size.z; ++pos.z) {
        const float3 p = volume.pos(pos);
        if(min_corner.x < p.x && min_corner.y < p.y && min_corner.z < p.z && 
           p.x < max_corner.x && p.y < max_corner.y && p.z < max_corner.z )
            volume.set(pos, make_float2(val, 0.0f));
    }
}

void initVolumeWrap( Volume volume, const float val ){
    dim3 block(32,16);
    initVolume<<<divup(dim3(volume.size.x, volume.size.y), block), block>>>(volume, make_float2(val, 0.0f));
}

void setBoxWrap(Volume volume, const float3 min_corner, const float3 max_corner, const float val ){
    dim3 block(32,16);
    setBox<<<divup(dim3(volume.size.x, volume.size.y), block), block>>>(volume, min_corner, max_corner, val);
}

void setSphereWrap(Volume volume, const float3 center, const float radius, const float val ){
    dim3 block(32,16);
    setSphere<<<divup(dim3(volume.size.x, volume.size.y), block), block>>>(volume, center, radius, val);
}

__global__ void renderNormals( Image<uchar3> out, const Image<float3> in ){
    float3 n = in.el();
    if(n.x == -2)
        out.el() = make_uchar3(0,0,0);
    else {
        n = normalize(n);
        out.el() = make_uchar3(n.x*128 + 128, n.y*128+128, n.z*128+128);
    }
}

void renderNormalMap( Image<uchar3> out, const Image<float3> & normal ){
    dim3 block(20,20);
    renderNormals<<<divup(normal.size, block), block>>>( out, normal );
}

__global__ void renderLightKernel( Image<uchar4> out, const Image<float3> vertex, const Image<float3> normal, const float3 light, const float3 ambient ){
    if(normal.el().x == -2.0f)
        out.el() = make_uchar4(0,0,0,255);
    else {
        const float3 diff = normalize(light - vertex.el());
        const float dir = fmaxf(dot(normal.el(), diff), 0.f);
        const float3 col = clamp(make_float3(dir) + ambient, 0.f, 1.f) * 255;
        out.el() = make_uchar4(col.x, col.y, col.z, 255);
    }
}

void renderLight( Image<uchar4> out, const Image<float3> & vertex, const Image<float3> & normal, const float3 light, const float3 ambient ){
    dim3 block(32,16);
    renderLightKernel<<<divup(out.size, block), block>>>( out, vertex, normal, light, ambient );
}

__global__ void renderTextureKernel( Image<uchar4> out, const Image<float3> vertex, const Image<float3> normal, const Image<uchar3> texture, const Matrix4 texproj, const float3 light){
    if(normal.el().x == -2.0f)
        out.el() = make_uchar4(0,0,0,255);
    else {
        const float3 proj = texproj * vertex.el();
        const float2 projPixel = make_float2( proj.x / proj.z + 0.5f, proj.y / proj.z + 0.5f);
        
        const float3 diff = normalize(light - vertex.el());
        const float dir = fmaxf(dot(normal.el(), diff), 0.f); // * 255;
        if(projPixel.x < 0 || projPixel.x > texture.size.x-1 || projPixel.y < 0 || projPixel.y > texture.size.y-1 ){
            out.el() = make_uchar4(dir*255,dir*255,dir*255,255);
        } else {
            const uchar3 texcol = texture[make_uint2(projPixel.x, projPixel.y)];
            out.el() = make_uchar4(texcol.x*dir, texcol.y*dir, texcol.z*dir, 255);
        }
    }
}

void renderTexture( Image<uchar4> out, const Image<float3> & vertex, const Image<float3> & normal, const Image<uchar3> & texture, const Matrix4 & texproj, const float3 light){
    dim3 block(32,16);
    renderTextureKernel<<<divup(out.size, block), block>>>( out, vertex, normal, texture, texproj, light);
}

__global__ void renderDepth( Image<uchar3> out, const Image<float> depth, const float nearPlane, const float farPlane){
    const float d = (clamp(depth.el(), nearPlane, farPlane) - nearPlane) / (farPlane - nearPlane);
    out.el() = make_uchar3(d * 255, d * 255, d * 255);
}

void renderDepthMap( Image<uchar3> out, const Image<float> & depth, const float nearPlane, const float farPlane ){
    dim3 block(32,16);
    renderDepth<<<divup(depth.size, block), block>>>( out, depth, nearPlane, farPlane );
}

__global__ void renderTrack( Image<uchar4> out, const Image<TrackData> data ){
    const uint2 pos = thr2pos2();
    switch(data[pos].result){
    case 1: out[pos] = make_uchar4(128, 128, 128,0);  // ok
        break;
    case -1: out[pos] = make_uchar4(0, 0, 0,0);      // no input
        break;
    case -2: out[pos] = make_uchar4(255,0,0,0);        // not in image
        break;
    case -3:  out[pos] = make_uchar4(0,255,0,0);        // no correspondence
        break;
    case -4: out[pos] = make_uchar4(0,0,255,0);        // to far away
        break;
    case -5: out[pos] = make_uchar4(255,255,0,0);     // wrong normal
        break;
    }
}

void renderTrackResult( Image<uchar4> out, const Image<TrackData> & data ){
    dim3 block(32,16);
    renderTrack<<<divup(out.size, block), block>>>( out, data );
}

__global__ void raycastLight( Image<uchar4> render, const Volume volume, const Matrix4 view, const float nearPlane, const float farPlane, const float step, const float largestep, const float3 light, const float3 ambient){
    const uint2 pos = thr2pos2();
    int2 posS = make_int2(pos.x, pos.y);
    
    float4 hit = raycast( volume, posS, view, nearPlane, farPlane, step, largestep);
    if(hit.w > 0){
        const float3 test = make_float3(hit);
        const float3 surfNorm = volume.grad(test);
        if(length(surfNorm) > 0){
            const float3 diff = normalize(light - test);
            const float dir = fmaxf(dot(normalize(surfNorm), diff), 0.f);
            const float3 col = clamp(make_float3(dir) + ambient, 0.f, 1.f) * 255;
            render.el() = make_uchar4(col.x, col.y, col.z,0);
        } else {
            render.el() = make_uchar4(0,0,0,0);
        }
    } else {
        render.el() = make_uchar4(0,0,0,0);
    }
}

void renderVolumeLight( Image<uchar4> out, const Volume & volume, const Matrix4 view, const float nearPlane, const float farPlane, const float largestep, const float3 light, const float3 ambient ){
    dim3 block(16,16);
    raycastLight<<<divup(out.size, block), block>>>( out,  volume, view, nearPlane, farPlane, volume.dim.x/volume.size.x, largestep, light, ambient );
}

__global__ void raycastInput( Image<float3> pos3D, Image<float3> normal, Image<float> depth, const Volume volume, const Matrix4 view, const float nearPlane, const float farPlane, const float step, const float largestep, const int2 outputSize){
    const uint2 pos = thr2pos2();

    int2 transPos = make_int2(pos.x, pos.y);
    transPos.x -= (outputSize.x - 640) / 2;
    transPos.y -= (outputSize.y - 480) / 2;
    
    float4 hit = raycast( volume, transPos, view, nearPlane, farPlane, step, largestep);
    if(hit.w > 0){
        pos3D[pos] = make_float3(hit);
        depth[pos] = hit.w;
        float3 surfNorm = volume.grad(make_float3(hit));
        if(length(surfNorm) == 0){
            normal[pos].x = -2;
        } else {
            normal[pos] = normalize(surfNorm);
        }
    } else {
        pos3D[pos] = make_float3(0);
        normal[pos] = make_float3(0);
        depth[pos] = 0;
    }
}

void renderInput( Image<float3> pos3D, Image<float3> normal, Image<float> depth, const Volume volume, const Matrix4 view, const float nearPlane, const float farPlane, const float step, const float largestep, const int2 outputSize){
    dim3 block(32,16);
    raycastInput<<<divup(pos3D.size, block), block>>>(pos3D, normal, depth, volume, view, nearPlane, farPlane, step, largestep, outputSize);
}

__global__ void OculusCam(Image<uchar4> out, const Volume volume, const Image<uchar3> texture, const Matrix4 view, const Matrix4 texproj, const float nearPlane, const float farPlane, const float step, const float largestep, const float3 light, const float3 ambient) {
    float f = 204.4f;
    float2 iResolution = make_float2(1280, 800);
    float DistortionScale = 1.71461;
    float pp_adjust = 48.62;
    //float IPD_by2 = 0.01;
    float3 K = make_float3(1.00f, 0.22f, 0.24f);
    //mat4 iT_bh = mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

    const uint2 pos = thr2pos2();
    float2 uv = make_float2(pos.x, pos.y);
    
    // Which eye?
    float i = uv.x <= iResolution.x/2.0 ? 0.0 : 1.0;

    // Compute Principle point
    float sx = i*2.0 - 1.0;
    float2 pp = make_float2( sx*(iResolution.x / 4.0 - pp_adjust) + iResolution.x / 2.0, iResolution.y / 2.0 );

    // Distort uv for Oculus (using res independant coords)
    float2 theta = make_float2(4.0 * (uv.x-pp.x) / iResolution.x, 4.0 * (uv.y-pp.y) / iResolution.x);
    float rSq= theta.x * theta.x + theta.y * theta.y;
    uv.x = pp.x + iResolution.x * theta.x * (K.x + rSq*(K.y + rSq*K.z ) ) / (4.0 * DistortionScale);
    uv.y = pp.y + iResolution.x * theta.y * (K.x + rSq*(K.y + rSq*K.z ) ) / (4.0 * DistortionScale);

    float2 uvMinusPp = make_float2(uv.x - pp.x, uv.y - pp.y);
    uvMinusPp.x /= f;
    uvMinusPp.y /= f;

    float3 ray = make_float3(uvMinusPp.x, uvMinusPp.y, 1.0f);
    float rayNorm = sqrt(ray.x*ray.x + ray.y*ray.y + ray.z*ray.z);
    ray.x /= rayNorm;
    ray.y /= rayNorm;
    ray.z /= rayNorm;

    float4 hit;
    float3 pos3D;
    float3 normal;

    hit = raycastDirPos(volume, view, nearPlane, farPlane, step, largestep, ray);

    if(hit.w > 0){
        pos3D = make_float3(hit);
        float3 surfNorm = volume.grad(make_float3(hit));
        if(length(surfNorm) == 0){
            normal.x = -2;
        } else {
            normal = normalize(surfNorm);
        }
    } else {
        pos3D = make_float3(0);
        normal = make_float3(0);
    }

    if(normal.x == -2.0f)
        out.el() = make_uchar4(0,0,0,255);
    else {
        const float3 proj = texproj * pos3D;
        const float2 projPixel = make_float2( proj.x / proj.z + 0.5f, proj.y / proj.z + 0.5f);
        
        const float3 diff = normalize(light - pos3D);
        const float dir = fmaxf(dot(normal, diff), 0.f); // * 255;
        if(projPixel.x < 0 || projPixel.x > texture.size.x-1 || projPixel.y < 0 || projPixel.y > texture.size.y-1 ){
            out.el() = make_uchar4(dir*255,dir*255,dir*255,255);
        } else {
            const uchar3 texcol = texture[make_uint2(projPixel.x, projPixel.y)];
            out.el() = make_uchar4(texcol.x*dir, texcol.y*dir, texcol.z*dir, 255);
        }
    }
}

void renderOculusCam(Image<uchar4> out, const Volume volume, const Image<uchar3> & texture, const Matrix4 view, const Matrix4 texproj, const float nearPlane, const float farPlane, const float step, const float largestep, const float3 light, const float3 ambient){
    dim3 block(32,16);
    OculusCam<<<divup(out.size, block), block>>>(out, volume, texture, view, texproj, nearPlane, farPlane, step, largestep, light, ambient);
}