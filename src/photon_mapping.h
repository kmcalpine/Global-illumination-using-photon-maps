/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "hit.h"
#include "sphere.h"

#include <vector>
#include <random>
#include <array>
using namespace std;

namespace constants
{
	const int GLOBAL_CONST_PHOTONS = 300000;
//	const int GLOBAL_CONST_PHOTONS = 20000;
	const int GLOBAL_CONST_CAUSTIC_PHOTONS = 200000;
	const int GLOBAL_CONST_PHOTON_ESTIMATES = 50;
	const int GLOBAL_CONST_BOUNCE = 7;
}

struct Photon {
    float x,y,z; //pos
    array<float, 3> power; //power
    char phi, theta; //direction
    short flag;
};

class Photon_mapping {
public:
	default_random_engine generator;
	int stored_photons;
	int causticStoredPhotons;
	bool causticPhoton;
	int specularSurfaceHit;
	int max_photons;
    vector<Photon> photons;
    vector<Photon> causticPhotons;
	float bounding_box_minx, bounding_box_maxx, bounding_box_miny, bounding_box_maxy, bounding_box_minz, bounding_box_maxz;

	Photon_mapping(int photons);
	bool photonHit(Ray pmap_ray, vector<Object *> &obj, Hit *hit, Photon &p);
	void randomDir(Vector &direction);
	void photon_trace(Vector direction, Vertex lpos, vector <Object*> &obj, int bounce, Photon &p, Hit *hit);
	void diffuseReflection(Vector &direction, Hit *hit, vector<Object *> &obj, int bounce, Photon &p,
						   float p_d);
	void specularReflection(Vector &direction, Vertex &lpos, vector<Object *> &obj, int bounce, Photon &p, Hit *hit, float p_s);
	void refraction(Vector &direction, Vertex &lpos, vector<Object *> &obj, int bounce, Photon &p, Hit *hit);
	void scalePower(int emittedPhotons);
};
