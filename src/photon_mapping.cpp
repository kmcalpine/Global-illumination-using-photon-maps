/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
comment
 */
#include <iostream>
#include <fstream>
#include <random>
#include <float.h>
#include <vector>
#include <numeric>
#include <array>

#include "framebuffer.h"
#include "polymesh.h"
#include "vector.h"
#include "ray.h"
#include "hit.h"
#include "photon_mapping.h"

using namespace std;
// Initialize a photon map
Photon_mapping::Photon_mapping(int n_photons) {
    stored_photons = 0;
    causticStoredPhotons = 0;
    max_photons = n_photons;
}
//scale photon powers uniformly by total number of photons
void Photon_mapping::scalePower(int emittedPhotons) {
    for (int i = 0; i < stored_photons; ++i) {
        photons[i].power[0] *= 1 / emittedPhotons;
        photons[i].power[1] *= 1 / emittedPhotons;
        photons[i].power[2] *= 1 / emittedPhotons;
    }
}

void Photon_mapping::randomDir(Vector &direction) {
    double x,y,z;
    do {
        x = ((double)rand() / (RAND_MAX))*2 - 1;
        y = ((double)rand() / (RAND_MAX))*2 - 1;
        z = ((double)rand() / (RAND_MAX))*2 - 1;
    } while (x*x + y*y + z*z > 1.0);
    direction.x = x;
    direction.y = y;
    direction.z = z;
}

// Traces the reflected photon after hitton adiffuse surface
void Photon_mapping::diffuseReflection(Vector &direction, Hit *hit, vector<Object *> &obj, int bounce, Photon &p,
                                       float p_d){

    array<float, 3> powerRefl = {p.power[0]*hit->what->material.diffuse.r*p_d,  //Power attenuation
                                 p.power[1]*hit->what->material.diffuse.g*p_d,
                                 p.power[2]*hit->what->material.diffuse.b*p_d};

    p.power = powerRefl;

    Photon newPhoton = p;

    Vertex new_pos = hit->position;

    newPhoton.power[0] = powerRefl[0];
    newPhoton.power[1] = powerRefl[1];
    newPhoton.power[2] = powerRefl[2];
    // Slightly shift from hit position
    new_pos.x += hit->normal.x*(1000*FLT_EPSILON);
    new_pos.y += hit->normal.y*(1000*FLT_EPSILON);
    new_pos.z += hit->normal.z*(1000*FLT_EPSILON);

    Vector reflectionDir;

    float incident = direction.dot(hit->normal);

    reflectionDir.x = direction.x - 2.0 * incident * hit->normal.x;
    reflectionDir.y = direction.y - 2.0 * incident * hit->normal.y;
    reflectionDir.z = direction.z - 2.0 * incident * hit->normal.z;
    direction = reflectionDir;
    direction.normalise();

    // If photon has previously hit/passed through specular surface, add to caustic map
    if (specularSurfaceHit >= 1 && causticPhoton) {
        Photon newPhoton = p;
        causticPhotons.push_back(newPhoton);
        causticStoredPhotons++;
    }
    else {
        Photon newPhoton = p;
        photons.push_back(newPhoton);
        stored_photons++;
    }
    photon_trace(direction, new_pos, obj, bounce, newPhoton, hit);
}

// Checks if an emitted photon from a light source hits a surface, then determines any shadow rays
// that should be stored from additional surface interactions beyond the initial hit point
bool Photon_mapping::photonHit(Ray pmap_ray, vector<Object *> &obj, Hit *hit, Photon &p) {
    bool hitFlag = false;
    Hit *shadowHit;
    float shadowDistance = 0.0;
    Photon shadowPhoton;
    int shadowCount = 0;
    float t = 1000000000.0;

    for (auto &object : obj) {
        object->intersection(pmap_ray, *hit);
        if (hit->flag) {
            hitFlag = true;
//            if ((hit->t > shadowDistance) || (hit->normal.dot(shadowHit->normal)) > 0.0) {
            if ((hit->t > shadowDistance)) {
                Hit *tempHit = new Hit(*hit);
                shadowHit = tempHit;
                shadowDistance = shadowHit->t;
                shadowPhoton = p;
                shadowCount++;
                delete tempHit;
            }
            if (hit->t < t) {
                t = hit->t;
            }
        }
    }
    // ** Fix this ** //
    if (shadowCount > 1 && !causticPhoton) {
        shadowPhoton.power[0] = -0.01 * shadowHit->what->material.diffuse.r;
        shadowPhoton.power[1] = -0.01 * shadowHit->what->material.diffuse.g;
        shadowPhoton.power[2] = -0.01 * shadowHit->what->material.diffuse.b;
        photons.push_back(shadowPhoton);
    }
    return hitFlag;
}

void Photon_mapping::photon_trace(Vector direction, Vertex lpos, vector<Object *> &obj, int bounce, Photon &p, Hit *hit) {


    // If bounce limit reached, return
    if (++bounce > constants::GLOBAL_CONST_BOUNCE) {
        return;
    }

    Ray pmap_ray;
    pmap_ray.direction = direction;
    pmap_ray.position = lpos;
    // Check if photon hits a surface, if true, use russian roulette to determine
    // if the photon is stored or reflected
    if(!photonHit(pmap_ray, obj, hit, p)) {
        return;
    }

    Vertex new_pos = hit->position;

    p.x = hit->position.x;
    p.y = hit->position.y;
    p.z = hit->position.z;

    //*******************************************************************************//
    // Create probabilities of reflection/absorb
    //max reflective value of surface
    float p_maxRef = max(hit->what->material.diffuse.r + hit->what->material.specular.r,
                    max(hit->what->material.diffuse.g + hit->what->material.specular.g,
                        hit->what->material.diffuse.b + hit->what->material.specular.b));
    //probability of absorption
    float p_absorb = 1 - p_maxRef;
    //probability of diffuse reflection
    float p_diffuse = ((hit->what->material.diffuse.r + hit->what->material.diffuse.g + hit->what->material.diffuse.b)/(
            hit->what->material.diffuse.r + hit->what->material.diffuse.g + hit->what->material.diffuse.b +
            hit->what->material.specular.r + hit->what->material.specular.g + hit->what->material.specular.b))*p_maxRef;
    //probability of specular reflection
    float p_specular = p_maxRef - p_diffuse;

    //*******************************************************************************//
    //random generator for russian roulette probability
    uniform_real_distribution<float> distribution(0,1);
    float prob = distribution(generator);

    //diffuse reflection
    if (prob >= 0 && prob < p_diffuse) {
        diffuseReflection(direction, hit, obj, bounce, p, p_diffuse);
    }
    //specular reflection
    else if (prob >= p_diffuse && prob < p_diffuse + p_specular) {

        //increment specular surface hit for caustic photon tracking
        if (causticPhoton && (hit->what->material.transparent || hit->what->material.reflect)) specularSurfaceHit++;

        // Copy variables to return to previous point in order to compute both
        // refraction and specular reflection
        Vector dir = direction;
        Vertex pos = new_pos;
        Hit *specHit = new Hit(*hit);
        Photon temp = p;

        if (hit->what->material.transparent) {
            refraction(direction, new_pos, obj, bounce, p, hit);
        }
        specularReflection(direction, new_pos, obj, bounce, p, hit, p_specular);
    }
    //absorption
    else if (prob >= p_specular + p_diffuse && prob < 1) {
//        array<float, 3> powerRefl = {p.power[0]*hit->what->material.ambient.r*p_a, //Power attenuation
//                                     p.power[1]*hit->what->material.ambient.g*p_a,
//                                     p.power[2]*hit->what->material.ambient.b*p_a};
//        p.power = powerRefl;
//
//        Photon newPhoton = p;
//
//        newPhoton.power[0] = powerRefl[0];
//        newPhoton.power[1] = powerRefl[1];
//        newPhoton.power[2] = powerRefl[2];


//      store absorbed photon
        photons.push_back(p);
        stored_photons++;
        return;
    }
}

void Photon_mapping::refraction(Vector &direction, Vertex &lpos, vector<Object *> &obj, int bounce, Photon &p, Hit *hit) {

    float eta, n1, n2, cosi;

    n1 = 1.56; //refractive index
    n2 = 1.0;
    //inside object (TIR)
    if (direction.dot(hit->normal) > 0) {
        hit->normal.negate();
        eta = n1/n2; //ratio of indices of refraction
    } else {
        eta = n2/n1; //ratio of indices of refraction
    }

    cosi = -hit->normal.dot(direction);

    float k = 1.0 - eta * eta * (1.0 - cosi * cosi);

    Vector dir;

    dir.x = direction.x * eta + hit->normal.x * (eta * cosi - sqrt(k));
    dir.y = direction.y * eta + hit->normal.y * (eta * cosi - sqrt(k));
    dir.z = direction.z * eta + hit->normal.z * (eta * cosi - sqrt(k));
    dir.normalise();

    lpos.x = hit->position.x - hit->normal.x * (1000*FLT_EPSILON);
    lpos.y = hit->position.y - hit->normal.y * (1000*FLT_EPSILON);
    lpos.z = hit->position.z - hit->normal.z * (1000*FLT_EPSILON);
    // trace refracted photon
    photon_trace(dir, lpos, obj, bounce, p, hit);
}

void Photon_mapping::specularReflection(Vector &direction, Vertex &new_pos, vector<Object *> &obj, int bounce, Photon &p, Hit *hit, float p_s) {

    Photon newPhoton = p;
    //adjust reflected power
    array<float, 3> powerRefl = {p.power[0]*hit->what->material.diffuse.r*p_s,  //Power attenuation
                                 p.power[1]*hit->what->material.diffuse.g*p_s,
                                 p.power[2]*hit->what->material.diffuse.b*p_s};

    newPhoton.power = powerRefl;

    new_pos.x += hit->normal.x * (1000 * FLT_EPSILON);
    new_pos.y += hit->normal.y * (1000 * FLT_EPSILON);
    new_pos.z += hit->normal.z * (1000 * FLT_EPSILON);

    Vector reflectionDir;

    float incident = direction.dot(hit->normal);

    //calculate reflective direction
    reflectionDir.x = direction.x - 2.0 * incident * hit->normal.x;
    reflectionDir.y = direction.y - 2.0 * incident * hit->normal.y;
    reflectionDir.z = direction.z - 2.0 * incident * hit->normal.z;
    direction = reflectionDir;
    direction.normalise();
    //trace reflected specular photon
    photon_trace(direction, new_pos, obj, bounce, newPhoton, hit);
}