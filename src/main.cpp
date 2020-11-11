/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "framebuffer.h"
#include "polymesh.h"
#include "kdtree.h"
#include "vector.h"
#include "ray.h"
#include "hit.h"
#include "material.h"
#include "photon_mapping.h"
#include "sphere.cpp"
#include <random>
#include <time.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <functional>
#include <numeric>
#include <algorithm>
#include <queue>
#include <float.h>
#include <vector>
#include <array>
using namespace std;

uniform_real_distribution<float> unif(0.0, 1.0);
mt19937_64 rng;

// Checks if a primary right hits surface
bool surfaceIntersect(Ray ray, Hit *h, vector<Object*> &obj) {
    float t = 1000000000.0;
    bool flag = false;
    for (auto &object : obj) {
        object->intersection(ray, *h);
        if (h->flag) {
            if (h->t < t) {
                t = h->t;
                flag = true;
            }
        }
    }
    return flag;
}

// Uniform distribution of random directions for diffuse point light
Vector randomDir() {
    double x,y,z, temp;
    Vector direction;
    do {
        x = ((double)rand() / (RAND_MAX))*2 - 1;
        y = ((double)rand() / (RAND_MAX))*2 - 1;
        z = ((double)rand() / (RAND_MAX))*2 - 1;
    } while (x*x + y*y + z*z > 1.0);
    direction.x = x;
    direction.y = y;
    direction.z = z;
    return direction;
}

// Estimates the radiance for a given point by locating surrounding neighbours in the kd-tree
Colour radianceCalc(Photon_mapping &pmap, KDTreeNode *root, KDTree kdTree, Hit *h) {

    Colour radiance;
    radiance.r = radiance.g = radiance.b = 0.0;

    // set search radius for nearest neighbours depending on if photon is caustic or
    // normal/shadow
    float radius;
    float estimates;
    // If caustic photon, minimise the search radius for sharp caustics
    if (pmap.causticPhoton) {
        estimates = 2000;
        radius = 0.09;
    }
    else {
        // If normal photon, use larger radius for smooth rendering to avoid too much noise
        estimates = 6950;
//        estimates = 1350;
//        radius = 0.14;
        radius = 0.11;
    }

    // nodeHeap stores the neighbours by distance to the point within the search radius
    priority_queue <KDTreeNode> nodeHeap;
    // Find nearest neighbours with respect to hit point
    array<float, 3> point = {h->position.x, h->position.y, h->position.z};
    // Search for the nearest neighbours from a distance to the given point
    kdTree.locatePhotons(root, point, radius, nodeHeap, estimates);
    // Holds node at top of heap
    KDTreeNode neighbour;
    // Ignore heaps with minimal photon neighbours
    if (nodeHeap.size() < 10) {
        return radiance;
    }

    float maxDistance = nodeHeap.top().distanceToPoint;


    // Keep popping the nearest neighbours from the heap (top of heap is farthest
    // distance from point) and sum their respective power values.
    // The radiance estimate values are scaled by normalized cone filter.
    while (nodeHeap.size()) {

        neighbour = nodeHeap.top();

        radiance.r += neighbour.ptr.power[0];
        radiance.g += neighbour.ptr.power[1];
        radiance.b += neighbour.ptr.power[2];

        nodeHeap.pop();
    }
    //Cone filter
    float k = 1;
    radiance.r *= M_PI*maxDistance*(1-2/(k*3));
    radiance.g *= M_PI*maxDistance*(1-2/(k*3));
    radiance.b *= M_PI*maxDistance*(1-2/(k*3));

    return radiance;
}

// If the hit point is on a dielectric material, calculate the new direction
// or Total Internal Reflection (TIR)
Ray refract(Ray r, Vertex pos, Vector norm, Hit *h) {

    float eta, n1, n2, cosi, n;
    Ray refract;

    n1 = 1.56;
    n2 = 1.0;
    // Snell's Law to determine refracted direction
    if (r.direction.dot(h->normal) > 0) {
        h->normal.negate();
        eta = n1/n2;
    } else {
        eta = n2/n1;
    }

    cosi = -h->normal.dot(r.direction);

    float k = 1 - eta * eta * (1.0 - cosi * cosi);

    Vector d;
    Vertex p;
    // Calculate refractive direction
    d.x = r.direction.x * eta + h->normal.x * (eta * cosi - sqrt(k));
    d.y = r.direction.y * eta + h->normal.y * (eta * cosi - sqrt(k));
    d.z = r.direction.z * eta + h->normal.z * (eta * cosi - sqrt(k));
    d.normalise();
    //shift point from its origin
    p.x = h->position.x - h->normal.x * (1000*FLT_EPSILON);
    p.y = h->position.y - h->normal.y * (1000*FLT_EPSILON);
    p.z = h->position.z - h->normal.z * (1000*FLT_EPSILON);

    refract.position = p;
    refract.direction = d;

    return refract;
}

Colour raytrace(Photon_mapping &pmap, Ray &ray, vector<Object*> &obj, KDTreeNode *root, KDTree kdTree, int bounce, Vertex lPos, Hit *h) {

    Colour radiance, colour;
    colour.r = colour.g = colour.b = 0.0;

    if (++bounce > 25) // Restrict the recursive calls
        return colour;

    // Check if primary ray hits surface
    if (surfaceIntersect(ray, h, obj)) {

        if (h->what->material.diffuseSurface) {

            // Calculate reflected radiance at the surface intersect point
            radiance = radianceCalc(pmap, root, kdTree, h);

            colour.r += radiance.r;
            colour.g += radiance.g;
            colour.b += radiance.b;

        }
        else if (h->what->material.transparent) {

            Ray refractRay;
            refractRay = refract(ray, ray.position, h->normal, h);

            // Recursive call to calculate radiance on diffuse surface
            radiance = raytrace(pmap, refractRay, obj, root, kdTree, bounce + 1, lPos, h);
            // Update the overall "colour" with the returned radiance values
            colour.r += radiance.r;
            colour.g += radiance.g;
            colour.b += radiance.b;
        }
        else if (h->what->material.reflect) {

            Ray reflectionRay;
            float incident = ray.direction.dot(h->normal);

            Vector reflectionDir;
            Vertex reflectionPos;
            // Calculate reflective direction
            reflectionDir.x = ray.direction.x - 2.0 * incident * h->normal.x;
            reflectionDir.y = ray.direction.y - 2.0 * incident * h->normal.y;
            reflectionDir.z = ray.direction.z - 2.0 * incident * h->normal.z;
            reflectionDir.normalise();
            // shift point from its origin
            reflectionPos.x = h->position.x + (1000 * FLT_EPSILON) * h->normal.x;
            reflectionPos.y = h->position.y + (1000 * FLT_EPSILON) * h->normal.y;
            reflectionPos.z = h->position.z + (1000 * FLT_EPSILON) * h->normal.z;

            reflectionRay.position = reflectionPos;
            reflectionRay.direction = reflectionDir;
            // Recursive call to calculate radiance on diffuse surface
            radiance = raytrace(pmap, reflectionRay, obj, root, kdTree, bounce + 1, lPos, h);
            // Update the overall "colour" with the returned radiance values
            colour.r += radiance.r;
            colour.g += radiance.g;
            colour.b += radiance.b;
        }
    }
    return colour;
}


int main(int argc, char *argv[])
{
    srand(time(NULL));

    // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
    Transform *transform = new Transform(0.25f, 0.0f, 0.0f,  0.0f,   //Transform teapot object
                                         0.0f, 0.0f, 0.25f, -0.84f,
                                         0.0f, 0.25f, 0.0f, -4.17f,
                                         0.0f, 0.0f, 0.0f, 1.0f);

    Transform *transform2 = new Transform(2.0f, 0.0f, 0.0f,  0.0f,  //Transform spheres and box object
                                          0.0f, 2.0f, 0.0f, 0.0f,
                                          0.0f, 0.0f, 2.0f, -4.0f,
                                          0.0f, 0.0f, 0.0f, 1.0f);

    Transform *transform3 = new Transform(0.3f, 0.0f, 0.0f,  0.0f,  //Transform spheres and box object
                                          0.0f, 0.3f, 0.0f, 3.0f,
                                          0.0f, 0.0f, 0.03f, -3.5f,
                                          0.0f, 0.0f, 0.0f, 1.0f);

    // Read in the teapot model.
    PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);
    //Set box
    PolyMesh *right = new PolyMesh((char *)"right.ply", transform2);
    PolyMesh *back = new PolyMesh((char *)"back.ply", transform2);
    PolyMesh *left = new PolyMesh((char *)"left.ply", transform2);
    PolyMesh *light = new PolyMesh((char *)"top.ply", transform3);
    PolyMesh *ceiling = new PolyMesh((char *)"top.ply", transform2);
    PolyMesh *floor1 = new PolyMesh((char *)"floor.ply", transform2);

    // create objects for scene
    Vertex *sp = new Vertex(-0.5, -0.6, -3.58);
    Vertex *sp2 = new Vertex(0.53, 0.4, -3.88);

    //create sphere objects
    Sphere *sphere = new Sphere(*sp, 0.3);
    Sphere *sphere2 = new Sphere(*sp2, 0.33);

    //*******************************************************************//
    // Set thedifferent materials

    Material bp1;
    Material bp2;
    Material bp3;
    Material bp4;
    Material bp5;
    Material bp6;
    Material bp7;
    Material floor;
    Material bp8;
    Material areaLight;

    bp1.reflect = true;
    bp1.transparent = false;
    bp1.diffuseSurface = false;
    bp1.ai = 1.0;
    bp1.ambient.r = 0.0f;
    bp1.ambient.g = 0.0f;
    bp1.ambient.b = 0.0f;
    bp1.diffuse.r = 0.0f;
    bp1.diffuse.g = 0.0f;
    bp1.diffuse.b = 0.0f;
    bp1.specular.r = 1.0f;
    bp1.specular.g = 1.0f;
    bp1.specular.b = 1.0f;
    bp1.power = 150.0f;
    bp1.ri = 1.56;

    bp2.reflect = false;
    bp2.transparent = true;
    bp2.diffuseSurface = false;
    bp2.ai = 1.0;
    bp2.ambient.r = 0.0f;
    bp2.ambient.g = 0.0f;
    bp2.ambient.b = 0.0f;
    bp2.diffuse.r = 0.5f;
    bp2.diffuse.g = 0.5f;
    bp2.diffuse.b = 0.5f;
    bp2.specular.r = 1.0f;
    bp2.specular.g = 1.0f;
    bp2.specular.b = 1.0f;
    bp2.power = 125.0f;
    bp2.ri = 1.56;

    bp3.reflect = false;
    bp3.transparent = false;
    bp3.diffuseSurface = true;
    bp3.ai = 0.9;
    bp3.ambient.r = 0.85;
    bp3.ambient.g = 0.2f;
    bp3.ambient.b = 0.8f;
    bp3.diffuse.r = 1.0f;
    bp3.diffuse.g = 1.0f;
    bp3.diffuse.b = 1.0f;
    bp3.specular.r = 0.1f;
    bp3.specular.g = 0.1f;
    bp3.specular.b = 0.1f;
    bp3.power = 10.0f;
    bp3.ri = 1.0;

    bp4.reflect = false;
    bp4.transparent = false;
    bp4.diffuseSurface = true;
    bp4.ai = 0.75;
    bp4.ambient.r = 0.1f;
    bp4.ambient.g = 0.9f;
    bp4.ambient.b = 1.0f;
    bp4.diffuse.r = 0.1f;
    bp4.diffuse.g = 0.9f;
    bp4.diffuse.b = 1.0f;
    bp4.specular.r = 0.0f;
    bp4.specular.g = 0.0f;
    bp4.specular.b = 0.0f;
    bp4.power = 300.0f;
    bp4.ri = 1.53;

    bp5.reflect = false;
    bp5.transparent = false;
    bp5.diffuseSurface = true;
    bp5.ai = 0.8;
    bp5.ambient.r = 0.9f;
    bp5.ambient.g = 0.9f;
    bp5.ambient.b = 0.9f;
    bp5.diffuse.r = 1.0f;
    bp5.diffuse.g = 1.0f;
    bp5.diffuse.b = 1.0f;
    bp5.specular.r = 0.0f;
    bp5.specular.g = 0.0f;
    bp5.specular.b = 0.0f;
    bp5.power = 0.0f;
    bp5.ri = 1.0;

    bp6.reflect = false;
    bp6.transparent = false;
    bp6.diffuseSurface = false;
    bp6.ai = 0.75;
    bp6.ambient.r = 0.04f;
    bp6.ambient.g = 0.12f;
    bp6.ambient.b = 0.73f;
    bp6.diffuse.r = 0.4f;
    bp6.diffuse.g = 0.12f;
    bp6.diffuse.b = 0.73f;
    bp6.specular.r = 0.f;
    bp6.specular.g = 0.f;
    bp6.specular.b = 0.f;
    bp6.power = 0.0f;
    bp6.ri = 1.53;

    bp7.reflect = false;
    bp7.transparent = false;
    bp7.diffuseSurface = true;
    bp7.ai = 0.75;
    bp7.ambient.r = 1.0f;
    bp7.ambient.g = 0.12f;
    bp7.ambient.b = 0.73f;
    bp7.diffuse.r = 1.0f;
    bp7.diffuse.g = 0.08f;
    bp7.diffuse.b = 0.08f;
    bp7.specular.r = 0.0f;
    bp7.specular.g = 0.0f;
    bp7.specular.b = 0.0f;
    bp7.power = 300.0f;
    bp7.ri = 1.53;

    floor.reflect = false;
    floor.transparent = false;
    floor.diffuseSurface = true;
    floor.ai = 0.75;
    floor.ambient.r = 0.9f;
    floor.ambient.g = 0.9f;
    floor.ambient.b = 0.9f;
    floor.diffuse.r = 0.9f;
    floor.diffuse.g = 0.9f;
    floor.diffuse.b = 0.9f;
    floor.specular.r = 0.3f;
    floor.specular.g = 0.3f;
    floor.specular.b = 0.3f;
    floor.power = 300.0f;
    floor.ri = 1.53;

    bp8.light = true;
    bp8.reflect = false;
    bp8.transparent = false;
    bp8.diffuseSurface = false;
    bp8.ai = 0.0;
    bp8.ambient.r = 0.1f;
    bp8.ambient.g = 0.1f;
    bp8.ambient.b = 0.1f;
    bp8.diffuse.r = 0.1f;
    bp8.diffuse.g = 0.1f;
    bp8.diffuse.b = 0.1f;
    bp8.specular.r = 0.5f;
    bp8.specular.g = 0.5f;
    bp8.specular.b = 0.5f;
    bp8.power = 30.0f;
    bp8.ri = 1.53;


    sphere->material = bp2;
    sphere2->material = bp1;

    pm->material = bp5;

    right->material = bp4;
    back->material = bp5;
    ceiling->material = bp8;
    left->material = bp7;
    floor1->material = bp5;
    light->material = bp8;

    std::vector <Object*> obj;

    // Box object
    obj.push_back(back);
    obj.push_back(ceiling);
    obj.push_back(right);
    obj.push_back(left);
    obj.push_back(floor1);
    obj.push_back(light);


    // Teampot object
    obj.push_back(pm);
//
    //spheres
    obj.push_back(sphere);

    obj.push_back(sphere2);


//
    //*******************************************************************//
    // Set diffuse point light position
    Vertex *lpos = new Vertex(0.0,0.35,-3.24);

    //*******************************************************************//
    // Initalize the photon map
    Photon_mapping pmap(constants::GLOBAL_CONST_PHOTONS);

    cout << "Emitting photons..." << endl;
    int emittedPhotons = 0;
    pmap.causticPhoton = false;
    //Emits photons into the scene
    Hit *phit = new Hit();
    while (emittedPhotons < constants::GLOBAL_CONST_PHOTONS) {
        Vector direction = randomDir(); //Generate random direction to emit photon in
        direction.normalise();
        Photon p;
        p.power[0] = p.power[1] = p.power[2] = 1.0;

        pmap.photon_trace(direction, *lpos, obj, 0, p, phit);

        emittedPhotons++;
    }

    cout << "Emitted photons: " << emittedPhotons << endl;
    cout << "Photon bounces stored: " << pmap.stored_photons << endl;

    //*******************************************************************//
    // A balanced kd tree is built to store the photon map, where each parent
    // node is the median of its child nodes for a given x,y,z plane, where the
    // plane is equal to (current depth) % 3, x=0, y=1, z=3

    KDTree kdTree;
    KDTreeNode *root = NULL;

    cout << endl;
    vector<vector<float>> photonPosition;
    vector<vector<float>> idx; // Holds all xyz, yzx and zxy indices
    vector<float> tempx, tempy, tempz; // // Holds their respective x,y,z positions of all photons

    Photon temp;
    // Loop through all stored photons in the map
    for (int i = 0; i < pmap.stored_photons; i++) {
        temp = pmap.photons[i];
        photonPosition.push_back({temp.x,temp.y,temp.z}); // Store the x,y,z position of a given photon

        tempx.push_back(temp.x); //Store x position of all photons
        tempy.push_back(temp.y); //Store y position of all photons
        tempz.push_back(temp.z); //Store z position of all photons
    }

    idx.push_back(tempx);
    idx.push_back(tempy);
    idx.push_back(tempz);

    vector<size_t> indexX, indexY, indexZ;
    //Sorts the indices of photons from smallest value to largest
    // value on their given plane.
    indexX = kdTree.sortIndex(idx, idx[0].size(), 0); // Sorts photons with respect to XYZ values (plane x)
    indexY = kdTree.sortIndex(idx, idx[0].size(), 1); // Sorts photons with respect to YZX values (plane y)
    indexZ = kdTree.sortIndex(idx, idx[0].size(), 2); // Sorts photons with respect to ZXY values (plane z)

    cout << "Building kd tree..." << endl;
    root = kdTree.buildTree(root, indexX, indexY, indexZ, idx, 0, pmap);

    // Power of all stored photons is scaled by the number of emitted photons
    // (scaling before sorting produces unwanted results)
    pmap.scalePower(emittedPhotons);


    //*******************************************************************//
    // Repeat for caustic photons
    cout << "Emitting photons..." << endl;
    int causticEmitted = 0;
    pmap.causticPhoton = true;
    //Emits photons into the scene
    while (causticEmitted < constants::GLOBAL_CONST_CAUSTIC_PHOTONS) {
        Vector direction = randomDir(); //Generate random direction to emit photon in
        direction.normalise();
        Photon p;
        pmap.specularSurfaceHit = 0;
        p.power[0] = p.power[1] = p.power[2] = 1.0;

        pmap.photon_trace(direction, *lpos, obj, 0, p, phit);

        causticEmitted++;
    }

    cout << "Caustic photons emitted: " << causticEmitted << endl;
    cout << "Caustic photon bounces stored: " << pmap.causticStoredPhotons << endl;

    //*******************************************************************//
    // A balanced kd tree is built to store the photon map, where each parent
    // node is the median of its child nodes for a given x,y,z plane, where the
    // plane is equal to (current depth) % 3, x=0, y=1, z=3

    KDTree kdTreeCaustic;
    KDTreeNode *rootCaustic = NULL;

    cout << endl;
    vector<vector<float>> photonPositionCaustic;
    vector<vector<float>> idxCaustic; // Holds all xyz, yzx and zxy indices
    vector<float> tempxCaustic, tempyCaustic, tempzCaustic; // // Holds their respective x,y,z positions of all photons

    Photon tempCaustic;
    // Loop through all stored photons in the map
    for (int i = 0; i < pmap.causticStoredPhotons; i++) {
        tempCaustic = pmap.causticPhotons[i];
        photonPositionCaustic.push_back({tempCaustic.x,tempCaustic.y,tempCaustic.z}); // Store the x,y,z position of a given photon

        tempxCaustic.push_back(tempCaustic.x); //Store x position of all photons
        tempyCaustic.push_back(tempCaustic.y); //Store y position of all photons
        tempzCaustic.push_back(tempCaustic.z); //Store z position of all photons
    }

    idxCaustic.push_back(tempxCaustic);
    idxCaustic.push_back(tempyCaustic);
    idxCaustic.push_back(tempzCaustic);

    vector<size_t> indexXCaustic, indexYCaustic, indexZCaustic;
    //Sorts the indices of photons from smallest value to largest
    // value on their given plane.
    indexXCaustic = kdTreeCaustic.sortIndex(idxCaustic, idxCaustic[0].size(), 0); // Sorts photons with respect to XYZ values (plane x)
    indexYCaustic = kdTreeCaustic.sortIndex(idxCaustic, idxCaustic[0].size(), 1); // Sorts photons with respect to YZX values (plane y)
    indexZCaustic = kdTreeCaustic.sortIndex(idxCaustic, idxCaustic[0].size(), 2); // Sorts photons with respect to ZXY values (plane z)

    cout << "Building caustic kd tree..." << endl;
    rootCaustic = kdTreeCaustic.buildTree(rootCaustic, indexXCaustic, indexYCaustic, indexZCaustic, idxCaustic, 0, pmap);

    // Power of all stored photons is scaled by the number of emitted photons
    // (scaling before sorting produces unwanted results)
    pmap.scalePower(causticEmitted);


    //*******************************************************************//
    //Image dimensions
    int width = 650;
    int height = 650;

    // Create a framebuffer
    FrameBuffer *fb = new FrameBuffer(width, height);

    float aspect_ratio = width/height;
    float xn, yn, sx, sy;   //Normalized pixels and Screen space
    float depth = 100000000.0;

    Hit *hit = new Hit();
    Vertex *camera = new Vertex(0.0, 0.0, 0.0);
    Vector *dir = new  Vector(0, 0, -1);
    Ray *ray = new Ray(*camera, *dir);

    float fov = 40.0;
    float angle = tan(0.5*fov*M_PI/180.0);

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            Colour radiance, outputColour;
            outputColour.r = outputColour.g = outputColour.b = 0.0;

            for (int x = 0; x < 2; x++ ) {
                for (int y = 0; y < 2; y++) {
                    xn = (i + unif(rng)) / width;    //Set to 1 of 4 quadrants of the pixel
                    yn = (j + unif(rng)) / height;

                    sx = (2 * xn - 1) * aspect_ratio * angle; //Adjust field of view
                    sy = (1 - 2 * yn) * angle;

                    dir->x = sx - camera->x;    //set direction and position of ray from viewpoint
                    dir->y = sy - camera->y;
                    dir->z = -1;
                    dir->normalise();
                    ray->position = *camera;
                    ray->direction = *dir;

                    radiance.r = radiance.g = radiance.b = 0.0;   //Reset radiance values for each sample

                    pmap.causticPhoton = false;
                    //normal and shadow photons
                    radiance = raytrace(pmap, *ray, obj, root, kdTree, 0, *lpos, hit);
                    // Sum sample radiance values for a given pixel
                    outputColour.r += radiance.r;
                    outputColour.g += radiance.g;
                    outputColour.b += radiance.b;

                    radiance.r = radiance.g = radiance.b = 0.0;   //Reset radiance values for each sample

                    dir->x = sx - camera->x;    //set direction and position of ray from viewpoint
                    dir->y = sy - camera->y;
                    dir->z = -1;
                    dir->normalise();
                    ray->position = *camera;
                    ray->direction = *dir;

                    pmap.causticPhoton = true;
                    //normal and shadow photons
                    radiance = raytrace(pmap, *ray, obj, rootCaustic, kdTreeCaustic, 0, *lpos, hit);
                    // Sum sample radiance values for a given pixel
                    outputColour.r += radiance.r;
                    outputColour.g += radiance.g;
                    outputColour.b += radiance.b;
                }
            }

            //divide by number of samples to reconstruct original pixel(box filter)
            outputColour.r /= 4;
            outputColour.g /= 4;
            outputColour.b /= 4;

            //gamma correction? (1^2.2)
            outputColour.r = pow(outputColour.r, 1/2.2);
            outputColour.g = pow(outputColour.g, 1/2.2);
            outputColour.b = pow(outputColour.b, 1/2.2);
            fb->plotPixel(i, j, outputColour.r, outputColour.g, outputColour.b);
        }
        cerr << "*" << flush;
    }

//   Output the framebuffer.
    fb->writeRGBFile((char *)"rgb_output.ppm");

    return 0;
}
