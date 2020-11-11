/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
comment
 */
#include "vector.h"
#include "photon_mapping.h"
#include "kdtree.h"
#include <array>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <functional>
#include <queue>

using namespace std;

// A method to create a node of K D tree
KDTreeNode* KDTree::newNode(int index, Photon_mapping &pmap, int cd){
    KDTreeNode *temp = new KDTreeNode();
    if (pmap.causticPhoton) {
        temp->ptr = pmap.causticPhotons[index];
    } else {
        temp->ptr = pmap.photons[index];
    }

    temp->ptr.flag = cd; //set flag for splitting plane equal to the current depth (0,1 or 2)
    temp->idx = index;

    temp->left = NULL;
    temp->right = NULL;

    return temp;
}

KDTreeNode* KDTree::insertNode(KDTreeNode *root, array<float, 3> point, int depth, int index, Photon_mapping &pmap) {

    // Calculate current dimension (cd) of comparison
    int cd = depth % 3;

    // Tree is empty?
    if (root == NULL) {
        nodes++;
        return newNode(index, pmap, cd);
    }
    // Compare the new point with root on current dimension 'cd'
    // and decide the left or right subtree
    array<float, 3> rootPos;
    rootPos[0] = root->ptr.x;
    rootPos[1] = root->ptr.y;
    rootPos[2] = root->ptr.z;
    if (rootPos==point) {
        return root;
    }
    if (point[cd] < (rootPos[cd])) {
        root->left = insertNode(root->left, point, depth + 1, index, pmap);
    }
    else {
        root->right = insertNode(root->right, point, depth + 1, index, pmap);
    }
    return root;
}


KDTreeNode* KDTree::buildTree(KDTreeNode *root, vector<size_t> idxX, vector<size_t> idxY, vector<size_t> idxZ, vector<vector<float>> indices, int axis, Photon_mapping &pmap) {
    vector<size_t> TempLowx;
    vector<size_t> TempHighx;

    vector<size_t> TempLowy;
    vector<size_t> TempHighy;

    vector<size_t> TempLowz;
    vector<size_t> TempHighz;
    array<float, 3> r;
    array<float, 3> rootPos;

    //use median of index array
    int d = axis%3;
    int index;

    // Determine current splitting plane axis,
    // in order to sort by XYZ,YXZ or ZXY
    if (idxX.size() <= 0) return root;
    if (d == 0) index = idxX[idxX.size()/2];
    if (d == 1) index = idxY[idxY.size()/2];
    if (d == 2) index = idxZ[idxZ.size()/2];

    //root node
    rootPos[0] = indices[0][index];
    rootPos[1] = indices[1][index];
    rootPos[2] = indices[2][index];


    root = insertNode(root, rootPos, 0, index, pmap);

    //if less than 2, child nodes are already split
    // Recursively split by the median, where the median becomes the parent node
    if (idxX.size() >= 2 ||idxY.size() >= 2 || idxZ.size() >= 2) {

        // Compare value to the super key value to determine if it should go to the left or right
        for (int i = 0; i < idxX.size(); i++) {
            r[0] = indices[0][idxX[i]];
            r[1] = indices[1][idxX[i]];
            r[2] = indices[2][idxX[i]];

            if (rootPos == r) {
            } else if (r[d] < rootPos[d]) {
                TempLowx.push_back(idxX[i]);

            } else {
                TempHighx.push_back(idxX[i]);
            }

            r[0] = indices[0][idxY[i]];
            r[1] = indices[1][idxY[i]];
            r[2] = indices[2][idxY[i]];

            if (rootPos == r) {
            } else if (r[d] < rootPos[d]) {
                TempLowy.push_back(idxY[i]);
            } else {
                TempHighy.push_back(idxY[i]);
            }

            r[0] = indices[0][idxZ[i]];
            r[1] = indices[1][idxZ[i]];
            r[2] = indices[2][idxZ[i]];

            if (rootPos == r) {
            } else if (r[d] < rootPos[d]) {
                TempLowz.push_back(idxZ[i]);
            } else {
                TempHighz.push_back(idxZ[i]);
            }
        }
        root = buildTree(root, TempLowx, TempLowy, TempLowz, indices, d+1, pmap);
        root = buildTree(root, TempHighx, TempHighy, TempHighz, indices, d+1, pmap);
    }
    return root;
}


void KDTree::locatePhotons(KDTreeNode *root, array<float, 3> point, float maxDistance, priority_queue<KDTreeNode> &nodeHeap, int estimates)
{
    float plane;
    if (root != NULL) {
        if (root->ptr.flag == 0) plane = root->ptr.x;
        if (root->ptr.flag == 1) plane = root->ptr.y;
        if (root->ptr.flag == 2) plane = root->ptr.z;

        if (nodeHeap.size() >= estimates) {
            maxDistance = nodeHeap.top().distanceToPoint;
            nodeHeap.pop(); //heap is full, pop to reduce max distance and prune unwanted
        }

        float distance = plane - point[root->ptr.flag]; //check distance from splitting plane
        // distance to splitting plane
        if (distance > 0) {
            // left side of plane
            locatePhotons(root->left, point, maxDistance, nodeHeap, estimates);
            if (distance*distance < maxDistance*maxDistance) {
                locatePhotons(root->right, point, maxDistance, nodeHeap, estimates);
            }
        } else {
            // right of plane
            locatePhotons(root->right, point, maxDistance, nodeHeap, estimates);
            if (distance*distance < maxDistance*maxDistance) {
                locatePhotons(root->left, point, maxDistance, nodeHeap, estimates);
            }
        }

        float dSquared = (root->ptr.x - point[0]) * (root->ptr.x - point[0]) +
                         (root->ptr.y - point[1]) * (root->ptr.y - point[1]) +
                         (root->ptr.z - point[2]) * (root->ptr.z - point[2]);
        root->distanceToPoint = dSquared;

        if (dSquared < maxDistance*maxDistance) {
            nodeHeap.push(*root);
        }
    }
}

// Returns vectors of indices of the sorted range of photon distances (close to far)
// One sorted vector for each splitting plane:
// plane 0 = XYZ, sorts by minimising x first, followed by y, then finally z
// plane 1 = YZX, sorts by minimising y first, followed by z, then finally x
// plane 2 = ZXY, sorts by minimising z first, followed by x, then finally y
vector<size_t> KDTree::sortIndex(vector<vector<float>> v, int u, int axis)
{
    // initialize original index locations
    vector<size_t> idx(u);
    iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v,&axis](size_t i1, size_t i2) {
            if(v[axis][i1] == v[axis][i2]) {
                if(v[(axis+1)%3][i1] == v[(axis+1)%3][i2]) {
                    return v[(axis+2)%3][i1] <= v[(axis+2)%3][i2];
                }  return v[(axis+1)%3][i1] <= v[(axis+1)%3][i2];
            }
        return v[axis][i1] <= v[axis][i2];});
    return idx;
}
