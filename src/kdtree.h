/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include "photon_mapping.h"
#include <queue>
using namespace std;

//A struct of kd tree node
struct KDTreeNode
{
	KDTreeNode *left, *right;
	Photon ptr;
	int idx;
	float distanceToPoint; //distance to a given point, used for creating max heap of neighbours
	bool operator<(const KDTreeNode &a) const {
        return a.distanceToPoint > distanceToPoint; //push new MAX distance from given point to top of heap
	}
};

class KDTree {
public:
	int nodes = 0;
	int treeDepth = 0;

	KDTreeNode* newNode(int index, Photon_mapping &pmap, int cd);
	KDTreeNode* insertNode(KDTreeNode *root, array<float, 3> point, int depth, int index, Photon_mapping &pmap);
	KDTreeNode* buildTree(KDTreeNode *root, vector<size_t> idxX, vector<size_t> idxY, vector<size_t> idxZ, vector<vector<float>> mergetest1, int axis, Photon_mapping &pmap);

	void locatePhotons(KDTreeNode *root, array<float, 3> testpoint2, float maxDistance, priority_queue<KDTreeNode> &nodeHeap, int estimates);
	vector<size_t> sortIndex(vector<vector<float>> v, int u, int axis);
};

#endif /* KDTREE_H */
