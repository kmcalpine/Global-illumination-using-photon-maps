/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "transform.h"
#include "object.h"

typedef int TriangleIndex[3];

class PolyMesh:public Object{
public:
	float xmin = 0;
	float xmax = 0;
	float ymin = 0;
	float ymax = 0;
	float zmin = 0;
	float zmax = 0;
	float height, w;
	int vertex_count;
	int triangle_count;
	bool light;
	Vertex *vertex;
	Vector *face_normal;
	Vector *vertex_normal;
	TriangleIndex *triangle;

	void do_construct(char *file, Transform *transform);
	float test_edge(Vector &normal, Vertex &p, Vertex &v1, Vertex &v0);
	void triangle_intersection(Ray ray, Hit &hit, int which_triangle);
	void intersection(Ray ray, Hit &hit);
	void compute_face_normal(int which_triangle, Vector &normal);
	void compute_vertex_normals(void);
	bool rayTriangleIntersect(const Ray& ray, const Vector &v0, const Vector &v1, const Vector &v2, float &t);
	PolyMesh(char *file);
	PolyMesh(char *file, Transform *transform);
	~PolyMesh(){}
};
