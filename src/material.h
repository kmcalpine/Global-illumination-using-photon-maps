/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Phong is a child class of Material and implement the simple Phong
// surface illumination model.

#pragma once

#include "colour.h"
#include "vector.h"

class Material {
public:
	bool light;
	bool reflect;
	bool transparent;
	bool diffuseSurface;
	Colour ambient;
	float ai;
	Colour diffuse;
	Colour specular;
	float power;
	float ri; //refractive index

};



