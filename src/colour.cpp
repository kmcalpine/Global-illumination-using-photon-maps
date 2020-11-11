/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
comment
 */

#include "colour.h"

Colour::Colour() {
    r = 0.0;
    g = 0.0;
    b = 0.0;
}

Colour::Colour(float red, float, green, float blue) {
    this.r = red;
    this.g = green;
    this.b = blue;
}
