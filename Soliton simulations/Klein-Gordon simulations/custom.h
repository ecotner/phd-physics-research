#ifndef CUSTOM_H_INCLUDED
#define CUSTOM_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* MATH STUFF */

/* computes the exponent of a double precision floating point number in scientific notation */
int exponent(double f);

/* computes the mantissa (significant digits) of a double precision floating point number in scientific notation */
double mantissa(double f);

/* prints a number in Mathematica readable format to a file */
void fprintM(FILE* file_name, double f);

/* STRING MANIPULATION */

/* replaces all instances of orig[] in in[] with rep[] */
void replace_string(char in[], char orig[], char rep[]);

#endif // CUSTOM_H_INCLUDED
