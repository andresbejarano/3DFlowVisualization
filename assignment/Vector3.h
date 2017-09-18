#pragma once

#include <string>

class Vector3
{
public:

	// The elements of the vector
	double x;
	double y;
	double z;

	/*
	Constructor of the class
	*/
	Vector3() : x(0), y(0), z(0) { }

	/*
	Constructor of the class
	*/
	Vector3(double x, double y, double z) : x(x), y(y), z(z) { }

	/*
	Constructor of the class
	*/
	Vector3(Vector3* V) : x(V->x), y(V->y), z(V->z) { }

	/*
	Returns the magnitude of the vector
	*/
	double magnitude();

	/*
	Normalizes the vector
	*/
	Vector3* normalize();

	/*
	Multiply the vector by the given scalar
	*/
	Vector3* multiply(double s);

	/*
	Clones the vector
	*/
	Vector3* clone();

	/*
	Adds the value of the given vector
	*/
	Vector3* add(Vector3* B);

	/*
	Subtracts the values of the given vector
	*/
	Vector3* sub(Vector3* B);

	/*
	Calculates the dot product with vector B
	*/
	double dot(Vector3* B);

	/*
	Calculates the cross product with vector B
	*/
	Vector3* cross(Vector3* B);

	/*
	Calculates the square euclidean distance to another vector
	*/
	double squareDistance(Vector3* B);

	/*
	Calculates the euclidean distance to another vector
	*/
	double distance(Vector3* B);

	/*
	Returns the intyerpolated point using parameter t and point B
	*/
	Vector3* interpolate(double t, Vector3* B);

	/*
	*/
	bool equal(Vector3* B);

	/*
	Assigns the (x, y, z) values to the vector
	*/
	Vector3* set(double x, double y, double z);

	/*
	Assigns the (x, y, z) values of V to the vector
	*/
	Vector3* set(Vector3* V);

	/*
	*/
	std::string toString();

};