#include "Vector3.h"
#include <math.h>
#include <sstream>

double Vector3::magnitude()
{
	// Calculate the magnitude of the vector and return it
	return sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
}

Vector3* Vector3::normalize()
{
	// Get the magnitude of the vector
	double mag = this->magnitude();

	// If it has a magnitude then normalize it
	if (mag > 0)
	{
		this->x /= mag;
		this->y /= mag;
		this->z /= mag;
		return this;
	}

	// Since the vector has no magnitude then return a zero vector
	return this;
}

Vector3* Vector3::multiply(double s)
{
	// Multiply each component by the scalar
	this->x *= s;
	this->y *= s;
	this->z *= s;
	return this;
}

Vector3* Vector3::clone()
{
	return new Vector3(this->x, this->y, this->z);
}

Vector3* Vector3::add(Vector3* B)
{
	this->x += B->x;
	this->y += B->y;
	this->z += B->z;
	return this;
}

Vector3* Vector3::sub(Vector3* B)
{
	this->x -= B->x;
	this->y -= B->y;
	this->z -= B->z;
	return this;
}

double Vector3::dot(Vector3* B)
{
	return (this->x * B->x) + (this->y * B->y) + (this->z * B->z);
}

Vector3* Vector3::cross(Vector3* B)
{
	// Calculate the cross product values
	double i = (this->y * B->z) - (this->z * B->y);
	double j = (this->z * B->x) - (this->x * B->z);
	double k = (this->x * B->y) - (this->y * B->x);

	this->x = i;
	this->y = j;
	this->z = k;
	return this;
}

double Vector3::squareDistance(Vector3* B)
{
	return ((B->x - this->x) * (B->x - this->x)) + 
		   ((B->y - this->y) * (B->y - this->y)) + 
		   ((B->z - this->z) * (B->z - this->z));
}

double Vector3::distance(Vector3* B)
{
	return sqrt(squareDistance(B));
}

Vector3* Vector3::interpolate(double t, Vector3* B)
{
	return this->multiply(1.0 - t)->add(B->multiply(t));
}

bool Vector3::equal(Vector3* B)
{
	return (this->x == B->x && this->y == B->y && this->z == B->z);
}

Vector3* Vector3::set(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	return this;
}

Vector3* Vector3::set(Vector3* V)
{
	this->x = V->x;
	this->y = V->y;
	this->z = V->z;
	return this;
}

std::string Vector3::toString()
{
	// Convert the vector to a string notation and return it
	std::stringstream ss;
	ss << "(" << this->x << ", " << this->y << ", " << this->z << ")";
	return ss.str();
}
