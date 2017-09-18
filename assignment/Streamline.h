#pragma once

#include <vector>
#include "Vector3.h"

/*
	A class representing a streamline in 3D
	A streamline is considered a sequen of points in a line (curve) through a vector field.

	NOTE: This class is for storing data purposes. No significant logic is implemented here.
*/
class Streamline 
{
	public:

		// The points in the streamline
		std::vector<Vector3*>* points;

		// The speed values for each point in the streamline
		std::vector<double>* speeds;

		// The location of each point in the streamline in the points object
		// NOTE: This is required for stream surface visualization
		std::vector<int>* indexInPoints;

		/*
			Constructor of the class
		*/
		Streamline();

		/*
			Destructor of the class
		*/
		~Streamline();

		/*
			Returns the size of the streamline (the number of points)
		*/
		size_t size();

		/*
			Returns the pointer to the vector at the given position
		*/
		Vector3* pointAt(int i);

		/*
			Returns the speed of the vector at the given position
		*/
		double speedAt(int i);

		/*
			Inserts the given point and speed value into the streamlines vectors
		*/
		void push_back(Vector3* P, double s);
};
