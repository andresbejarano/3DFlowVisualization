#include "Streamline.h"

Streamline::Streamline()
{
	this->points = new std::vector<Vector3*>();
	this->speeds = new std::vector<double>();
	this->indexInPoints = new std::vector<int>();
}

Streamline::~Streamline()
{
	this->points->clear();
	this->speeds->clear();
	delete this->points;
	delete this->speeds;
}

size_t Streamline::size()
{
	return points->size();
}

Vector3* Streamline::pointAt(int i)
{
	return points->at(i);
}

double Streamline::speedAt(int i)
{
	return speeds->at(i);
}

void Streamline::push_back(Vector3* P, double s)
{
	this->points->push_back(P);
	this->speeds->push_back(s);
}
