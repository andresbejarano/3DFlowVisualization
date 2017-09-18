#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <iostream>
#include "Vector3.h"
#include "Streamline.h"


/*
	Loads the file with the given filename and returns the reference to the generated vtkStructuredPoints object
	@param filename The name of the file
	@return A vtkSmartPointer to the vtkStructuredPoints object
*/
vtkSmartPointer<vtkStructuredPoints> loadStructuredPointsFile(char* filename)
{
	// Define the reader and load the file
	vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName(filename);
	reader->Update();

	// Return the generated object
	return reader->GetOutput();
}


/*
	Calculates a new point using trilinear interpolation
	@param P1
	@param P2
	@param P3
	@param P4
	@param P5
	@param P6
	@param P7
	@param P8
	@param u
	@param v
	@param w
	@return A Vector3 pointer to the object with the interpolated values
*/
Vector3* trilinearInterpolation(Vector3* P1, Vector3* P2, Vector3* P3, Vector3* P4, Vector3* P5, Vector3* P6, Vector3* P7, Vector3* P8, double u, double v, double w) 
{
	// Calculate the required terms for the interpolation
	Vector3* term1 = P1->clone();
	Vector3* term2 = P2->clone()->sub(P1)->multiply(u);
	Vector3* term3 = P4->clone()->sub(P1)->multiply(v);
	Vector3* term4 = P5->clone()->sub(P1)->multiply(w);
	Vector3* term5 = P1->clone()->sub(P2)->add(P3)->sub(P4)->multiply(u * v);
	Vector3* term6 = P1->clone()->sub(P2)->add(P6)->sub(P5)->multiply(u * w);
	Vector3* term7 = P1->clone()->sub(P4)->add(P8)->sub(P5)->multiply(v * w);
	Vector3* term8 = P1->clone()->sub(P2)->add(P3)->sub(P4)->add(P5)->sub(P6)->add(P7)->sub(P8)->multiply(u * v * w);

	// Calculate the interpolated value
	Vector3* result = term1->clone()->add(term2)->add(term3)->add(term4)->add(term5)->add(term6)->add(term7)->add(term8);

	// Delete pointers
	delete term1;
	delete term2;
	delete term3;
	delete term4;
	delete term5;
	delete term6;
	delete term7;
	delete term8;

	// Return the interpolated value
	return result;
}


/*
	Calculates the interpolated 3D vector for a (x, y, z) physical point using the data in the vector field
	@param vectorField The pointer to the vtkStructuredPoints object with the information of the vector field
	@param T A physical point in the vector field
	@return A Vector3 pointer with the interpolated value. This is actually a 2D vector with z = 0
*/
Vector3* getVoxel3DVector(vtkSmartPointer<vtkStructuredPoints> vectorField, Vector3* T)
{
	// Get the dimensions, spacing and origin of the vector field
	int* dimensions = vectorField->GetDimensions();
	double* spacing = vectorField->GetSpacing();
	double* origin = vectorField->GetOrigin();

	// Get the dimensions of the data structure
	double width = (double)dimensions[0];
	double height = (double)dimensions[1];
	double depth = (double)dimensions[2];

	// Get the equivalent point in terms of the data structure
	// NOTE: This is a point between (0,0,0) and (width,height,depth)
	Vector3* P = new Vector3((T->x - origin[0]) / spacing[0], (T->y - origin[1]) / spacing[1], (T->z - origin[2]) / spacing[2]);

	// Calculate the (u, v, w) values of the voxel
	// NOTE: These are not the (u, v, w) values for the trilinear interpolation but for identifying in wich voxel 
	// the current point is
	double u = ((P->x * (width - 1.0))) / width;
	double v = ((P->y * (height - 1.0))) / height;
	double w = ((P->z * (depth - 1.0))) / depth;

	// Get the coordinates for the eigth closest points in the data grid
	double u_min = floor(u);
	double v_min = floor(v);
	double w_min = floor(w);
	double u_max = ceil(u);
	double v_max = ceil(v);
	double w_max = ceil(w);

	// Get the array of points from the vector field
	// NOTE: Since there is a single sequence of points then it is safe asking for the first one only
	vtkSmartPointer<vtkDataArray> array = vectorField->GetPointData()->GetArray(0);

	// Get the eight points for the trilinear interpolation
	// NOTE: Keep in mind the axes are aligned different for the data in the array and for the points for the trilinear 
	// interpolation. For data points u is horizontal, v is vertical and w is depth. For trilineat interpolation
	// u is horizontal, v is depth and w is vertical. So, when asking for the points in the array it is follow the
	// first approach. However, for the order of the points that will be sent to the trilinear interpolation function
	// they follow the second approach. Try not to get confused with this part

	// Define P1
	int i = u_min + (v_min * width) + (w_min * width * height);
	double* data = array->GetTuple(i);
	Vector3* P1 = new Vector3(data[0], data[1], data[2]);

	// Define P2
	i = u_max + (v_min * width) + (w_min * width * height);
	data = array->GetTuple(i);
	Vector3* P2 = new Vector3(data[0], data[1], data[2]);

	// Define P3
	i = u_max + (v_max * width) + (w_min * width * height);
	data = array->GetTuple(i);
	Vector3* P3 = new Vector3(data[0], data[1], data[2]);

	// Define P4
	i = u_min + (v_max * width) + (w_min * width * height);
	data = array->GetTuple(i);
	Vector3* P4 = new Vector3(data[0], data[1], data[2]);

	// Define P5
	i = u_min + (v_min * width) + (w_max * width * height);
	data = array->GetTuple(i);
	Vector3* P5 = new Vector3(data[0], data[1], data[2]);

	// Define P6
	i = u_max + (v_min * width) + (w_max * width * height);
	data = array->GetTuple(i);
	Vector3* P6 = new Vector3(data[0], data[1], data[2]);

	// Define P7
	i = u_max + (v_max * width) + (w_max * width * height);
	data = array->GetTuple(i);
	Vector3* P7 = new Vector3(data[0], data[1], data[2]);

	// Define P8
	i = u_min + (v_max * width) + (w_max * width * height);
	data = array->GetTuple(i);
	Vector3* P8 = new Vector3(data[0], data[1], data[2]);

	// Get (u, v, w) for the trilinear interpolation
	// NOTE: Values between [0, 1]
	u = u - u_min;
	v = v - v_min;
	w = w - w_min;

	// Calculate the interpolated vector using trilinear interpolation
	Vector3* result = trilinearInterpolation(P1, P2, P3, P4, P5, P6, P7, P8, u, v, w);

	// Delete pointers
	delete P;
	delete P1;
	delete P2;
	delete P3;
	delete P4;
	delete P5;
	delete P6;
	delete P7;
	delete P8;

	// Return the result of the trilinear interpolation
	return result;
}


/*
	Fourth order Runge-Kutta (RK4) for numerical integration. Given a physical point, it returns the next physical 
	point in the stream line following the information of the vector field
	@param vectorField The pointer to the vtkStructuredPoints object with the information of the vector field
	@param T A physical point
	@param h The step value on the streamline. Negative values traverses the streamline in the negative section
	@return A Vector3 pointer with the coordinates of the next point in the streamline
*/
Vector3* RungeKutta4(vtkSmartPointer<vtkStructuredPoints> vectorField, Vector3* T, double h)
{
	// Get the dimensions, spacing and origin of the vector field
	int* dimensions = vectorField->GetDimensions();
	double* spacing = vectorField->GetSpacing();
	double* origin = vectorField->GetOrigin();

	// Calculate the physical boundary of the vector field
	double width = (dimensions[0] * spacing[0]) + origin[0];
	double height = (dimensions[1] * spacing[1]) + origin[1];
	double depth = (dimensions[2] * spacing[2]) + origin[2];

	// Get K1
	Vector3* K1 = getVoxel3DVector(vectorField, T);

	// Get K2
	Vector3* K2 = NULL;
	double xn = T->x + ((h / 2.0) * K1->x);
	double yn = T->y + ((h / 2.0) * K1->y);
	double zn = T->z + ((h / 2.0) * K1->z);
	if (xn >= origin[0] && xn < width && yn >= origin[1] && yn < height && zn >= origin[2] && zn < depth) 
	{
		Vector3* Xn = new Vector3(xn, yn, zn);
		K2 = getVoxel3DVector(vectorField, Xn);
		delete Xn;
	}
	else 
	{
		// Delete pointers and return NULL
		delete K1;
		return NULL;
	}

	// Get K3
	Vector3* K3 = NULL;
	xn = T->x + ((h / 2.0) * K2->x);
	yn = T->y + ((h / 2.0) * K2->y);
	zn = T->z + ((h / 2.0) * K2->z);
	if (xn >= origin[0] && xn < width && yn >= origin[1] && yn < height && zn >= origin[2] && zn < depth)
	{
		Vector3* Xn = new Vector3(xn, yn, zn);
		K3 = getVoxel3DVector(vectorField, Xn);
		delete Xn;
	}
	else 
	{
		// Delete pointers and return NULL
		delete K1;
		delete K2;
		return NULL;
	}

	// Get K4
	Vector3* K4 = NULL;
	xn = T->x + (h * K3->x);
	yn = T->y + (h * K3->y);
	zn = T->z + (h * K3->z);
	if (xn >= origin[0] && xn < width && yn >= origin[1] && yn < height && zn >= origin[2] && zn < depth)
	{
		Vector3* Xn = new Vector3(xn, yn, zn);
		K4 = getVoxel3DVector(vectorField, Xn);
		delete Xn;
	}
	else 
	{
		// Delete pointers and return NULL
		delete K1;
		delete K2;
		delete K3;
		return NULL;
	}

	// Calculate the next coordinate
	Vector3* term1 = K1->clone();
	Vector3* term2 = K2->clone()->multiply(2.0);
	Vector3* term3 = K3->clone()->multiply(2.0);
	Vector3* term4 = K4->clone();
	Vector3* Xnext = term1->clone()->add(term2)->add(term3)->add(term4)->multiply(h / 6.0)->add(T);

	// Delete pointers
	delete K1;
	delete K2;
	delete K3;
	delete K4;
	delete term1;
	delete term2;
	delete term3;
	delete term4;

	// Return the result
	return Xnext;
}


/*
	Returns a section (positive or negative) of the streamline starting on the given point
	@param vectorField The pointer to the vtkStructuredPoints object with the information of the vector field
	@param T The origin point of the streamline section. This is a physical point
	@param h The step value on the streamline. Negative values traverses the streamline in the negative section
	@return The streamline object with the information of the streamline section
*/
Streamline* getStreamlineSection(vtkSmartPointer<vtkStructuredPoints> vectorField, Vector3* T, double h) 
{
	// Get the dimensions, spacings and origin of the vector field
	int* dimensions = vectorField->GetDimensions();
	double* spacing = vectorField->GetSpacing();
	double* origin = vectorField->GetOrigin();

	// Get the physical boundary values of the vector field
	double width = (dimensions[0] * spacing[0]) + origin[0];
	double height = (dimensions[1] * spacing[1]) + origin[1];
	double depth = (dimensions[2] * spacing[2]) + origin[2];

	// Initialize the object where the streamline section will be stored
	Streamline* streamlineSection = new Streamline();

	// Define the max length for a streamline section
	int n = 600;
	int length = 0;

	// Store the coordinates of the initial point
	int x0 = T->x;
	int y0 = T->y;
	int z0 = T->z;

	// Traverse through the streamline until it reaches a boundary or its maximum length
	while (length < n) 
	{
		// If the current point lies inside the physical boundary of the vector field then proceed; otherwise, 
		// stop the streamline section generation
		if (x0 >= origin[0] && x0 <= width && y0 >= origin[1] && y0 <= height && z0 >= origin[2] && z0 <= depth) 
		{
			// The current point is valid. Generate a Vector3 object for it
			Vector3* current = new Vector3(x0, y0, z0);

			// The speed at this point is the interpolation of the current point with respect of the eight 
			// neighboring points in the data structure. Calculate such interpolated vector, get its magnitude
			// and delete the pointer to the vector (it is not longer needed)
			Vector3* speedVector = getVoxel3DVector(vectorField, current);
			double speed = speedVector->magnitude();
			delete speedVector;

			// Insert the current point and its speed in the streamline section
			streamlineSection->push_back(current, speed);
			length += 1;

			// Get the next point in the vector field
			Vector3* next = RungeKutta4(vectorField, current, h);

			// If a next location is given then process is; otherwise, stop the cycle
			if (next != NULL) 
			{
				if (x0 == next->x && y0 == next->y && z0 == next->z) 
				{
					//std::cout << "+++++++++++++++++++++++++++ EQUAL POINT +++++++++++++++++++++++++++" << std::endl;

					// Delete pointer and stop the cycle
					delete next;
					break;
				}
				else 
				{
					// Store the coordinates for the next point in the streamline
					x0 = next->x;
					y0 = next->y;
					z0 = next->z;

					// Delete pointer
					delete next;
				}
			}
			else 
			{
				//std::cout << "+++++++++++++++++++++++++++ NULL POINTER +++++++++++++++++++++++++++" << std::endl;

				// Delete pointer and stop the cycle
				delete next;
				break;
			}
		}
		else 
		{
			//std::cout << "+++++++++++++++++++++++++++ OUT OF BOUNDS +++++++++++++++++++++++++++" << std::endl;

			break;
		}
	}

	// Return the streamline section
	return streamlineSection;
}


/*
	Returns the indicated streamline with origin in the given point
	@param vectorField The pointer to the vtkStructuredPoints object with the information of the vector field
	@param T The origin point of the streamline section. This is a physical point
	@param h The length of the step
	@return The streamline object
*/
Streamline* getStreamline(vtkSmartPointer<vtkStructuredPoints> vectorField, Vector3* T, double h)
{
	// Get the positive and negative sections of the streamline
	Streamline* positiveSection = getStreamlineSection(vectorField, T, h);
	Streamline* negativeSection = getStreamlineSection(vectorField, T, -h);

	// Define the streamline object
	Streamline* streamline = new Streamline();

	// Get the size of the negative section of the streamline
	int n = negativeSection->size();

	// Put the elements of the negative streamline section into the main streamline object
	for (int i = n - 1; i >= 0; i -= 1) 
	{
		// Copy the current element into the main streamline object
		// NOTE: Here the point vector is cloned instead of just storing its pointer. The reason is the negative 
		// streamline section object is later destroyed.
		streamline->push_back(negativeSection->pointAt(i)->clone(), negativeSection->speedAt(i));
	}

	// get the size of the positive section of the streamline
	n = positiveSection->size();

	// Put the elements of the positive section into the main streamline object
	// NOTE: The first element is not inserted since it already came from the negative section
	for (int i = 1; i < n; i += 1) 
	{
		// Copy the current element into the main streamline object
		// NOTE: Here the point vector is cloned instead of just storing its pointer. The reason is the positive 
		// streamline section object is later destroyed.
		streamline->push_back(positiveSection->pointAt(i)->clone(), positiveSection->speedAt(i));
	}
	
	// Delete pointers
	delete positiveSection;
	delete negativeSection;

	// Returns the streamline
	return streamline;
}


/*
	Returns the streamlines generated along the main diagonal of the vector field
	@param vectorField The pointer to the vtkStructuredPoints object with the information of the vector field
	@param n The number of seeds (points) on the main diagonal where the streamlines will be generated
	@param h The length of the step
	@return A vector with the information of the generated streamlines
*/
std::vector<Streamline*>* getStreamlines(vtkSmartPointer<vtkStructuredPoints> vectorField, Vector3* P, Vector3* Q, int n, double h)
{
	// Define the direction of the diagonal segment (from start to end)
	Vector3* D = Q->clone()->sub(P);

	// Initialize the vector where the streamlines will be stored
	std::vector<Streamline*>* streamlines = new std::vector<Streamline*>();

	// Initialize the auxiliary vector pointers for the streamlines generation
	Vector3* tD = new Vector3();
	Vector3* T = new Vector3();

	// Calculate the step for the seed values along the diagonal line
	double step = 1.0 / (double)n;

	// Initialize the parameter of the diagonal
	double t = 0.0;

	// Generate the streamline for each one of the points in the diagonal
	for (int i = 0; i < n; i += 1) 
	{
		// Calculate the current point
		// NOTE: This is a physical point, not a data structure point!!
		tD->set(D)->multiply(t);
		T->set(P)->add(tD);

		// Get the streamline and push it into the streamlines vector
		streamlines->push_back(getStreamline(vectorField, T, h));

		// Set the value of t for the next point in the diagonal
		t += step;
	}

	// Delete pointers
	delete P;
	delete Q;
	delete D;
	delete tD;
	delete T;

	// Returns the streamlines
	return streamlines;
}


/*
	Returns the total number of points in all of the streamlines
	@param streamlines A pointer to the vector where the streamlines are stored
	@return The number of points in all of the streamlines
*/
int getNumPoints(std::vector<Streamline*>* streamlines) 
{
	// Get the number of stored streamlines
	int n = streamlines->size();

	// Initialize the counter of points
	int count = 0;

	// Traverse through the streamlines and accumulate their size
	for (int i = 0; i < n; i += 1) 
	{
		count += streamlines->at(i)->size();
	}

	// Return the counted number of points
	return count;
}


/*
	Generate the stream surface based on the information of the streamlines. Each triangle
	ribbon is generated using the greedy approach described in [Hulquist 1992]
	@param streamlines The pointer to the vector with the streamlines
	@return A pointer to the vtkCellArray object with the information of the triangles
*/
vtkSmartPointer<vtkCellArray> buildStreamSurface(std::vector<Streamline*>* streamlines) 
{
	// Initialize the object where the stream surface triangles will be stored
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	// Get the number of stored streamlines
	int n = streamlines->size();

	// Get the pointer to the first streamline
	// NOTE: This is the pointer to the left streamline
	Streamline* L = streamlines->at(0);

	// Traverse through the streamlines, starting from index 1
	for (int i = 1; i < n; i += 1) 
	{
		// Let's build the triangle ribbon between the current streamline and the previous one

		// Get the pointer to the second streamline
		// NOTE: This is the pointer to the right streamline
		Streamline* R = streamlines->at(i);

		// Get the size of each streamline
		int nLeft = L->size();
		int nRight = R->size();

		// Get the first point for each streamline
		// NOTE: It is assumed the streamlines have at least one point
		Vector3* L0 = L->pointAt(0);
		Vector3* R0 = R->pointAt(0);

		// Define the index values in the respective streamlines for L0, L1, R0 and R1 points
		int iL0 = 0;
		int iL1 = 1;
		int iR0 = 0;
		int iR1 = 1;

		// Get the second point on each streamline
		// NOTE: If it doesn't exist then set the pointer to NULL
		Vector3* L1 = (nLeft == 1) ? NULL : L->pointAt(1);
		Vector3* R1 = (nRight == 1) ? NULL : R->pointAt(1);

		// Generate the ribbons until one of the second points is NULL
		// NOTE: This means the end of one streamline has been reached; therefore, L0 or R0 has 
		// the pointer to the last point in the streamline
		while (L1 != NULL && R1 != NULL) 
		{
			// Calculate the distance of the diagonals (L0-R1 and L1-R0)
			// NOTE: The square distance is faster. Calculating the square root is not actually necessary for 
			// determining which distance is the longest one
			double dL0R1 = L0->squareDistance(R1);
			double dL1R0 = L1->squareDistance(R0);

			// If the first distance is shorter or equal than the second distance then define triangle L0 R1 R0
			if (dL0R1 <= dL1R0) 
			{
				// Define the triangle L0 R1 R0
				vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
				triangle->GetPointIds()->SetId(0, L->indexInPoints->at(iL0));
				triangle->GetPointIds()->SetId(1, R->indexInPoints->at(iR1));
				triangle->GetPointIds()->SetId(2, R->indexInPoints->at(iR0));

				// Add triangle to triangles object
				triangles->InsertNextCell(triangle);

				// Update pointer and indexes for the right streamline
				R0 = R1;
				R1 = (iR1 + 1 < nRight) ? R->pointAt(iR1 + 1) : NULL;
				iR0 += 1;
				iR1 += 1;
			}
			else 
			{
				// Define the triangle R0 L1 L0
				vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
				triangle->GetPointIds()->SetId(0, R->indexInPoints->at(iR0));
				triangle->GetPointIds()->SetId(1, L->indexInPoints->at(iL1));
				triangle->GetPointIds()->SetId(2, L->indexInPoints->at(iL0));

				// Add triangle to triangles object
				triangles->InsertNextCell(triangle);

				// Update pointer and indexes for the left streamline
				L0 = L1;
				L1 = (iL1 + 1 < nLeft) ? L->pointAt(iL1 + 1) : NULL;
				iL0 += 1;
				iL1 += 1;
			}
		}

		// If both L1 and R1 are NULL then nothing to do; otherwise, a triangle fan has to be generated
		// NOTE: Actually, this would happen if both of the streamlines have only one point, a very weird case
		if (!(L1 == NULL && R1 == NULL)) 
		{
			// Finish the ribbon for the remaining points
			// NOTE: This is technically generating a triangle fan

			// If L1 is NULL it means L0 is in the last point in the left streamline. Then, connect all of the 
			// remaining points in the right streamline to this point
			if (L1 == NULL) 
			{
				// Traverse through the right streamline and finish the ribbon with a triangle fan
				while (R1 != NULL) 
				{
					// Define the triangle L0 R1 R0
					vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
					triangle->GetPointIds()->SetId(0, L->indexInPoints->at(iL0));
					triangle->GetPointIds()->SetId(1, R->indexInPoints->at(iR1));
					triangle->GetPointIds()->SetId(2, R->indexInPoints->at(iR0));

					// Add triangle to triangles object
					triangles->InsertNextCell(triangle);

					// Update pointer and indexes for the right streamline
					R0 = R1;
					R1 = (iR1 + 1 < nRight) ? R->pointAt(iR1 + 1) : NULL;
					iR0 += 1;
					iR1 += 1;
				}
			}
			else 
			{
				// It means R1 is NULL and R0 is in the last point in the right streamline. Then, connect all of
				// the remaining points in the left streamline to this point

				// Traverse through the left streamline and finish the ribbon with a triangle fan
				while (L1 != NULL) 
				{
					// Define the triangle R0 L1 L0
					vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
					triangle->GetPointIds()->SetId(0, R->indexInPoints->at(iR0));
					triangle->GetPointIds()->SetId(1, L->indexInPoints->at(iL1));
					triangle->GetPointIds()->SetId(2, L->indexInPoints->at(iL0));

					// Add triangle to triangles object
					triangles->InsertNextCell(triangle);

					// Update pointer and indexes for the right streamline
					L0 = L1;
					L1 = (iL1 + 1 < nLeft) ? L->pointAt(iL1 + 1) : NULL;
					iL0 += 1;
					iL1 += 1;
				}
			}
		}

		// Set the current right streamline as the next left streamline
		L = R;
	}

	// Returns the generated triangles
	return triangles;
}


/*
	The main function
*/
int main(int argc, char** argv)
{
	// Set the length of the step
	double h = 0.15;
	std::cout << "h = ";
	std::cin >> h;

	// Set the number of seeds along the main diagonal for finding the streamlines
	int nSeeds = 100;
	std::cout << "Seeds: ";
	std::cin >> nSeeds;

	// Indicate whether the stream surface has to be generated
	bool generateSurface = true;
	std::cout << "Stream Surface? (0=No, 1=Yes): ";
	std::cin >> generateSurface;

	// Load the vector field information
	vtkSmartPointer<vtkStructuredPoints> vectorField = loadStructuredPointsFile("delta.vtk");

	// Get the dimensions, spacing and origin of the physical vector field
	int* dimensions = vectorField->GetDimensions();
	double* spacing = vectorField->GetSpacing();
	double* origin = vectorField->GetOrigin();

	// Find the physical start and end point of the vector field
	// NOTE: This two points correspond to a diagonal that passes through the vector field from origin corner to opposite corner
	Vector3* P = new Vector3(origin[0], origin[1], origin[2]);
	Vector3* Q = new Vector3((dimensions[0] * spacing[0]) + origin[0], (dimensions[1] * spacing[1]) + origin[1], (dimensions[2] * spacing[2]) + origin[2]);

	std::cout << std::endl;
	std::cout << "Origin point found: " << P->toString() << std::endl;
	std::cout << "Opposite point found: " << Q->toString() << std::endl;

	// Indicate whether the diagonal will be used for generating the streamlines
	bool useDiagonal = true;

	// Ask whether the diagonal will be used or not
	std::cout << std::endl;
	std::cout << "Do you want to use the diagonal for finding the streamlines? (0=No, 1=Yes): ";
	std::cin >> useDiagonal;

	// If the diagonal is not going to be used then ask for the start and end point
	if (!useDiagonal) 
	{
		// Initialize the coordinates for the starting point
		double px;
		double py;
		double pz;

		// Initialize the coordinates for the ending point
		double qx;
		double qy;
		double qz;

		// Ask for the start point coordinates
		std::cout << "Start point coordinates (e.g., 20, 31, 42): ";
		std::cin >> px >> py >> pz;

		// Ask for the end point coordinates
		std::cout << "End point coordinates (e.g., 56, 78, 90): ";
		std::cin >> qx >> qy >> qz;

		// Update the start and end points
		P->set(px, py, pz);
		Q->set(qx, qy, qz);
	}

	// Initialize the object where the streamlines points will be stored
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Initialize the object where the streamlines lines will be stored
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	// Set an indicator for processing a new streamline
	// NOTES: This is required for generating the lines of a streamline. If it is set to true it means 
	// only the first points has been inserted into the visual line; otherwise, we have at least two 
	// points and therefore we can define a line segment in the streamline
	bool newStreamline = true;

	// Get the streamlines
	std::vector<Streamline*>* streamlines = getStreamlines(vectorField, P, Q, nSeeds, h);

	// Initialize the variables where the max and min speeds will be stored
	double maxSpeed = 0;
	double minSpeed = 0;

	// Indicate the next speed value to be processes is the first one
	// NOTE: This is done for obtaining the max and min speed values in the vector field
	bool firstSpeed = true;

	// Initialize the index of the points
	// NOTE: This is required for defining the lines
	int pointIndex = 0;

	// Initialize the array where the color of the points will be stored. Then, set its size with the 
	// total number of points in all of the streamlines plus eight (for the boundary lines)
	vtkSmartPointer<vtkDoubleArray> pointColors = vtkSmartPointer<vtkDoubleArray>::New();
	pointColors->SetNumberOfValues(getNumPoints(streamlines) + 8);
	
	// Traverse each one of the streamlines and generate them
	int nStreamlines = streamlines->size();
	for (int i = 0; i < nStreamlines; i += 1) 
	{
		// Get the pointer to the current streamline
		Streamline* streamline = streamlines->at(i);

		// Get the number of points in the streamline
		int nPoints = streamline->size();

		// Indicate a new streamline is going to be processed
		newStreamline = true;

		// Add the points in the visual streamline, generating the lines in the process
		for (int j = 0; j < nPoints; j += 1) 
		{
			// Get the pointer to the current point in the streamline and store it in the points object
			Vector3* P = streamline->pointAt(j);
			points->InsertNextPoint(P->x, P->y, P->z);

			// Store its speed for coloring purposes
			pointColors->SetValue(pointIndex, streamline->speedAt(j));

			// Store the index of the point in the points object
			// NOTE: This is required for building the stream surface
			streamline->indexInPoints->push_back(pointIndex);

			// Increase the index for the next point
			pointIndex += 1;

			// If this is the first stored speed then set it for both max and min; otherwise, update variables respectively
			if (firstSpeed) 
			{
				maxSpeed = streamline->speedAt(j);
				minSpeed = streamline->speedAt(j);
				firstSpeed = false;
			}
			else 
			{
				if (streamline->speedAt(j) > maxSpeed) 
				{
					maxSpeed = streamline->speedAt(j);
				}
				else if (streamline->speedAt(j) < minSpeed) 
				{
					minSpeed = streamline->speedAt(j);
				}
			}

			// If it is not a new streamline then add the line respecting the last two points
			if (!newStreamline) 
			{
				// Initialize the current line, define it and store it in the lines object
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, pointIndex - 2);
				line->GetPointIds()->SetId(1, pointIndex - 1);
				lines->InsertNextCell(line);
			}

			// Indicate it is no longer a new streamline since it already has a point in the points object
			newStreamline = false;
		}
	}

	// Add the boundary points and generate the boundary lines
	points->InsertNextPoint(origin[0], origin[1], origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint((dimensions[0] * spacing[0]) + origin[0], origin[1], origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint((dimensions[0] * spacing[0]) + origin[0], origin[1], (dimensions[2] * spacing[2]) + origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint(origin[0], origin[1], (dimensions[2] * spacing[2]) + origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint(origin[0], (dimensions[1] * spacing[1]) + origin[1], origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint((dimensions[0] * spacing[0]) + origin[0], (dimensions[1] * spacing[1]) + origin[1], origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint((dimensions[0] * spacing[0]) + origin[0], (dimensions[1] * spacing[1]) + origin[1], (dimensions[2] * spacing[2]) + origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	points->InsertNextPoint(origin[0], (dimensions[1] * spacing[1]) + origin[1], (dimensions[2] * spacing[2]) + origin[2]);
	pointColors->SetValue(pointIndex, maxSpeed);
	pointIndex += 1;

	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 8);
	line->GetPointIds()->SetId(1, pointIndex - 7);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 7);
	line->GetPointIds()->SetId(1, pointIndex - 6);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 6);
	line->GetPointIds()->SetId(1, pointIndex - 5);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 5);
	line->GetPointIds()->SetId(1, pointIndex - 8);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 4);
	line->GetPointIds()->SetId(1, pointIndex - 3);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 3);
	line->GetPointIds()->SetId(1, pointIndex - 2);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 2);
	line->GetPointIds()->SetId(1, pointIndex - 1);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 1);
	line->GetPointIds()->SetId(1, pointIndex - 4);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 8);
	line->GetPointIds()->SetId(1, pointIndex - 4);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 5);
	line->GetPointIds()->SetId(1, pointIndex - 1);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 7);
	line->GetPointIds()->SetId(1, pointIndex - 3);
	lines->InsertNextCell(line);

	line = vtkSmartPointer<vtkLine>::New();
	line->GetPointIds()->SetId(0, pointIndex - 6);
	line->GetPointIds()->SetId(1, pointIndex - 2);
	lines->InsertNextCell(line);

	// Add the diagonal line (TO BE REMOVED)
	//line = vtkSmartPointer<vtkLine>::New();
	//line->GetPointIds()->SetId(0, pointIndex - 2);
	//line->GetPointIds()->SetId(1, pointIndex - 1);
	//lines->InsertNextCell(line);

	// If the stream surface has to be displayed then generate it; otherwise, set this object to NULL
	vtkSmartPointer<vtkCellArray> triangles = generateSurface ? buildStreamSurface(streamlines) : NULL;
	
	// Initialize the object where all of the geometry will be added. Then, insert points and lines
	vtkSmartPointer<vtkPolyData> linesPolydata = vtkSmartPointer<vtkPolyData>::New();
	linesPolydata->SetPoints(points);
	linesPolydata->SetLines(lines);
	linesPolydata->GetPointData()->SetScalars(pointColors);

	// If a stream surface was built then insert the triangles
	if (generateSurface)
	{
		linesPolydata->SetPolys(triangles);
	}

	// Define the color function for visualizing the speed on the vector field
	vtkSmartPointer<vtkColorTransferFunction> colorFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	colorFunction->AddRGBPoint(maxSpeed, 1.0, 1.0, 1.0);	// White for max speed
	colorFunction->AddRGBPoint(minSpeed, 0.0, 0.0, 1.0); 	// Blue for low speed

	// Initialize the mapper and define it
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(linesPolydata);
	mapper->SetLookupTable(colorFunction);

	// Initialize the actor and define it
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// Initialize the renderer and define it
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->SetBackground(0.1, 0.1, 0.1);

	// Initialize the render window and define it
	vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
	window->AddRenderer(renderer);
	window->SetSize(600, 600);
	window->Render();

	// Initialize the render window interactor and define it
	vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(window);
	interactor->Initialize();
	interactor->Start();

	return EXIT_SUCCESS;
}
