#ifndef CIRCLE_H_
#define CIRCLE_H_

#include <iomanip>
#include <iostream>
#include <math.h>

using namespace std;

template<typename T> inline T SQR(T t) { return (t)*(t); };


/************************************************************************
			DECLARATION OF THE CLASS DATA
************************************************************************/
// Class for Data
// A data has 5 fields:
//       n (of type int), the number of data points
//       X and Y (arrays of type double), arrays of x- and y-coordinates
//       meanX and meanY (of type double), coordinates of the centroid (x and y sample means)

class Data
{
public:

	int n;
	double *X;		//space is allocated in the constructors
	double *Y;		//space is allocated in the constructors
	double meanX, meanY;

	// constructors
	Data();
	Data(int N);
	Data(int N, double X[], double Y[]);

	// routines
	void means(void);
	void center(void);
	void scale(void);
	void print(void);

	// destructors
	~Data();
};


/************************************************************************
			DECLARATION OF THE CLASS CIRCLE
************************************************************************/
// Class for Circle
// A circle has 7 fields:
//     a, b, r (of type reals), the circle parameters
//     s (of type reals), the estimate of sigma (standard deviation)
//     g (of type reals), the norm of the gradient of the objective function
//     i and j (of type int), the iteration counters (outer and inner, respectively)

class Circle
{
public:

	// The fields of a Circle
	double a, b, r, s, g, Gx, Gy, variance;
	int i, j;

	// constructors
	Circle();
	Circle(double aa, double bb, double rr);

	// routines
	void print(void);

	// no destructor we didn't allocate memory by hand.
};


#endif /* CIRCLE_H_ */
