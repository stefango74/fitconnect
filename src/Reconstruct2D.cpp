// Author      : Stefan Ohrhallinger
// Version     :
// Copyright   : GPL v3
// Description : Reconstruct a curve from noisy 2D points
//============================================================================

#include <list>
#include <map>
#include <fstream>
#include "Reconstruct2D.h"

using namespace std;
using namespace casa;

const bool DEBUG = false;
//const bool DEBUG = true;

#define SHARPLEAF

//****************** Sigma ************************************
//
//   estimate of Sigma = square root of RSS divided by N
//   gives the root-mean-square error of the geometric circle fit

double Sigma (Data& data, Circle& circle)
{
	double sum=0.,dx,dy;

    for (int i=0; i<data.n; i++)
    {
        dx = data.X[i] - circle.a;
        dy = data.Y[i] - circle.b;
        sum += SQR(sqrt(dx*dx+dy*dy) - circle.r);
    }
    return sqrt(sum/data.n);
}

Circle CircleFitByHyper (Data& data)
/*
      Circle fit to a given set of data points (in 2D)

      This is an algebraic fit based on the journal article

      A. Al-Sharadqah and N. Chernov, "Error analysis for circle fitting algorithms",
      Electronic Journal of Statistics, Vol. 3, pages 886-911, (2009)

      It is an algebraic circle fit with "hyperaccuracy" (with zero essential bias).
      The term "hyperaccuracy" first appeared in papers by Kenichi Kanatani around 2006

      Input:  data     - the class of data (contains the given points):

	      data.n   - the number of data points
	      data.X[] - the array of X-coordinates
	      data.Y[] - the array of Y-coordinates

     Output:
               circle - parameters of the fitting circle:

	       circle.a - the X-coordinate of the center of the fitting circle
	       circle.b - the Y-coordinate of the center of the fitting circle
 	       circle.r - the radius of the fitting circle
 	       circle.s - the root mean square error (the estimate of sigma)
 	       circle.j - the total number of iterations

     This method combines the Pratt and Taubin fits to eliminate the essential bias.

     It works well whether data points are sampled along an entire circle or
     along a small arc.

     Its statistical accuracy is theoretically higher than that of the Pratt fit
     and Taubin fit, but practically they all return almost identical circles
     (unlike the Kasa fit that may be grossly inaccurate).

     It provides a very good initial guess for a subsequent geometric fit.

       Nikolai Chernov  (September 2012)

*/
{
    int i,iter,IterMAX=99;

    double Xi,Yi,Zi;
    double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
    double A0,A1,A2,A22;
    double Dy,xnew,x,ynew,y;
    double DET,Xcenter,Ycenter;

    Circle circle;

    data.means();   // Compute x- and y- sample means (via a function in the class "data")

//     computing moments

    Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

    for (i=0; i<data.n; i++)
    {
        Xi = data.X[i] - data.meanX;   //  centered x-coordinates
        Yi = data.Y[i] - data.meanY;   //  centered y-coordinates
        Zi = Xi*Xi + Yi*Yi;

        Mxy += Xi*Yi;
        Mxx += Xi*Xi;
        Myy += Yi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
        Mzz += Zi*Zi;
    }
    Mxx /= data.n;
    Myy /= data.n;
    Mxy /= data.n;
    Mxz /= data.n;
    Myz /= data.n;
    Mzz /= data.n;

//    computing the coefficients of the characteristic polynomial

    Mz = Mxx + Myy;
    Cov_xy = Mxx*Myy - Mxy*Mxy;
    Var_z = Mzz - Mz*Mz;

    A2 = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    A1 = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    A22 = A2 + A2;

//    finding the root of the characteristic polynomial
//    using Newton's method starting at x=0
//     (it is guaranteed to converge to the right root)

	for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
        Dy = A1 + x*(A22 + 16.*x*x);
        xnew = x - y/Dy;
        if ((xnew == x)||(!isfinite(xnew))) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
        if (fabs(ynew)>=fabs(y))  break;
        x = xnew;  y = ynew;
    }

//    computing parameters of the fitting circle

    DET = x*x - x*Mz + Cov_xy;

    // test for degenerate case (near-collinear points)
    if (DET < 1e-20)
    {
    	// construct circle with radius ~variance
    	Xcenter = (data.Y[1] - data.Y[0])/sqrt(Zi);
    	Ycenter = -(data.X[1] - data.X[0])/sqrt(Zi);
    }
    else
    {
		Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
		Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;
    }

//       assembling the output

    circle.a = Xcenter + data.meanX;
    circle.b = Ycenter + data.meanY;
    circle.r = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz - x - x);
    circle.s = Sigma(data,circle);
    circle.i = 0;
    circle.j = iter;  //  return the number of iterations, too

    return circle;
}

Reconstruct2D::Reconstruct2D(const vector<Point> &p_points, int p_mode)
{
	int i;

	points = p_points;
	mode = p_mode;
	maxIter = -1;

	for (i = 0; i < (int)points.size(); i++)
		pClasses.push_back(CONFORMING);
}

Reconstruct2D::Reconstruct2D(const vector<Point> &p_points, const vector<float> &p_noise, int p_mode)
{
	int i;

	points = p_points;
	noise = p_noise;
	mode = p_mode;
	maxIter = -1;

	for (i = 0; i < (int)points.size(); i++)
		pClasses.push_back(CONFORMING);
}

void Reconstruct2D::setMaxIter(int p_maxIter)
{
	maxIter = p_maxIter;
}

void Reconstruct2D::getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtime)
{
	p_output = outputSize;
	p_iter = iterations;
	p_fit = handledFitCount;
	p_point = handledPointCount;
	p_squared = squaredFitCount;
	p_runtime = runtime;
}

enum class PState: char { INITIAL, OUTLIER, CONFORM, NONCONFORM, UNFITTED, FITTED0, FITTED1, MANIFOLD, SHARP, LEAF, ELIMINATED };

/*
 * intersect circle with line (p, vec) and return intersections in result0 and result1
 */
bool intersectCircleLine(Circle &circle, Point &p, Point &vec, Point &result0, Point &result1)
{
	float a = -vec[1];
	float b = vec[0];
	float c = a*(p[0] - circle.a) + b*(p[1] - circle.b);	// translate point into coordinate system centered at circle center
	float d2 = SQR(a) + SQR(b);
	float det = SQR(circle.r)*d2 - SQR(c);

	if (det < 0.0)
		return false;

	float temp = b*sqrt(det);
	result0[0] = circle.a + (a*c + temp)/d2;
	result1[0] = circle.a + (a*c - temp)/d2;
	temp = a*sqrt(det);
	result0[1] = circle.b + (b*c - temp)/d2;
	result1[1] = circle.b + (b*c + temp)/d2;

	return true;
}

/*
 * convert a vector into a degree angle (0-360)
 */
float computeDegreesFromVector(Point vec)
{
	float angle = atan2(vec.y(), vec.x());

	if (angle < -PI)
		angle += 2*PI;
	else
	if (angle > PI)
		angle -= 2*PI;

	return angle;
}

void projectPoints(vector<Point> &points, vector<PState> &pStates, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, vector<Point> &normals, vector<Circle> &circles, vector<Point> &projPoints)
{
	int i, j, l;

	if (DEBUG)
		cout << "blending points by projections on referencing fits" << endl;

	// blend points along their normal with all referencing fitted circles
	vector<int> weights(points.size());

	for (i = 0; i < (int)points.size(); i++)
		if ((pStates[i] == PState::SHARP) || (pStates[neighbors[i][0]] == PState::SHARP) || (pStates[neighbors[i][1]] == PState::SHARP))
		{
			// keep original point (will not be modified)
			weights[i] = 1;
			projPoints[i] = points[i];
		}
		else
		if (pStates[i] == PState::MANIFOLD)
		{
			// initialize with projected point
			weights[i] = 1;
			projPoints[i] = Point(circles[i].a, circles[i].b) + normals[i]*circles[i].r;
		}
		else
		{
			weights[i] = 0;
			projPoints[i] = Point(0.0, 0.0);
		}

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] == PState::MANIFOLD)
		{
			for (j = 0; j < 2; j++)
			{
				for (auto k:nhood[j][i])
				{
					if ((pStates[k] != PState::SHARP) && (pStates[neighbors[k][0]] != PState::SHARP) && (pStates[neighbors[k][1]] != PState::SHARP))
					{
						// intersect normal of point k with circle fit of point i
						Point pp[2];

						// DEBUG
						if (k == 4)
							cout << "";

						// if there are intersections
						if (intersectCircleLine(circles[i], points[k], normals[k], pp[0], pp[1]))
						{
							// test intersection with arc - determine points on circle in degrees
							const int COUNT = 6;
							float deg[COUNT];	// start/end/i/k/pp0/pp1
							Point center(circles[i].a, circles[i].b);
							Point pNormal[COUNT] = { points[nhood[0][i].back()] - center, points[nhood[1][i].back()] - center,
												points[i] - center, points[k] - center,
												pp[0] - center, pp[1] - center };

							for (l = 0; l < COUNT; l++)
							{
								pNormal[l].normalize();
								deg[l] = ((computeDegreesFromVector(pNormal[l])+PI)/PI*180);
							}

							// reverse start/end of arc if not oriented in positive direction
							if (fmod(deg[2] + 360.0 - deg[0], 360.0) > fmod(deg[1] + 360.0 - deg[0], 360.0))
								swap(deg[0], deg[1]);

							// standardize angles (arc starting with degree 0)
							for (l = 1; l < COUNT; l++)
								deg[l] = fmod(deg[l] + 360.0 - deg[0], 360.0);

							deg[0] = 0.0;

							// select intersection between projection of k and i (on side of k and closer to projection of k)
							bool firstSide = (deg[3] < deg[2]);
							int ppIndex = -1;	// assume no intersections inside arc

							for (l = 0; l < 2; l++)
								if ((deg[4 + l] <= deg[1]) && (!firstSide ^ (deg[4 + l] < deg[2])))
								{
									if (ppIndex == -1)
										ppIndex = l;
									else
										ppIndex = (abs(deg[4] - deg[3]) < abs(deg[5] - deg[3])) ? 0 : 1;
								}

							if (ppIndex != -1)
							{
								// add to centroid
								projPoints[k] = projPoints[k] + pp[ppIndex];
								weights[k]++;
							}
						}
					}
				}
			}
		}

	if (DEBUG)
		cout << endl;

	// weigh to get centroid
	for (i = 0; i < (int)points.size(); i++)
	{
		if (weights[i] > 0)
			projPoints[i] = projPoints[i]*(1.0/weights[i]);
		else
			projPoints[i] = points[i];	// keep original value (also outliers, nonfitted)
	}
}

/*
 * denoise points: minimize energy function that removes noise from features
 */
void denoisePointsLinearFunction(vector<Point> &points, vector<PState> &pStates,
		vector<vector<int> > *nhood, vector<int[2]> &neighbors,
		vector<Point> &normals, vector<Circle> &circles, vector<Point> &denoisedPoints)
{
	const int FACTOR = 1;
	int i, j, k, iter;
	denoisedPoints = points;
	vector<Point> newPoints = points;

	// iteratively move each sample along its normal to minimize the energy functional
	for (iter = 0; iter < 10; iter++)
	{
		// DEBUG
		cout << "iter #" << iter << ":";

		for (i = 0; i < (int)points.size(); i++)
		{
			if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
			{
				// DEBUG
				if (i == 40)
					cout << "";

				int tIndex[2], t[2] = { neighbors[i][0], neighbors[i][1] };
				Point currP = denoisedPoints[i];

				// to locate minimum, evaluate function values only at each neighbor |f(n)|=0 since linear in-between
				vector<int> nVec;
				nVec.push_back(i);

				for (j = 0; j < 2; j++)
					nVec.insert(++nVec.begin(), nhood[j][i].begin(), nhood[j][i].end());

				for (j = 0; j < (int)nVec.size(); j++)
					for (k = 0; k < 2; k++)
						if (nVec[j] == t[k])
							tIndex[k] = j;

				float yd[nVec.size()];

				for (j = 0; j < (int)nVec.size(); j++)
				{
					Point p;

					if ((j == 0) || (j == tIndex[0]) || (j == tIndex[1]))
						p = denoisedPoints[nVec[j]];
					else
						p = points[nVec[j]];

					yd[j] = normals[i].dot(p - currP);
				}

				// assume tangent line through each point in nhood
				float fMin = numeric_limits<float>::max();
				float ydMin = 0.0;

				for (j = 0; j < (int)nVec.size(); j++)
				{
					// evaluate F=a*(|t0,s'|+|t1,s'|+2*|s,s'|)+sum_N(s)|n(s),s'|; s' is variable and s, t[0|1] change iteratively
					float dt[2], ds = abs(yd[0] - yd[j]);

					for (k = 0; k < 2; k++)
						dt[k] = abs(yd[tIndex[k]] - yd[j]);

					float f = FACTOR*(2*ds + dt[0] + dt[1]);

					for (k = 0; k < (int)nVec.size(); k++)
						f += abs(yd[k] - yd[j]);

					if (f < fMin)
					{
						fMin = f;
						ydMin = yd[j];
					}
				}

				// DEBUG
				if (abs(ydMin) > 0.001)
				{
					cout << " " << i;

					if (iter > 10)
						cout << "";
				}

				newPoints[i] = denoisedPoints[i] + normals[i]*(ydMin - yd[0]);

				// TODO: test quadratic function (allow movement inside circular disc, not restricted to L1)
				// TODO: test circular instead of linear fit?

				// TODO: clamp at limits inside PDF, avoid self-intersections
			}
		}

		// DEBUG
		cout << endl;

		denoisedPoints = newPoints;
	}
}

/*
 * denoise points: minimize energy function that removes noise from features
 */
void denoisePointsCentroid(vector<Point> &points, vector<PState> &pStates,
		vector<vector<int> > *nhood, vector<int[2]> &neighbors,
		vector<Point> &normals, vector<Circle> &circles, vector<Point> &denoisedPoints)
{
	int i, j, k, iter;
	denoisedPoints = points;
	vector<Point> newPoints = points;

	// iteratively move each sample along its normal to minimize the energy functional
	for (iter = 0; iter < 1; iter++)
	{
		// DEBUG
		cout << "iter #" << iter << ":";

		for (i = 0; i < (int)points.size(); i++)
		{
			if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
			{
				// DEBUG
				if (i == 73)
					cout << "";

				int tIndex[2], t[2] = { neighbors[i][0], neighbors[i][1] };
				Point currP = denoisedPoints[i];

				// to locate minimum, evaluate function values only at each neighbor |f(n)|=0 since linear in-between
				vector<int> nVec;
				nVec.push_back(i);

				for (j = 0; j < 2; j++)
					nVec.insert(++nVec.begin(), nhood[j][i].begin(), nhood[j][i].end());

				for (j = 0; j < (int)nVec.size(); j++)
					for (k = 0; k < 2; k++)
						if (nVec[j] == t[k])
							tIndex[k] = j;

				nVec.push_back(i);	// double current point to avoid shrinking (of conforming points)

				Point centroid(0.0, 0.0);

				for (j = 0; j < (int)nVec.size(); j++)
				{
					Point p;

					if ((j == 0) || (j == tIndex[0]) || (j == tIndex[1]))
						p = denoisedPoints[nVec[j]];
					else
						p = points[nVec[j]];

					centroid = centroid + p;
				}

				centroid = centroid*(1.0/nVec.size());

#ifdef NOTGOOD
				// TEST: new function
				Point side[2];

				for (j = 0; j < 2; j++)
				{
					side[j] = Point(0.0, 0.0);

					for (auto n:nhood[j][i])
						side[j] = side[j] + points[n];

					side[j] = side[j]*(1.0/nhood[j][i].size());
				}

				centroid = currP*0.5 + side[0]*0.25 + side[1]*0.25;
#endif

				// DEBUG
				if (currP.distance(centroid) > 0.001)
					cout << " " << i;

				newPoints[i] = centroid;

				// TODO: test circular instead of linear fit?

				// TODO: clamp at limits inside PDF, avoid self-intersections
			}
		}

		// DEBUG
		cout << endl;

		denoisedPoints = newPoints;
	}
}

#ifdef NEW
/*
 * construct line touching the top or bottom of both discs (p0, r0) and (p1, r1)
 */
bool constructLine(Point p0, Point p1, float r0, float r1, Point &tp, Point &tv, bool isTop0, bool isTop1)
{
	float centerDist = p0.distance(p1);

	// ensure r0 >= r1
	if (r0 < r1)
	{
		swap(p0, p1);
		swap(r0, r1);
	}

	if (centerDist <= r1)
		return false;

	if ((r0 == r1) && (isTop0 == isTop1))
	{
		// exterior tangent line is parallel to p0-p1
		tv = p1 - p0;
		tv.normalize();
		Point n(tv[1], -tv[0]);
		n = n*r0;

		if (isTop0 ^ (n[1] < 0.0))
			tp = p0 + n;
		else
			tp = p0 - n;
	}
	else
	{
		if (isTop0 == isTop1)
			tp = Point((p1[0]*r0 - p0[0]*r1)/(r0 - r1), (p1[1]*r0 - p0[1]*r1)/(r0 - r1));	// exterior tangent
		else
			tp = Point((p1[0]*r0 + p0[0]*r1)/(r0 + r1), (p1[1]*r0 + p0[1]*r1)/(r0 + r1));	// interior tangent

		float dx = tp[0] - p0[0];
		float dy = tp[1] - p0[1];
		float denom = SQR(dx) + SQR(dy);
		float root = sqrt((denom) - SQR(r0));
		Point tv0((SQR(r0)*dx + r0*dy*root)/denom, (SQR(r0)*dy - r0*dx*root)/denom);
		Point tv1((SQR(r0)*dx - r0*dy*root)/denom, (SQR(r0)*dy + r0*dx*root)/denom);
		tv = (isTop0 ^ (tv0[1] < tv1[1])) ? tv0 : tv1;
	}

	return true;
}

/*
 * returns distance of point p to line(p0, v) with normalized vector v
 */
float distancePointToLine(Point p, Point p0, Point v)
{
	// ASSERT
	Point v2 = v;
	v2.normalize();
	assert((v[0] == v2[0]) && (v[1] == v2[1]));

	Point n(v[1], -v[0]);

	// make normal point upwards
	if (n[1] < 0.0)
	{
		n[0] = -n[0];
		n[1] = -n[1];
	}

	return n.dot(p - p0);
}
#endif

/*
 * denoise points: minimize energy function that removes noise from features
 */
void denoisePointsL0L2(vector<Point> &points, vector<PState> &pStates,
		vector<vector<int> > *nhood, vector<int[2]> &neighbors,
		vector<Point> &normals, vector<Circle> &circles, vector<float> &noiseIn, vector<Point> &denoisedPoints)
{
#ifdef NEW
	int i, j;
	vector<float> noise, noiseDe;	// noise extent determined by connectivity reconstruction
	vector<bool> straight;
	noiseDe.resize(points.size());
	noise.resize(points.size());
	straight.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
		noiseDe[i] = 0.0;

	// determine noise extent at points as largest distance of samples from any referencing neighborhood
	for (i = 0; i < (int)points.size(); i++)
	{
		for (j = 0; j < 2; j++)
			for (auto n:nhood[j][i])
			{
				if (circles[n].variance > noiseDe[i])
					noiseDe[i] = circles[n].variance;
			}
	}

	// determine noise extent as max of modeled input and connectivity-determined noise
	for (i = 0; i < (int)points.size(); i++)
		noise[i] = (noiseIn[i] > noiseDe[i]) ? noiseIn[i] : noiseDe[i];

	// determine whether line fits neighborhood within noise extent
	for (i = 0; i < (int)points.size(); i++)
	{
		vector<int> samples;
		samples.insert(samples.end(), nhood[0][i].rend(), nhood[0][i].rbegin());
		samples.push_back(i);
		samples.insert(samples.end(), nhood[1][i].begin(), nhood[1][i].end());

		// construct constraining lines for both top and bottom of line bundle passing through discs
		for (j = 0; j < 2; j++)
		{
			// adjust line from bottom to pass through discs formed by samples + radius of noise extent
			vector<int>::iterator sIter = samples.begin();
			sIter++;
			sIter++;

			// construct starting line with first two samples
			int firstP = samples[0], lastP = samples[1];
			Point p, v;
			bool isTop = (j == 0);

			if (!constructLine(points[firstP], points[lastP], noise[firstP], noise[lastP], p, v, isTop, isTop))
			{
				// circles are contained within each other - no single tangent exists, choose horizontal tangent to disc with smaller radius
				bool isFirstSmaller = (noise[firstP] < noise[lastP]);
				p = points[firstP];
				v = points[lastP] - points[firstP];
				Point n = Point(v[1], -v[0]);

				if (isTop ^ (n[1] < 0.0))
				{
					n[0] = -n[0];
					n[1] = -n[1];
				}

				if (isFirstSmaller)
				{
					p = p + n*noise[firstP];
					firstP = lastP;
				}
				else
				{
					p = p + n*noise[lastP];
					lastP = firstP;
				}
			}
			else	// exterior tangent to two circles
			{
				straight[i] = true;

				// adjust for subsequent samples
				while ((sIter != samples.end()) && straight[i])
				{
					int currP = *sIter;

					// determine signed distance from sample to line in upward direction
					float dist = distancePointToLine(points[currP], p, v);

					// test whether next sample is inside
					if (abs(dist) > noise[currP])
					{
						// adjust line accordingly
						if (isTop ^ (dist < 0.0))
						{
							// is below disc: ext-tangent from first limit to current sample and set last=current
							straight[i] = (constructLine(points[firstP], points[currP], noise[firstP], noise[currP], p, v, isTop, isTop));
							lastP = currP;
						}
						else
						{
							// is above disc: in-tangent from last limit to current sample and set first=last
							straight[i] = (constructLine(points[firstP], points[currP], noise[firstP], noise[currP], p, v, isTop, !isTop));
							firstP = lastP;
						}

						// check all other samples if inside disc
						vector<int>::iterator s2Iter = samples.begin();

						while ((s2Iter != sIter) && straight[i])
							straight[i] = (distancePointToLine(points[*s2Iter], p, v) <= noise[*s2Iter]);
					}
				}
			}
		}

		// draw lines
		// TODO
	}
#endif
}

/*
 * cross product of 2D vectors
 */
float crossVec(Point a, Point b)
{
	return a[0]*b[1] - a[1]*b[0];
}

/*
 * tests if edge p0-p1 intersects edge q0-q1
 */
bool intersectsEdge(Point p0, Point p1, Point q0, Point q1, Point x)
{
	Point r = p1 - p0;
	Point s = q1 - q0;
	float crossRS = crossVec(r, s);
	float t = crossVec(q0 - p0, s)/crossRS;
	float u = crossVec(q0 - p0, r)/crossRS;
	x = p0 + r*(t/sqrt(r.squared_length()));

	return ((crossRS == 0.0) || ((t >= 0.0) && (t <= 1.0) && (u >= 0.0) && (u <= 1.0)));
}

/*
 * return k-nearest neighbor of currP
 */
int getKNeighbor(ANNkd_tree *kdTree, ANNpointArray ann_points, int currP, int pos)
{
	// TODO: cache previously computed neighbors

	ANNidxArray nnIdx = new ANNidx[pos];
	ANNdistArray distances = new ANNdist[pos];
	kdTree->annkSearch(ann_points[currP], pos, nnIdx, distances);

	return nnIdx[pos - 1];
}

/*
 * fit circle and order points according to normal
 */
bool fitCircleAndOrderNeighbors(vector<Point> &points, vector<int> &indices, vector<int> *orderedNeighbors, Circle &circle)
{
	int j, k = indices.size();
	map<int, Point> ppoints;
	Data data(k);

	for (j = 0; j < k; j++)
	{
		Point p = points[indices[j]];
		data.X[j] = p[0];
		data.Y[j] = p[1];
	}

	// TODO
//	if (indices.size() == 3)
//		circle = CircleFit3Points(data);	// degenerate case at least for 3 points on vertical line
//	else
		circle = CircleFitByHyper(data);

	Point center(circle.a, circle.b);

	// project neighbor points and order them (left, right separately)
	multimap<float, int> angleMap;
	Point centerNormal(points[indices[0]] - center);
	centerNormal.normalize();
	circle.variance = 0.0;

	for (j = 0; j < k; j++)
	{
		Point p = points[indices[j]];
		Point normal = (p - center);
		normal.normalize();
		ppoints[indices[j]] = center + normal*circle.r;

		float angle = atan2(normal.y(), normal.x()) - atan2(centerNormal.y(), centerNormal.x());

		if (angle < -PI)
			angle += 2*PI;
		else
		if (angle > PI)
			angle -= 2*PI;

		angleMap.insert(pair<float, int>(angle, indices[j]));

		float var = abs(p.distance(center) - circle.r);

		if (var > circle.variance)
			circle.variance = var;
	}

	int sideIndex = 0;

	for (auto entry:angleMap)
	{
		if (entry.second == indices[0])
			sideIndex++;
		else
			orderedNeighbors[sideIndex].push_back(entry.second);
	}

	reverse(orderedNeighbors[0].begin(), orderedNeighbors[0].end());

	// test if feature condition fulfilled
	bool featureFitted = false;

	if ((orderedNeighbors[0].size() > 0) && (orderedNeighbors[1].size() > 0))
	{
		Point currPP = ppoints[indices[0]];
		int lastNN[2] = { orderedNeighbors[0].back(), orderedNeighbors[1].back() };
		Point lastPP0 = ppoints[lastNN[0]];
		Point lastPP1 = ppoints[lastNN[1]];
		float nnDist = lastPP0.distance(lastPP1);
		featureFitted = ((nnDist > currPP.distance(lastPP0)) && (nnDist > currPP.distance(lastPP1)));
	}

	return featureFitted;
}

/*
 * returns whether neighbor is CW-oriented on circle fitted to point index
 */
bool isFitCW(vector<Point> &points, vector<int[2]> &neighbors,
		vector<Circle> &circles, int index)
{
	int i;
	Point center = Point(circles[index].a, circles[index].b);
	Point n[2] = { points[neighbors[index][0]] - center,
					points[index] - center };
	float angle[2];

	for (i = 0; i < 2; i++)
		angle[i] = (computeDegreesFromVector(n[i])+PI)/PI*180;

	float arcAngle = fmod(angle[1] - angle[0] + 360.0, 360.0);

	return (arcAngle > 180);
}

/*
 * orient manifold points, returns true if single closed manifold
 */
bool orientManifoldPoints(vector<Point> &points, vector<PState> &pStates, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, vector<Point> &normals, vector<Circle> &circles)
{
	int i, j;
	bool isAllClosed = true;

	// orient neighbors by their projections and propagating along neighbors, using a list for O(N) time complexity, assuming single connected component
	vector<bool> isOriented(points.size());

	for (i = 0; i < (int)points.size(); i++)
		isOriented[i] = false;

	stack<int> unorientedStack;
	i = 0;
	int ccCount = 0;

	// process all points in list
	while (i < (int)points.size())
	{
		// seek to a yet unoriented point which could be fitted (to process all connected components)
		while ((isOriented[i] || ((pStates[i] != PState::MANIFOLD) && (pStates[i] != PState::SHARP))) && (i < (int)points.size()))
			i++;

		bool isStartCW;

		if (i < (int)points.size())
		{
			unorientedStack.push(i);
			isStartCW = isFitCW(points, neighbors, circles, i);
			ccCount++;

			if (DEBUG)
				cout << "starting connected component #" << ccCount << ":" << isStartCW << endl;
		}

		list<int> cc;

		// process all points in stack
		while (unorientedStack.size() > 0)
		{
			int currP = unorientedStack.top();
			unorientedStack.pop();
			isOriented[currP] = true;
			cc.push_back(currP);

			if (DEBUG)
				cout << " " << currP;

			for (j = 0; j < 2; j++)
			{
				int n = neighbors[currP][j];

				if ((n != -1) && (!isOriented[n]))
				{
					if (neighbors[n][1 - j] != currP)
					{
						swap(neighbors[n][0], neighbors[n][1]);
						swap(nhood[0][n], nhood[1][n]);

						if (DEBUG)
							cout << " (flip #" << n << ")";
					}

					if ((pStates[n] == PState::LEAF) || (pStates[n] == PState::SHARP))
					{
						// since point n not fitted, assign neighbor's orientation by dot product of normals
						float dot = normals[currP][0]*normals[n][0] + normals[currP][1]*normals[n][1];

						if (dot < 0.0)
						{
							normals[n][0] = -normals[n][0];
							normals[n][1] = -normals[n][1];

							if (DEBUG)
								cout << "(flip normal for leaf or sharp point)";
						}
					}
					else
					{
						bool isCurrCW = isFitCW(points, neighbors, circles, n);

						if (DEBUG)
							cout << "(" << isCurrCW << ")";

						if (isCurrCW != isStartCW)
						{
							normals[n][0] = -normals[n][0];
							normals[n][1] = -normals[n][1];

							if (DEBUG)
								cout << "(flip normal)";
						}
					}

					isOriented[n] = true;

					if ((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::SHARP))
						unorientedStack.push(n);	// no points get pushed twice on stack since isOriented flag has been set
				}
			}

			if (DEBUG)
				cout << endl;
		}

		// test whether closed
		if (isAllClosed && (cc.size() > 0))
		{
			int lastP = cc.back();

			if ((neighbors[lastP][0] != cc.front()) && (neighbors[lastP][1] != cc.front()))
				isAllClosed = false;
		}

		// DEBUG
		if (cc.size() > 0)
		{
			int lastP = cc.back();

			if (DEBUG)
			{
				if ((neighbors[lastP][0] == cc.front()) || (neighbors[lastP][1] == cc.front()))
					cout << "Closed";
				else
					cout << "Open";

				cout << " CC with " << cc.size() << " points (+2 end points if open)" << endl;
			}
		}
	}

	if (DEBUG)
		cout << ccCount << " connected components oriented consistently" << endl;

	return (isAllClosed && (ccCount == 1));
}

class IntPair
{
public:
	int a, b;

	IntPair(int p_a, int p_b)
	{
		a = p_a;
		b = p_b;

		if (a > b)
			swap(a, b);
	}

	bool operator<(const IntPair& ip) const
	{
		return (a < ip.a) || ((a == ip.a) && (b < ip.b));
	}

	bool operator==(const IntPair& ip) const
	{
		return ((a == ip.a) && (b == ip.b));
	}
};

/*
 * add next neighbor point and return point state
 */
bool increaseNHood(ANNkd_tree *kdTree, ANNpointArray ann_points, vector<Point> &points,
		vector<PState> &pStates, vector<Circle> &circles, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, int currP, Point &normal, map<int, list<int> > &nhoodRefsMap)
{
	int i;
	vector<int> indices(2 + nhood[2][currP].size());

	if (pStates[currP] != PState::INITIAL)
	{
		int pos = 2 + nhood[2][currP].size();

		if (pos > (int)points.size())
			return false;

		int k = getKNeighbor(kdTree, ann_points, currP, pos);
		nhood[2][currP].push_back(k);
		indices[0] = currP;

		for (i = 0; i < (int)nhood[2][currP].size(); i++)
			indices[1 + i] = nhood[2][currP][i];

		nhoodRefsMap[k].push_back(currP);
	}
	else
	{
		indices.resize(3);
		indices[0] = currP;

		// DEBUG
		if (currP == 1)
			cout << "";

		for (i = 0; i < 2; i++)
		{
			int k = getKNeighbor(kdTree, ann_points, currP, 2 + i);
			indices[1 + i] = k;
			nhood[2][currP].push_back(k);
			nhoodRefsMap[k].push_back(currP);
		}
	}

	vector<int> orderedNeighbors[2];
	orderedNeighbors[0].clear();
	orderedNeighbors[1].clear();
	Circle circle;
	bool featureFitted = fitCircleAndOrderNeighbors(points, indices, orderedNeighbors, circle);
	circles[currP] = circle;

	for (i = 0; i < 2; i++)
		nhood[i][currP] = orderedNeighbors[i];

	// compute normal
	Point center(circle.a, circle.b);
	normal = points[currP] - center;
	normal.normalize();

	pStates[currP] = (featureFitted ? PState::CONFORM : PState::NONCONFORM);

	return true;
}

/*
 * determine redundant points in overlap between two conforming points and add to set
 */
void addRedundantPointsInOverlap(int currP, int side, vector<vector<int> > *nhood, vector<int[2]> &neighbors, set<int> &unhandledPointSet, list<int> &redundantPointsInOverlap)
{
	int currN = neighbors[currP][side];
	set<int> set0, set1;

	// add all active points in nhood of currP up to currN to set0
	vector<int>::iterator nIter = nhood[side][currP].begin();

	while ((nIter != nhood[side][currP].end()) && (*nIter != currN))
	{
		if (unhandledPointSet.find(*nIter) != unhandledPointSet.end())
			set0.insert(*nIter);

		nIter++;
	}

	assert(nIter != nhood[side][currP].end());	// currN must be in neighbors

	// add all active points in nhood of currN up to currP to set1
	int side2 = ((neighbors[currN][0] == currP) ? 0 : 1);
	nIter = nhood[side2][currN].begin();

	while ((nIter != nhood[side2][currN].end()) && (*nIter != currP))
	{
		if (unhandledPointSet.find(*nIter) != unhandledPointSet.end())
			set1.insert(*nIter);

		nIter++;
	}

	assert(nIter != nhood[side2][currN].end());	// currP must be in neighbors

	// determine intersection of the two sets as shared points and add to redundant set
	vector<int> redundantVec(set0.size() + set1.size());
	vector<int>::iterator endIter = set_intersection(set0.begin(), set0.end(), set1.begin(), set1.end(), redundantVec.begin());
	redundantPointsInOverlap.insert(redundantPointsInOverlap.begin(), redundantVec.begin(), endIter);
}

/*
 * check rules for overlap of point p with neighbor n
 */
void assertOverlap(int p, int n, int side, vector<PState> &pStates, vector<vector<int> > *nhood, vector<int[2]> &neighbors)
{
	assert((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::CONFORM) ||
			(pStates[n] == PState::LEAF) || (pStates[n] == PState::SHARP));

	// verify that overlap is mutual from n to p
	assert((neighbors[n][0] == p) || (neighbors[n][1] == p));

	if ((pStates[p] == PState::SHARP) || (pStates[n] == PState::SHARP))
		return;

	// follow neighbors of p from n while inside nhood side
	set<int> visitedSet, invSet;
	int prevP = p;
	int currP = n;
	int nextP = neighbors[currP][(neighbors[currP][0] == prevP) ? 1 : 0];
	visitedSet.insert(currP);

	// follow and mark neighbor outwards sequentially inside nhood as visited
	while (find(nhood[side][p].begin(), nhood[side][p].end(), nextP) != nhood[side][p].end())
	{
		prevP = currP;
		currP = nextP;
		visitedSet.insert(currP);
		nextP = neighbors[currP][(neighbors[currP][0] == prevP) ? 1 : 0];
	}

	// assert that only unhandled points remain in nhood side (CONFORM, NONCONFORM, REDUNDANT or OUTLIER)
	for (auto n2:nhood[side][p])
		if ((pStates[n2] != PState::CONFORM) && (pStates[n2] != PState::NONCONFORM) &&
			(pStates[n2] != PState::ELIMINATED) && (pStates[n2] != PState::OUTLIER) &&
			(visitedSet.find(n2) == visitedSet.end()))
			assert(false);
}

/*
 * check rules for point
 */
void checkRuleForPoint(int p, vector<PState> &pStates, vector<vector<int> > *nhood, vector<int[2]> &neighbors)
{
	int i;

	if (pStates[p] == PState::ELIMINATED)
	{
		for (i = 0; i < 2; i++)
			assert(neighbors[p][i] == -1);
	}
	else
	if (pStates[p] == PState::NONCONFORM)
	{
		for (i = 0; i < 2; i++)
			assert(neighbors[p][i] == -1);
	}
	else
	if (pStates[p] == PState::SHARP)
	{
		for (i = 0; i < 2; i++)
		{
			int n = neighbors[p][i];
			assert(n != -1);
			assert((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::CONFORM));
			assertOverlap(p, n, i, pStates, nhood, neighbors);
		}
	}
	else
	if (pStates[p] == PState::LEAF)
	{
		int n = -1, side = -1;

		if (neighbors[p][0] != -1)
		{
			n = neighbors[p][0];
			side = 0;
		}

		if (neighbors[p][1] != -1)
		{
			assert(n == -1);
			n = neighbors[p][1];
			side = 1;
		}

		assert(n != -1);
		assertOverlap(p, n, side, pStates, nhood, neighbors);	// can contain other MANIFOLD points in nhood, from intersecting closed curve
	}
	else
	if (pStates[p] == PState::CONFORM)
	{
		for (i = 0; i < 2; i++)
		{
			int n = neighbors[p][i];

			if ((n != -1) && (pStates[n] != PState::SHARP))
				assertOverlap(p, n, i, pStates, nhood, neighbors);
		}
	}
	else
	if (pStates[p] == PState::MANIFOLD)
	{
		for (i = 0; i < 2; i++)
		{
			int n = neighbors[p][i];
			assert(n != -1);

			if (pStates[n] != PState::SHARP)
				assertOverlap(p, n, i, pStates, nhood, neighbors);
		}
	}
}

/*
 * DEBUG output
 */
void outputStateStatistics(int iterations, set<int> &unhandledPointSet, vector<Point> &points, vector<PState> &pStates)
{
	int i;

	cout << "ITER #" << iterations << ": " << unhandledPointSet.size() << " unhandled points:";

	for (auto n:unhandledPointSet)
		cout << " " << n;

	cout << ", SHARP:";

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] == PState::SHARP)
			cout << " " << i;

	cout << ", LEAF:";

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] == PState::LEAF)
			cout << " " << i;

	cout << ", MANIFOLD:";

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] == PState::MANIFOLD)
			cout << " " << i;

	cout << ", ELIMINATED:";

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] == PState::ELIMINATED)
			cout << " " << i;

	cout << endl;
}

/*
 * DEBUG output neighborhood
 */
void outputNeighborhood(int currP, vector<vector<int> > *nhood)
{
	int i;

	cout << "nhood:";

	for (i = 0; i < 2; i++)
	{
		for (auto p:nhood[i][currP])
			cout << p << " ";

		cout << "; ";
	}
}

/*
 * determines overlap for point's neighborhood, per side
 */
int determinePointOverlapSide(int currP, int side, vector<vector<int> > *nhood, vector<int[2]> &neighbors,
		vector<PState> &pStates)
{
	int j;
	bool overlap = false;

	// check for complete overlap between available neighbors: n(n(p)_1)_0 = p
	vector<int>::iterator nIter = nhood[side][currP].end();

	while (!overlap && (nIter != nhood[side][currP].begin()))
	{
		nIter--;
		int currN = *nIter;

		if ((pStates[currN] == PState::CONFORM) || (pStates[currN] == PState::LEAF) || (pStates[currN] == PState::MANIFOLD))
		{
			j = 0;

			while (!overlap && (j < 2))
			{
				for (auto currN2:nhood[j][currN])
					// if neighbor's nhood includes currP: potential overlap
					if (currN2 == currP)
					{
						// if neighbor is already manifold
						if (pStates[currN] == PState::MANIFOLD)
						{
							// then one of its neighbors must be currP to be an overlap
							if ((neighbors[currN][0] == currP) || (neighbors[currN][1] == currP))
								overlap = true;
						}
						else	// if neighbor conforming, instead of manifold
						{
							// neighbor's neighbor must be either empty or currP, to avoid creating dangling pointers to neighbors
							if ((neighbors[currN][j] == -1) || (neighbors[currN][j] == currP))
								overlap = true;
						}

						if (overlap)
							return currN;
					}

				j++;
			}
		}
	}

	return -1;
}

/*
 * change point state consistently
 */
void changePointState(int p, PState state, vector<PState> &pStates, set<int>::iterator &pIter, set<int> &unhandledPointSet)
{
	// each state change (except CONFORM <-> NONCONFORM) is between handled and unhandled:
	if ((state == PState::MANIFOLD) || (state == PState::ELIMINATED) || (state == PState::LEAF) || (state == PState::SHARP))
	{
		if (pStates[p] != PState::LEAF)
		{
			// if already handled:
			assert(unhandledPointSet.find(p) != unhandledPointSet.end());

			// advance pointer to remain at position in loop
			if (*pIter == p)
				pIter++;

			unhandledPointSet.erase(p);
		}
	}
	else
	{
		// if yet unhandled:
		assert(unhandledPointSet.find(p) == unhandledPointSet.end());
		unhandledPointSet.insert(p);
	}

	pStates[p] = state;

	if (DEBUG)
	{
		cout << "#" << p << "->";

		switch (state)
		{
			case PState::MANIFOLD: cout << "M"; break;
			case PState::CONFORM: cout << "C"; break;
			case PState::NONCONFORM: cout << "N"; break;
			case PState::ELIMINATED: cout << "E"; break;
			case PState::SHARP: cout << "S"; break;
			case PState::LEAF: cout << "L"; break;
			default: assert(false); break;
		}

		cout << " ";
	}
}

/*
 * remove overlap: do actual work
 */
void doRemoveOverlap(int currP, int side, vector<int[2]> &neighbors, vector<PState> &pStates, set<int>::iterator &pIter, set<int> &unhandledPointSet)
{
	if (neighbors[currP][side] != -1)
	{
		int currN = neighbors[currP][side];
		int nSide = (neighbors[currN][0] == currP) ? 0 : 1;
		neighbors[currP][side] = -1;
		neighbors[currN][nSide] = -1;

		// change state of neighbor
		if ((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::LEAF))
			changePointState(currN, PState::CONFORM, pStates, pIter, unhandledPointSet);
		else
		if (pStates[currN] == PState::SHARP)
		{
			changePointState(currN, PState::NONCONFORM, pStates, pIter, unhandledPointSet);

			// also remove other overlap of SHARP point
			int currN2 = neighbors[currN][1 - nSide];
			int n2Side = (neighbors[currN2][0] == currN) ? 0 : 1;
			neighbors[currN][1 - nSide] = -1;
			neighbors[currN2][n2Side] = -1;

			// neighbor of SHARP can only be MANIFOLD (except when removing inconsistencies)
			if (pStates[currN2] == PState::MANIFOLD)
				changePointState(currN2, PState::CONFORM, pStates, pIter, unhandledPointSet);
		}
	}
}

/*
 * remove one overlap of point
 */
void removeOverlap(int currP, int side, vector<int[2]> &neighbors, vector<PState> &pStates,
		set<int>::iterator &pIter, set<int> &unhandledPointSet, set<int> &affectedPointSet)
{
	assert((pStates[currP] == PState::MANIFOLD) || (pStates[currP] == PState::CONFORM));

	doRemoveOverlap(currP, side, neighbors, pStates, pIter, unhandledPointSet);

	if (neighbors[currP][side] != -1)
		affectedPointSet.insert(neighbors[currP][side]);

	affectedPointSet.insert(currP);

	// change state of current point
	if (pStates[currP] == PState::MANIFOLD)
		changePointState(currP, PState::CONFORM, pStates, pIter, unhandledPointSet);
}

/*
 * remove overlaps of point
 */
void removeOverlaps(int currP, vector<int[2]> &neighbors, vector<PState> &pStates,
		set<int>::iterator &pIter, set<int> &unhandledPointSet, set<int> &affectedPointSet)
{
	assert((pStates[currP] == PState::MANIFOLD) || (pStates[currP] == PState::CONFORM) ||
			(pStates[currP] == PState::LEAF) || (pStates[currP] == PState::SHARP));

	for (int side = 0; side < 2; side++)
	{
		doRemoveOverlap(currP, side, neighbors, pStates, pIter, unhandledPointSet);

		if (neighbors[currP][side] != -1)
			affectedPointSet.insert(neighbors[currP][side]);
	}

	affectedPointSet.insert(currP);

	// change state of current point
	if ((pStates[currP] == PState::MANIFOLD) || (pStates[currP] == PState::LEAF))
		changePointState(currP, PState::CONFORM, pStates, pIter, unhandledPointSet);
	else
	if (pStates[currP] == PState::SHARP)
		changePointState(currP, PState::NONCONFORM, pStates, pIter, unhandledPointSet);
}

/*
 * determine points contained in nhood side consistent with an overlap (followable)
 */
void determineConsistentPointSetInNhood(int currP, int side, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, vector<PState> &pStates, set<int> &consistentPointSet)
{
	int p = currP;
	int prevP = -1, nextP = neighbors[p][side];

	// follow neighbor outwards sequentially inside nhood
	while (find(nhood[side][currP].begin(), nhood[side][currP].end(), nextP) != nhood[side][currP].end())
	{
		prevP = p;
		p = nextP;
		consistentPointSet.insert(p);
		nextP = neighbors[p][(neighbors[p][0] == prevP) ? 1 : 0];
	}
}

/*
 * remove inconsistent points contained in nhood side of an overlap
 */
void removeInconsistentPointSetInNhood(int currP, int side, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, vector<PState> &pStates,
		set<int>::iterator &pIter, set<int> &unhandledPointSet, set<int> &affectedPointSet)
{
	// determine points followed from subsequent overlaps in nhood
	set<int> consistentPointSet;

	determineConsistentPointSetInNhood(currP, side, nhood, neighbors, pStates, consistentPointSet);

	// for all remaining handled points (MANIFOLD, SHARP or LEAF)
	for (auto currN:nhood[side][currP])
	{
		if ((consistentPointSet.find(currN) == consistentPointSet.end()) &&
			((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::LEAF) ||
			(pStates[currN] == PState::SHARP)))
		{
			// remove overlaps
			removeOverlaps(currN, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);

			// mark point as unhandled
			if (pStates[currN] == PState::SHARP)
				changePointState(currN, PState::NONCONFORM, pStates, pIter, unhandledPointSet);
			else
			if ((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::LEAF))
				changePointState(currN, PState::CONFORM, pStates, pIter, unhandledPointSet);
		}
	}
}

/*
 * determine if point is inside an overlap (also with SHARP point)
 */
bool isInsideEdge(int currP, vector<int[2]> &neighbors, vector<vector<int> > *nhood, vector<PState> &pStates)
{
	// if not SHARP or MANIFOLD (or neighbor to one of these)
	if (((pStates[currP] != PState::MANIFOLD) || (pStates[currP] != PState::SHARP)) &&
		((neighbors[currP][0] == -1) || ((pStates[neighbors[currP][0]] != PState::MANIFOLD))) &&
		((neighbors[currP][1] == -1) || ((pStates[neighbors[currP][1]] != PState::MANIFOLD))))
	{
		// test if inside an overlap (also with SHARP point) = a manifold edge of boundary
		list<int> manifoldPoints[2];

		for (int side = 0; side < 2; side++)
		{
			for (auto n:nhood[side][currP])
				if ((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::SHARP))
					manifoldPoints[side].push_back(n);
		}

		for (auto n:manifoldPoints[0])
		{
			int n2 = neighbors[n][0];

			if (find(manifoldPoints[1].begin(), manifoldPoints[1].end(), n2) == manifoldPoints[1].end())
			{
				n2 = neighbors[n][1];

				if (find(manifoldPoints[1].begin(), manifoldPoints[1].end(), n2) == manifoldPoints[1].end())
					n2 = -1;
			}

			if (n2 != -1)
			{
				// verify that no handled points in between in ordering of nhood
				vector<int>::iterator nIter = find(nhood[0][currP].begin(), nhood[0][currP].end(), n);
				vector<int>::iterator iter = nhood[0][currP].begin();
				bool isHandled = false;

				while (!isHandled && (iter != nIter))
				{
					PState state = pStates[*iter];
					isHandled = (state == PState::MANIFOLD) || (state == PState::LEAF) ||
							(state == PState::SHARP) || (state == PState::CONFORM) || (state == PState::NONCONFORM);
					iter++;
				}

				if (!isHandled)
				{
					iter = nhood[1][currP].begin();
					vector<int>::iterator n2Iter = find(nhood[1][currP].begin(), nhood[1][currP].end(), n2);

					while (!isHandled && (iter != n2Iter))
					{
						PState state = pStates[*iter];
						isHandled = (state == PState::MANIFOLD) || (state == PState::LEAF) ||
								(state == PState::SHARP) || (state == PState::CONFORM) || (state == PState::NONCONFORM);
						iter++;
					}

					if (!isHandled)
						return true;
				}
			}
		}
	}

	return false;
}

/*
 * returns distance of point p to line(p0, v) with normalized vector v
 */
float distancePointToLine(Point p0, Point p, Point v)
{
	v.normalize();
	Point n(v[1], -v[0]);

	// make normal point downwards
	if (n[1] > 0.0)
	{
		n[0] = -n[0];
		n[1] = -n[1];
	}

	return n.dot(p - p0);
}

/*
 * calculate tangent point x for line through point p tangent to circle (c, r), upper tangent if isTop==true, else lower tangent, return whether tangent exists
 */
bool tangentPToDisc(Point p, Point c, float r, bool isLeft, bool isTop, Point &x)
{
	double det = SQR(p[0] - c[0]) + SQR(p[1] - c[1]);

	if (det < 0.0)
		return false;

	double d = sqrt(det);
	double det2 = SQR(d) - SQR(r);

	if (det2 < 0.0)
		return false;

	double r1 = sqrt(det2);
	double a = (SQR(r) - SQR(r1) + SQR(d))/(2*d);
	double det3 = SQR(r) - SQR(a);

	if (det3 < 0.0)
		return false;

	double h = sqrt(det3);
	Point v = p - c;
	Point n(v[1], -v[0]);
	n = n*(h/d);
	x = c + v*(a/d);
/*
	if (isTop ^ (n[1] > 0.0))
		x = x + n;
	else
		x = x - n;
*/
	// p is in halfspace of c-n depending on its side in nhood
	if (isTop ^ isLeft)
		x = x + n;
	else
		x = x - n;

	return true;
}

/*
 * intersects lines p0+s*v0 and p1+t*v1 and returns false if parallel, else intersection point x
 */
bool intersectLines(Point p0, Point v0, Point p1, Point v1, Point &x)
{
	double det = (double)v0[1]*v1[0] - (double)v0[0]*v1[1];

	if (det == 0.0)
		return false;

	double s0 = (p0[0]*((double)p0[1] + v0[1]) - p0[1]*((double)p0[0] + v0[0]));
	double s1 = (p1[0]*((double)p1[1] + v1[1]) - p1[1]*((double)p1[0] + v1[0]));
	x[0] = ((double)s0*v1[0] - s1*v0[0])/det;
	x[1] = ((double)s0*v1[1] - s1*v0[1])/det;

	return true;
}

/*
 * intersect circle with line
 */
bool intersectCircleLine2(Circle &circle, Point &p, Point &vec, Point &result0, Point &result1)
{
	float a = -vec[1];
	float b = vec[0];
	float c = a*(p[0] - circle.a) + b*(p[1] - circle.b);	// translate point into coordinate system centered at circle center
	float d2 = SQR(a) + SQR(b);
	float det = SQR(circle.r)*d2 - SQR(c);

	if (det < 0.0)
		return false;

	float temp = b*sqrt(det);
	result0[0] = circle.a + (a*c + temp)/d2;
	result1[0] = circle.a + (a*c - temp)/d2;
	temp = a*sqrt(det);
	result0[1] = circle.b + (b*c - temp)/d2;
	result1[1] = circle.b + (b*c + temp)/d2;

	return true;
}

/*
 * compares points to their y coordinate
 */
bool ySort(Point p0, Point p1)
{
	return (p0[1] < p1[1]);
}

/*
 * reposition vertex along its normal s.t. it minimizes least squares distances from samples and angles of it and its neighbors, constrained to samples discs
 */
Point minimizeVertex(vector<Point> &points, vector<Point> &ppoints, vector<vector<int> > *nhood,
		vector<int[2]> &neighbors, vector<Point> &normals, vector<float> &noise, int currP)
{
	if (noise[currP] == 0.0)
		return points[currP];

	const float EPSILON = 0.0001;
	int i, j, sampleCount[2];
	vector<int> samplesSide[2];
	Point xx[2];

	for (i = 0; i < 2; i++)
	{
		samplesSide[i].push_back(currP);
		samplesSide[i].insert(samplesSide[i].end(), nhood[i][currP].begin(), nhood[i][currP].end());
		sampleCount[i] = samplesSide[i].size();
	}

	// calculate range along normal for which edges and their intersection are inside discs:
	int prev = neighbors[currP][0];
	int next = neighbors[currP][1];
	Point p[2] = { ppoints[prev], ppoints[next] }, x[2];
	Point v(normals[currP][1], -normals[currP][0]);
	bool isFirstLeft = (v.dot(p[0] - points[currP]) > 0.0);	// is left to normal direction

	for (int k = 0; k < 2; k++)	// top/bottom range tangent
	{
		bool isTop = (k == 1);
		bool hasTangent[2] = { false, false };

		for (i = 0; i < 2; i++)	// per side
		{
			bool isLeft = (isFirstLeft ^ (i == 1));

			// from center point to edge point per side, compute highest/lowest tangent to current point if outside its disc
			j = sampleCount[i] - 1;

			while ((j >= 0) && (samplesSide[i][j] != neighbors[currP][i]))
				j--;

			assert(j > 0);

			j--;

			while (j >= 1)	// do not consider tangent to center point
			{
				int curr2P = samplesSide[i][j];

				if (!hasTangent[i] || (abs(distancePointToLine(points[curr2P], p[i], x[i] - p[i])) > noise[curr2P] + EPSILON))
				{
					if (tangentPToDisc(p[i], points[curr2P], noise[curr2P], isLeft, isTop, x[i]))
						hasTangent[i] = true;
				}

				j--;
			}
		}

		// finally, locate highest/lowest center point along normal inside disc and tangents, range is between the intersection points (in terms of height/base ratio)
		vector<Point> pp;
		Point cx[2], tx;
		Circle circle(points[currP][0], points[currP][1], noise[currP]);
		intersectCircleLine2(circle, points[currP], normals[currP], cx[0], cx[1]);

		// select intersection with circle with y value towards desired direction
		if (isTop ^ (normals[currP][1] < 0.0) ^ (cx[0][1] < cx[1][1]))
			pp.push_back(cx[0]);
		else
			pp.push_back(cx[1]);

		for (i = 0; i < 2; i++)	// per side
		{
			bool isLeft = (isFirstLeft ^ (i == 1));

			// only relevant if tangent exists (else all samples on side inside disc)
			if (hasTangent[i])
			{
				intersectLines(p[i], x[i] - p[i], points[currP], normals[currP], tx);

				// only consider if in correct halfspace (top/down)
				Point v = points[currP] - p[i];
				Point n(v[1], -v[0]);

				if (isTop ^ isLeft ^ (n.dot(tx - p[i]) > 0.0))
					pp.push_back(tx);
			}
		}

		// sort points to y-values and select limit value
		sort(pp.begin(), pp.end(), ySort);
		xx[k] = (isTop ^ (normals[currP][1] < 0.0)) ? pp.front() : pp.back();
	}

#ifdef OLD2
	// minimize squared distance to samples together with angle (by its weighted distance to straight line)
	double sumDW = 0.0, sumWW = 0.0;
	Point n = xx[1] - xx[0];
	double length = sqrt(n.squared_length());
	n = n*(1/length);	// normalize normal

	// add distance and weights for samples
	int sampleWeight = 1;	// start with center point, count samples without neighbors (equalizes zigzag line)

	for (i = 0; i < 2; i++)	// per side
	{
//		for (j = 1; j < sampleCount[i]; j++)	// not for center point, 0 impact for neighbor points
		j = sampleCount[i] - 1;	// for samples in nhood closer than neighbor

		while ((j >= 0) && (samplesSide[i][j] != neighbors[currP][i]))
			j--;

		assert(j > 0);

		j--;

		while (j >= 1)	// do not consider tangent to center point
		{
			int currS = samplesSide[i][j];
			Point projS;
			intersectLines(p[i], xx[0] - p[i], ppoints[currS], normals[currP], projS);	// point on bottom line
			double dist = n.dot(projS - ppoints[currS]);
			double weight = 1.0 - projS.distance(xx[0])/p[i].distance(xx[0]);

			if (weight < 0.0)
				weight = 0.0;

			sumDW += dist*weight;
			sumWW += SQR(weight);
			sampleWeight++;
			j--;
		}
	}
/*
	// add distance for center point
	double dist = n.dot(xx[0] - points[currP]);
	double weight = 1.0;
	sumDW += dist*weight;
//	sumDW = 0.0; // TEST
	sumWW += SQR(weight);
*/
#endif

#ifdef OLD
	// TEST
	sumDW = 0.0;
	sumWW = 0.0;

	if (sumWW == 0.0)
		sumWW = 1.0;

	// add weighted distance to straight line (angle approximation)
	Point baseP;
	intersectLines(p[0], p[1] - p[0], ppoints[currP], normals[currP], baseP);	// point on bottom line
	double baseLen = p[1].distance(p[0]);
	double dist = n.dot(xx[0] - baseP);
	float sumDist = /*0.5**/dist;
//	sampleWeight = 1.0;

	// add weighted distances to neighbor baselines (angle approximation)
	int prev2 = neighbors[prev][0];

	if (prev2 == currP)
		prev2 = neighbors[prev][1];

	if (prev2 != -1)
	{
		Point prev2P = ppoints[prev2];
		Point vec = ppoints[prev] - prev2P;
		vec.normalize();
		double dot = v.dot(vec);

		if (dot < 0.0)
		{
			// has intersection with zero angle at prevP
			intersectLines(prev2P, vec, ppoints[currP], normals[currP], baseP);
			double baseLenPrev = prev2P.distance(ppoints[currP]);
			dist = n.dot(xx[0] - baseP);
			sumDist += /*0.25**/abs(dot)*dist*baseLen/baseLenPrev;
//			sampleWeight += abs(dot);
		}
	}

	int next2 = neighbors[next][0];

	if (next2 == currP)
		next2 = neighbors[next][1];

	if (next2 != -1)
	{
		Point next2P = ppoints[next2];
		Point vec = ppoints[next] - next2P;
		vec.normalize();
		double dot = v.dot(vec);

		if (dot > 0.0)
		{
			// has intersection with zero angle at nextP
			intersectLines(next2P, vec, ppoints[currP], normals[currP], baseP);
			double baseLenNext = next2P.distance(ppoints[currP]);
			dist = n.dot(xx[0] - baseP);
			sumDist += /*0.25**/abs(dot)*dist*baseLen/baseLenNext;
//			sampleWeight += abs(dot);
		}
	}

//	sampleWeight = sqrt(sumWW);	// add same weight as samples
	sampleWeight = 3.0;	// TEST
	sumDW += sampleWeight*sumDist;
	sumWW += SQR(sampleWeight);

	// compute distance along normal and point
	double dd = -sumDW/sumWW;

	if ((currP == 87) || (currP == 88) || (currP == 89))
		cout << " #" << currP << ":" << dd/length;

	// clamp to [0..length]
	if (dd < 0.0)
		dd = 0.0;
	else
	if (dd > length)
j		dd = length;

	return xx[0] + n*dd;
#endif

#ifdef OLD2
	// DEBUG
	if (currP == 88)
		cout << "";

	// compute prev, curr and next clockwise angles
	int prev2 = neighbors[prev][0];

	if (prev2 == currP)
		prev2 = neighbors[prev][1];

	int next2 = neighbors[next][0];

	if (next2 == currP)
		next2 = neighbors[next][1];

	Point vec0 = ppoints[prev] - ppoints[prev2];
	Point vec1 = ppoints[currP] - ppoints[prev];
	Point vec2 = ppoints[next] - ppoints[currP];
	Point vec3 = ppoints[next2] - ppoints[next];
	vec0.normalize();
	vec1.normalize();
	vec2.normalize();
	vec3.normalize();
	float prevAngle = (1.0 - vec0.dot(vec1))*PI/2;
	float currAngle = (1.0 - vec1.dot(vec2))*PI/2;
	float nextAngle = (1.0 - vec2.dot(vec3))*PI/2;
	Point n0(vec0[1], -vec0[0]);
	Point n1(vec1[1], -vec1[0]);
	Point n2(vec2[1], -vec2[0]);

	if (n0.dot(vec1) > 0.0)
		prevAngle = -prevAngle;

	if (n1.dot(vec2) > 0.0)
		currAngle = -currAngle;

	if (n2.dot(vec3) > 0.0)
		nextAngle = -nextAngle;

	// determine averaged angle
	float avgAngle = (prevAngle + currAngle + nextAngle)/3.0;

	// determine intersection point with normal as new point
	Point newVec0(vec0[0]*cos(avgAngle) - vec0[1]*sin(avgAngle), vec0[0]*sin(avgAngle) + vec0[1]*cos(avgAngle));
	Point newVec1(vec0[0]*cos(2*avgAngle) - vec0[1]*sin(2*avgAngle), vec0[0]*sin(2*avgAngle) + vec0[1]*cos(2*avgAngle));
	Point prev2P = ppoints[prev2];
	Point prevP = ppoints[prev];
	Point baseP, basePrevP;
	intersectLines(prev2P, newVec0, prevP, normals[prev], basePrevP);
	intersectLines(basePrevP, newVec1, ppoints[currP], normals[currP], baseP);

	// limit with constraints
	Point nn = xx[1] - xx[0];
	float len = sqrt(nn.squared_length());
	nn.normalize();
	float dist = nn.dot(baseP - xx[0]);

	if (dist < 0.0)
		baseP = xx[0];
	else
	if (dist > len)
		baseP = xx[1];

	return baseP;
#endif

	// project nhood on baseline, order inside neighbors [0..1]
	multimap<float, pair<int, int> > samplesMMap;
	Point baseV = ppoints[next] - ppoints[prev];
	float baseLen = sqrt(baseV.squared_length());
	Point baseVN = baseV;
	baseVN.normalize();

	for (i = 0; i < 2; i++)	// per side
		for (auto n:nhood[i][currP])
		{
			// compute signed distance from left neighbor with dot product
			float dist = baseVN.dot(points[n] - ppoints[prev])/baseLen;
			samplesMMap.insert(pair<float, pair<int, int> >(dist, pair<int, int>(n, i)));
		}

	float dist = baseVN.dot(points[currP] - ppoints[prev])/baseLen;
	samplesMMap.insert(pair<float, pair<int, int> >(dist, pair<int, int>(currP, 0)));

	// compute voronoi weights per samples
	map<int, float> voronoiW;
	multimap<float, pair<int, int> >::iterator mmapIter = samplesMMap.begin();
	pair<float, pair<int, int> > prevEntry(-numeric_limits<float>::max(), pair<int, int>(-1, -1));

	// compute voronoi boundaries between all entries
	while (mmapIter != samplesMMap.end())
	{
		pair<float, pair<int, int> > currEntry = *mmapIter;
		mmapIter++;

		// weight is halved distances to neighbor samples
		float rBegin, rEnd;

		if (prevEntry.second.first == -1)
			rBegin = -numeric_limits<float>::max();
		else
			rBegin = 0.5*(prevEntry.first + currEntry.first);

		if (mmapIter == samplesMMap.end())
			rEnd = numeric_limits<float>::max();
		else
			rEnd = 0.5*(currEntry.first + mmapIter->first);

		if (((rBegin >= 0.0) && (rBegin <= 1.0)) || ((rEnd >= 0.0) && (rEnd <= 1.0)))
		{
			if (rBegin < 0.0)
				rBegin = 0.0;

			if (rEnd > 1.0)
				rEnd = 1.0;

			voronoiW[currEntry.second.first] = rEnd - rBegin;
		}

		prevEntry = currEntry;
	}

	// compute point along normal s.t. weighted samples distances zero out
	Point nn = xx[1] - xx[0];
	float len = sqrt(nn.squared_length());
	nn.normalize();
	double sumDV = 0.0, sumV = 0.0, sumW = 0.0, sumWV = 0.0;

	// TEST
	if (currP == 89)
		cout << "";

	for (auto entry:samplesMMap)
	{
		// compute weighted distances + change and weights
		int currS = entry.second.first;
		int side = entry.second.second;
		Point projS;
		intersectLines(p[side], ppoints[currP] - p[side], ppoints[currS], normals[currP], projS);	// point on edge
		double dist = nn.dot(projS - ppoints[currS]);
		double weight = 1.0 - projS.distance(ppoints[currP])/p[side].distance(ppoints[currP]);

		if (weight < 0.0)
			weight = 0.0;

		if ((voronoiW.find(currS) != voronoiW.end()) && (weight > 0.0))
		{
			sumDV += dist*voronoiW[currS];
//			sumDV += dist;
//			sumDW += dist*weight;
			sumW += weight;
			sumV += voronoiW[currS];
			sumWV += weight*voronoiW[currS];
//			sumV = 1.0;

			if (currP == 102)
				cout << " #" << currS << ":" << weight << "," << voronoiW[currS] << "," << dist;
		}
	}

	double dd = 0;

	if ((sumW != 0.0) && (sumV != 0.0))
//		dd = -sumDV/sumV/sumW;
		dd = -sumDV/sumWV;

	Point baseP = ppoints[currP] + nn*dd;

	// TEST: check if weighted Voronoi distances sum up to zero
	double sum = 0.0;

	for (auto entry:samplesMMap)
	{
		int currS = entry.second.first;
		int side = entry.second.second;
		Point projS;
		intersectLines(p[side], ppoints[currP] - p[side], ppoints[currS], normals[currP], projS);	// point on edge
		double dist = nn.dot(projS - ppoints[currS]);
		double weight = 1.0 - projS.distance(ppoints[currP])/p[side].distance(ppoints[currP]);

		if (weight < 0.0)
			weight = 0.0;

		if ((voronoiW.find(currS) != voronoiW.end()) && (weight > 0.0))
		{
			sum += (dist + dd*weight)*voronoiW[currS];
//			sum += dist + dd*weight;

			// DEBUG
//			if (currP == 2)
//				cout << "#" << currP << ": #" << currS << ": d=" << dist << ", w=" << weight << ", v=" << voronoiW[currS] << endl;
		}
	}

	// TEST
	if (currP == 2)
		cout << "";

	// check if neighbor angles mix concave and convex
	int prev2 = neighbors[prev][0];

	if (prev2 == currP)
		prev2 = neighbors[prev][1];

	int next2 = neighbors[next][0];

	if (next2 == currP)
		next2 = neighbors[next][1];

	Point vec0 = ppoints[prev] - ppoints[prev2];
	Point vec1 = ppoints[currP] - ppoints[prev];
	Point vec2 = ppoints[next] - ppoints[currP];
	Point vec3 = ppoints[next2] - ppoints[next];
	vec0.normalize();
	vec1.normalize();
	vec2.normalize();
	vec3.normalize();
	Point n0(vec0[1], -vec0[0]);
	Point n1(vec1[1], -vec1[0]);
	Point n2(vec2[1], -vec2[0]);
//	bool convex[3] = { (n0.dot(vec1) > 0.0), (n1.dot(vec2) > 0.0), (n2.dot(vec3) > 0.0) };

//	if ((convex[0] != convex[1]) || (convex[2] != convex[1]))
	if (false)
	{
		// is mixed-angle: compute intersection of normal with baseline and mix
		Point baseX;
		intersectLines(ppoints[currP], normals[currP], ppoints[prev], baseV, baseX);

		// DEBUG
		if (currP == 102)
//			cout << "";
			cout << " #" << currP << ": P:" << baseP[0] << "/" << baseP[1] << ", X:" << baseX[0] << "/" << baseX[1] << endl;

		baseP = baseP*0.5 + baseX*0.5;
	}

	// clamp inside disk
	dist = nn.dot(baseP - xx[0]);

	if (dist < 0.0)
		baseP = xx[0];
	else
	if (dist > len)
		baseP = xx[1];

	return baseP;
}

/*
 * denoise points: minimize energy function that removes noise from features
 */
void denoisePointsL2(vector<Point> &points, vector<PState> &pStates,
		vector<vector<int> > *nhood, vector<int[2]> &neighbors, vector<Point> &normals,
		vector<Circle> &circles, vector<float> &noise, vector<Point> &denoisedPoints)
{
//	return;

	const float TOLERANCE = 0.01;
	int i;
	bool convergedP[points.size()];
	vector<Point> outPoints = points;

	for (i = 0; i < (int)points.size(); i++)
		convergedP[i] = false;

	// determine noise extent as max of modeled input and connectivity-determined noise
	// TODO

	// iteratively adjust points
	bool converged = false;
	int iter = 0;

	while ((iter < 100) && !converged)
	{
		// DEBUG
		cout << "iteration #" << iter << ", converged:";

		converged = true;

		for (i = 0; i < (int)points.size(); i++)
		{
			if (pStates[i] == PState::MANIFOLD)	// TODO: also handle LEAF vertices
			{
				if (!convergedP[i])
				{
					Point oldP = denoisedPoints[i];
					outPoints[i] = minimizeVertex(points, denoisedPoints, nhood, neighbors, normals, noise, i);

					if (oldP.distance(outPoints[i]) < TOLERANCE)
					{
						convergedP[i] = true;

						// DEBUG
						cout << " " << i;
					}
					else
					{
						converged = false;

						// DEBUG
						if (i == 97)
							cout << " " << i << ": " << oldP.distance(outPoints[i]) << ", (" << outPoints[i][0] << "/"  << outPoints[i][1] << ")";
					}
				}
			}
		}

		denoisedPoints = outPoints;
		iter++;

		// DEBUG
		cout << endl;
	}

	cout << "not converged:";

	for (i = 0; i < (int)points.size(); i++)
		if ((pStates[i] == PState::MANIFOLD) && !convergedP[i])	// TODO: also handle LEAF vertices
			cout << " #" << i;

	cout << endl;
}

/*
 * collect all segments in nhood outside currP = connected neighbors which include MANIFOLD points
 */
void collectIntersectingSegmentsInNHood(int p, vector<vector<int> > *nhood, vector<PState> &pStates,
		vector<int[2]> &neighbors, vector<list<int> > &segments, bool &existsClosedCurve, int &openCurveCount)
{
	int i;
	set<int> manifoldPointSet, followPointSet;

	existsClosedCurve = false;
	openCurveCount = 0;

	// collect all MANIFOLD/CONF/LEAF points
	for (i = 0; i < 2; i++)
		for (auto n:nhood[i][p])
		{
			if (pStates[n] == PState::MANIFOLD)
				manifoldPointSet.insert(n);

			if ((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::CONFORM) || (pStates[n] == PState::LEAF))
				followPointSet.insert(n);
		}

	// merge nhood sets
	set<int> nhoodSet;

	for (i = 0; i < 2; i++)
		nhoodSet.insert(nhood[i][p].begin(), nhood[i][p].end());

	// from each MANIFOLD point, follow both directions to traverse CC, until all MANIFOLD/CONF/LEAF visited
	while (manifoldPointSet.size() > 0)
	{
		// get an initial MANIFOLD
		set<int>::iterator pIter = manifoldPointSet.begin();
		int startP = *pIter;
		manifoldPointSet.erase(pIter);
		followPointSet.erase(startP);
		list<int> segment;
		segment.push_back(startP);
		bool intersects[2] = { false, false };

		// follow neighbors of starting point in both directions
		for (i = 0; i < 2; i++)
		{
			int prevP = startP;
			int currP = neighbors[startP][i];

			// traverse neighbors while inside nhood
			while ((currP != -1) && (currP != startP) && (nhoodSet.find(currP) != nhoodSet.end()))
			{
				// add in correct order to list
				if (i == 0)
					segment.push_back(currP);
				else
					segment.push_front(currP);

				if (pStates[currP] == PState::MANIFOLD)
					manifoldPointSet.erase(currP);

				followPointSet.erase(currP);

				int nextP = neighbors[currP][(neighbors[currP][0] == prevP) ? 1 : 0];
				prevP = currP;
				currP = nextP;
			}

			if (currP == startP)
				cout << "";

			intersects[i] = ((currP != -1) && (currP != p) && (currP != startP) && (nhoodSet.find(currP) == nhoodSet.end()));
		}

		if (intersects[0] || intersects[1])
		{
			if (intersects[0] && intersects[1])
			{
				// if segment intersects at both ends, both end points must be MANIFOLD
				assert(((pStates[segment.front()] == PState::MANIFOLD) || (pStates[segment.front()] == PState::SHARP)) &&
						((pStates[segment.back()] == PState::MANIFOLD) || (pStates[segment.back()] == PState::SHARP)));
				existsClosedCurve = true;
			}
			else
			{
				// make sure that last point's neighbor intersects
				if (intersects[1])
					segment.reverse();

				// if only end intersects, last point must be MANIFOLD (first point is MANIFOLD iff it connects to currP)
				assert((pStates[segment.back()] == PState::MANIFOLD) || (pStates[segment.back()] == PState::SHARP));
				openCurveCount++;
			}

			segments.push_back(segment);
		}
	}
}

/*
 * test whether sharp condition fulfilled if exactly two (open) curves exist
 */
bool isSharpPoint(int curr, vector<Point> &points, vector<list<int> > &segments,
		vector<PState> &pStates, vector<int[2]> &neighbors, vector<vector<int> > *nhood)
{
	const float EPS = 0.0001;
	int i;
	float r = 0.0, angle[2];
	Point c = points[curr], x[2];

	// merge nhood sets including current point
	set<int> nhoodSet;

	for (i = 0; i < 2; i++)
		nhoodSet.insert(nhood[i][curr].begin(), nhood[i][curr].end());

	// compute radius of nhood circle
	for (auto n:nhoodSet)
	{
		float dist = sqrt(points[n].squared_distance(c));

		if (dist > r)
			r = dist;
	}

	nhoodSet.insert(curr);

	// test if open curves intersect nhood circle in a <60 degree angle
	for (i = 0; i < 2; i++)
	{
		int p0Index = segments[i].back();
		int p1Index = neighbors[p0Index][0];

		if (nhoodSet.find(p1Index) != nhoodSet.end())
			p1Index = neighbors[p0Index][1];

		Point p0 = points[p0Index];
		Point p1 = points[p1Index];
		Point v = p1 - p0;

		// intersect edge p0 + t*v with circle (c, r), s.t. t=[0..1]:
		// (p0x + t*v0x - c0x)^2 + (p0y + t*v0y - c0y)^2 = r^2
		float a = SQR(v[0]) + SQR(v[1]);
		float b = 2*(v[0]*(p0[0] - c[0]) + v[1]*(p0[1] - c[1]));
		Point pc = p0 - c;
		float cc = SQR(pc[0]) + SQR(pc[1]) - SQR(r);
		float det = SQR(b) - 4*a*cc;

//		if ((det < 0.0) && (det > -EPS))
//			det = 0.0;

		float t[2] = { (-b - sqrt(det))/(2*a), (-b + sqrt(det))/(2*a) };

		if ((t[0] < -EPS) || (t[0] > 1.0 + EPS))
			t[0] = t[1];

//		assert((t[0] >= -EPS) && (t[0] <= 1.0 + EPS));
		x[i] = p0 + v*t[0];

		// compute angle
		angle[i] = atan2(x[i][1] - c[1], x[i][0] - c[0])*180.0/PI;
	}

	float angleDiff = abs(angle[0] - angle[1]);

	if (angleDiff > 180.0)
		angleDiff = 360.0 - angleDiff;

	if (angleDiff > 60.0)
		return false;

	// test whether edges to new SHARP point intersect any edges incident to MANIFOLD/CONF/LEAF points in nhood
	set<pair<int, int> > edgeSet;

	for (auto n:nhoodSet)
		if ((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::CONFORM) || (pStates[n] == PState::LEAF))
		{
			for (i = 0; i < 2; i++)
			{
				pair<int, int> edge(n, neighbors[n][i]);

				if ((edge.second != -1) && ((pStates[edge.second] == PState::MANIFOLD) || (pStates[edge.second] == PState::CONFORM) || (pStates[edge.second] == PState::LEAF)))
				{
					if (edge.first > edge.second)
						swap(edge.first, edge.second);

					edgeSet.insert(edge);
				}
			}
		}

	Point xx;

	for (auto segment:segments)
	{
		int last = segment.front();
		Point lastP = points[last];

		for (auto edge:edgeSet)
		{
			if ((edge.first != last) && (edge.second != last) && intersectsEdge(c, lastP, points[edge.first], points[edge.second], xx))
				return false;
		}
	}

/*
	// test if half space opposite of intersection points contains points
	// compute half space line parallel to intersection points' line
	Point v = x[1] - x[0];
	Point n(-v[1], v[0]);

	if (n.dot(x[0] - c) < 0.0)
		n = Point(0.0, 0.0) - n;

	for (i = 0; i < 2; i++)
		for (auto currN:nhood[i][curr])
		{
			// test all points in nhood if in other half space
			Point currV = points[currN] - c;

			if (n.dot(currV) < 0.0)
				return true;
		}

	return false;
*/
	return true;
}

/*
 * reconstruct points assuming they are noisy
 */
bool Reconstruct2D::reconstructNoisy()
{
	int i, j, side;
	vector<vector<int> > nhood[3];

	for (i = 0; i < 3; i++)
		nhood[i].resize(points.size());

	vector<int[2]> neighbors(points.size());
	circles.resize(points.size());
	arcs.resize(points.size());
	normals.resize(points.size());
	projPoints.resize(points.size());
	denoisedPoints.resize(points.size());

	PrecTimer timer;

	timer.start();
	// create kd-tree
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points = annAllocPts(points.size(), 2);

	for (i = 0; i < (int)points.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = points[i][0];
		p[1] = points[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, points.size(), 2);

	set<int> unhandledPointSet;
	map<int, int> eliminatedMap;
	vector<PState> pStates(points.size());
	map<int, list<int> > nhoodRefsMap;

	for (i = 0; i < (int)points.size(); i++)
	{
		pStates[i] = PState::OUTLIER;
		neighbors[i][0] = -1;
		neighbors[i][1] = -1;

		// get 2 nearest neighbors per point
		ANNidxArray nnIdx = new ANNidx[3];
		ANNdistArray distances = new ANNdist[3];
		kdTree->annkSearch(ann_points[i], 3, nnIdx, distances);
		nhood[0][i].push_back(nnIdx[1]);
		nhood[1][i].push_back(nnIdx[2]);
	}

	// unmark non-outliers
	for (i = 0; i < (int)points.size(); i++)
	{
		pStates[nhood[0][i][0]] = PState::INITIAL;
		pStates[nhood[1][i][0]] = PState::INITIAL;
	}

	// TEST: unmark all outliers
	for (i = 0; i < (int)points.size(); i++)
		pStates[i] = PState::INITIAL;

	for (i = 0; i < (int)points.size(); i++)
		if (pStates[i] != PState::OUTLIER)
			unhandledPointSet.insert(i);

	if (DEBUG)
	{
		cout << "outliers:";

		for (i = 0; i < (int)points.size(); i++)
			if (pStates[i] == PState::OUTLIER)
				cout << " " << i;

		cout << endl;
	}

//	return;

	// process unhandled points (handled = OUTLIER, NONCONFORM, CONFORM, MANIFOLD, SHARP or LEAF) until none remain
	// unhandled points states can oscillate between NONCONFORM, CONFORM
	iterations = 0;
	handledFitCount = 0;
	handledPointCount = 0;
	squaredFitCount = 0;

	while ((unhandledPointSet.size() > 0) && (iterations < ((maxIter == -1) ? 99999 : maxIter)))	// TODO: limit of iterations
	{
		if (DEBUG)
			outputStateStatistics(iterations, unhandledPointSet, points, pStates);

		// while points exist which are NONCONF or CONF with 0-1 overlaps
		// NOTE: set is traversed in order while removing/adding arbitrary elements
		set<int>::iterator pIter = unhandledPointSet.begin();

		while (pIter != unhandledPointSet.end())
		{
			int currP = *pIter++;	// advance to next point in set
			PState oldState = pStates[currP];
			set<int> affectedPointSet;
			handledFitCount++;

			if (DEBUG)
			{
				cout << "handle #" << currP << ": ";

				if ((iterations == 1) && (currP == 47))
					cout << "";
			}

			// store old fit
			Circle oldCircle = circles[currP];
			Point oldNormal = normals[currP];
			vector<int> oldNhood[2];

			for (i = 0; i < 2; i++)
				oldNhood[i] = nhood[i][currP];

			// DEBUG
			if ((iterations == 10) && (currP == 99))
				cout << "";

#ifdef SHARPLEAF
			// LEAF candidate: check if point is CONFORM and has exactly 1 overlap, with >=3 MANIFOLD points
			bool leafOverlap = false;
			int leafOverlapSide = -1;

			if (pStates[currP] == PState::CONFORM)
			{
				if ((neighbors[currP][0] != -1) ^ (neighbors[currP][1] != -1))
				{
					// determine overlap side
					int side = (neighbors[currP][0] != -1) ? 0 : 1;
					leafOverlap = true;

					// test for M_COUNT successive MANIFOLD neighbors in nhood (= M_COUNT+2 edges)
					const int M_COUNT = 2;
					int count = M_COUNT;
					int currN = neighbors[currP][side];
					int prevN = currP;

					while (leafOverlap && (count > 0))
					{
						int nSide = (neighbors[currN][0] == prevN) ? 0 : 1;
						int currN2 = neighbors[currN][1 - nSide];
						leafOverlap = (pStates[currN2] == PState::MANIFOLD) && (find(nhood[side][currP].begin(), nhood[side][currP].end(), currN2) != nhood[side][currP].end());
						prevN = currN;
						currN = currN2;
						count--;
					}

					if (leafOverlap)
					{
						leafOverlapSide = side;

						if (DEBUG)
							cout << "OVERLAP ";
					}
				}
			}
#endif
			if (!increaseNHood(kdTree, ann_points, points, pStates, circles, nhood, neighbors, currP, normals[currP], nhoodRefsMap))
			{
				if (DEBUG)
					cout << "no more kNN ";

				if ((pStates[currP] == PState::CONFORM) && ((neighbors[currP][0] != -1) || (neighbors[currP][1] != -1)))
				{
#ifdef SHARPLEAF
					// test if currP has a single NONCONFORM point in non-overlap nhood
					int side = (neighbors[currP][0] != -1) ? 0 : 1;
					list<int> ncPoints;

					for (auto n:nhood[1 - side][currP])
						if (pStates[n] == PState::NONCONFORM)
							ncPoints.push_back(n);

					if (ncPoints.size() == 1)
					{
						// make that point the leaf instead
						int currN = ncPoints.front();
						neighbors[currP][1 - side] = currN;
						int nSide = (find(nhood[0][currN].begin(), nhood[0][currN].end(), currP) != nhood[0][currN].end()) ? 0 : 1;
						neighbors[currN][nSide] = currP;
						changePointState(currN, PState::LEAF, pStates, pIter, unhandledPointSet);
						affectedPointSet.insert(currN);
						changePointState(currP, PState::MANIFOLD, pStates, pIter, unhandledPointSet);
						affectedPointSet.insert(currP);
					}
					else
					{
						changePointState(currP, PState::LEAF, pStates, pIter, unhandledPointSet);
						affectedPointSet.insert(currP);
					}
#endif
				}
				else
				{
					bool isLeaf = false;

					// test if point is NONCONFORM and its next non-eliminated neighbor is a LEAF candidate
					if ((pStates[currP] == PState::NONCONFORM) && ((nhood[0][currP].size() == 0) ^ (nhood[1][currP].size() == 0)))
					{
						int side = (nhood[0][currP].size() == 0) ? 1 : 0;
						vector<int>::iterator nIter = nhood[side][currP].begin();

						while ((nIter != nhood[side][currP].end()) && (pStates[*nIter] == PState::ELIMINATED))
							nIter++;

						if (nIter != nhood[side][currP].end())
						{
							int currN = *nIter;

							isLeaf = ((pStates[currN] == PState::CONFORM) && ((neighbors[currN][0] != -1) ^ (neighbors[currN][1] != -1)));

							if (isLeaf)
							{
#ifdef SHARPLEAF
								int nSide = (neighbors[currN][0] == -1) ? 1 : 0;
								neighbors[currN][1 - nSide] = currP;
								neighbors[currP][side] = currN;
								changePointState(currN, PState::MANIFOLD, pStates, pIter, unhandledPointSet);
								affectedPointSet.insert(currN);
								changePointState(currP, PState::LEAF, pStates, pIter, unhandledPointSet);
								affectedPointSet.insert(currP);
#endif
							}
						}
					}

					if (!isLeaf)
						changePointState(currP, PState::ELIMINATED, pStates, pIter, unhandledPointSet);
				}
			}
#ifdef SHARPLEAF
			else	// neighborhood increased
			if (leafOverlapSide != -1)
			{
				bool isLeaf = false;

				// check if overlap consistency is lost
				set<int> consistentPointSet;

				determineConsistentPointSetInNhood(currP, leafOverlapSide, nhood, neighbors, pStates, consistentPointSet);

				// for all remaining handled points (MANIFOLD, SHARP or LEAF)
				for (auto currN:nhood[leafOverlapSide][currP])
				{
					if ((consistentPointSet.find(currN) == consistentPointSet.end()) &&
						((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::LEAF) ||
						(pStates[currN] == PState::SHARP)))
						isLeaf = true;
				}

				if (isLeaf && (nhood[1 - leafOverlapSide][currP].size() == 0))
				{
					// revert to old nhood and make it either a LEAF or SHARP
					circles[currP] = oldCircle;
					normals[currP] = oldNormal;

					for (i = 0; i < 2; i++)
						nhood[i][currP] = oldNhood[i];

					// make it a LEAF: test if currP has a single NONCONFORM point in non-overlap nhood
					list<int> ncPoints;

					for (auto n:nhood[1 - leafOverlapSide][currP])
						if (pStates[n] == PState::NONCONFORM)
							ncPoints.push_back(n);

					changePointState(currP, PState::LEAF, pStates, pIter, unhandledPointSet);
					affectedPointSet.insert(currP);

					// make neighbors consistent with nhood
					if (((neighbors[currP][0] != -1) && find(nhood[0][currP].begin(), nhood[0][currP].end(), neighbors[currP][0]) == nhood[0][currP].end()) ||
						((neighbors[currP][1] != -1) && find(nhood[1][currP].begin(), nhood[1][currP].end(), neighbors[currP][1]) == nhood[1][currP].end()))
						swap(neighbors[currP][0], neighbors[currP][1]);
				}
			}
#endif

			if (DEBUG)
				outputNeighborhood(currP, nhood);

			handledPointCount += 1 + nhood[0][currP].size() + nhood[1][currP].size();

			if (pStates[currP] == PState::CONFORM)
			{
				if (DEBUG)
					cout << "C ";

				// determine largest overlap in neighborhood
				int newN[2];

				for (side = 0; side < 2; side++)
					newN[side] = determinePointOverlapSide(currP, side, nhood, neighbors, pStates);

				// remove redundant overlaps
				for (side = 0; side < 2; side++)
				{
					// if previous overlap not included anymore
					int currN = neighbors[currP][side];

					// remove overlap
					if ((currN != -1) && (currN != newN[0]) && (currN != newN[1]))
						removeOverlap(currP, side, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
				}

				// add new overlaps
				for (i = 0; i < 2; i++)
				{
					// if overlap is new
					int nextN = newN[i];

					if ((nextN != -1) && (nextN != neighbors[currP][0]) && (nextN != neighbors[currP][1]))
					{
						if (DEBUG)
							cout << "+" << currP << "-" << nextN << " ";

						// add overlap
						if (neighbors[currP][0] == -1)
							side = 0;
						else
						{
							assert(neighbors[currP][1] == -1);
							side = 1;
						}

						neighbors[currP][side] = nextN;

						// locate side of nhood
						int nextSide = find(nhood[0][nextN].begin(), nhood[0][nextN].end(), currP) != nhood[0][nextN].end() ? 0 : 1;

						neighbors[nextN][nextSide] = currP;

						if (neighbors[currP][1 - side] != -1)
						{
							changePointState(currP, PState::MANIFOLD, pStates, pIter, unhandledPointSet);
							affectedPointSet.insert(currP);
						}

						if (neighbors[nextN][1 - nextSide] != -1)
						{
							changePointState(nextN, PState::MANIFOLD, pStates, pIter, unhandledPointSet);
							affectedPointSet.insert(nextN);
						}
					}
				}

				// make neighbors consistent with nhood
				if (((neighbors[currP][0] != -1) && find(nhood[0][currP].begin(), nhood[0][currP].end(), neighbors[currP][0]) == nhood[0][currP].end()) ||
					((neighbors[currP][1] != -1) && find(nhood[1][currP].begin(), nhood[1][currP].end(), neighbors[currP][1]) == nhood[1][currP].end()))
					swap(neighbors[currP][0], neighbors[currP][1]);

				if (DEBUG)
					cout << ", n=" << neighbors[currP][0] << "/" << neighbors[currP][1] << " ";

				// remove redundant points in overlaps
				for (side = 0; side < 2; side++)
				{
					// per overlap
					if (neighbors[currP][side] != -1)
					{
						// determine points shared inside overlaps of nhoods
						list<int> redundantPoints;
						addRedundantPointsInOverlap(currP, side, nhood, neighbors, unhandledPointSet, redundantPoints);

						// for all redundant points
						for (auto currN:redundantPoints)
						{
							// remove overlaps
							if ((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::CONFORM) ||
								(pStates[currN] == PState::SHARP))
								removeOverlaps(currN, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);

							// remove point
							changePointState(currN, PState::ELIMINATED, pStates, pIter, unhandledPointSet);
						}
					}
				}

				// remove inconsistent points in nhood for overlaps
				for (side = 0; side < 2; side++)
				{
					// per overlap
					if (neighbors[currP][side] != -1)
					{
						removeInconsistentPointSetInNhood(currP, side, nhood, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);

						// also for nhood of overlap neighbor
						int currN = neighbors[currP][side];
						int nSide = (neighbors[currN][0] == currP) ? 0 : 1;
						removeInconsistentPointSetInNhood(currN, nSide, nhood, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
					}
				}

				// detect 1-point obstacle to second overlap (projections of overlaps overlap and cannot be resolved otherwise)
				if ((neighbors[currP][0] != -1) ^ (neighbors[currP][1] != -1))
				{
					bool removePoint = false;

					// determine side with overlap
					int side = (neighbors[currP][0] != -1) ? 0 : 1;

					// locate CONFORM point between currP and its overlap neighbor
					vector<int>::iterator iter = nhood[side][currP].begin();

					while ((iter != nhood[side][currP].end()) && (*iter != neighbors[currP][side]) && (pStates[*iter] != PState::CONFORM))
						iter++;

					// check if it has one overlap, with neighbor in other nhood
					if ((iter != nhood[side][currP].end()) && (*iter != neighbors[currP][0]) && (pStates[*iter] == PState::CONFORM) && ((neighbors[*iter][0] != -1) ^ (neighbors[*iter][1] != -1)))
					{
						int nSide = (neighbors[*iter][0] != -1) ? 0 : 1;
						removePoint = (find(nhood[1 - side][currP].begin(), nhood[1 - side][currP].end(), neighbors[*iter][nSide]) != nhood[1 - side][currP].end());
					}

					if (removePoint)
					{
						// remove currP
						removeOverlap(currP, side, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
						changePointState(currP, PState::ELIMINATED, pStates, pIter, unhandledPointSet);
					}
				}
			}
			else
			if (pStates[currP] == PState::NONCONFORM)
			{
				if (DEBUG)
					cout << "N ";
			}

			// remove points projecting inside boundary edge
			if (((pStates[currP] == PState::CONFORM) || (pStates[currP] == PState::NONCONFORM)) && isInsideEdge(currP, neighbors, nhood, pStates))
			{
				// if was CONFORM before, set state temporarily to remove overlaps
				if ((pStates[currP] == PState::NONCONFORM) && (oldState == PState::CONFORM))
					pStates[currP] = PState::CONFORM;

				if ((pStates[currP] == PState::MANIFOLD) || (pStates[currP] == PState::CONFORM) || (pStates[currP] == PState::SHARP))
					removeOverlaps(currP, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);

				// remove point
				changePointState(currP, PState::ELIMINATED, pStates, pIter, unhandledPointSet);
				affectedPointSet.insert(currP);
			}
#ifdef SHARPLEAF
//#ifdef NEW
			// NEW: for NONCONFORM points, check if whether curve segment exist in nhood -> eliminate or connect
			if (pStates[currP] == PState::NONCONFORM)
			{
				// determine intersecting segments (MANIFOLD point at at least one end) inside nhood
				vector<list<int> > segments;
				bool existsClosedCurve = false;
				int openCurveCount = 0;

				// DEBUG
//				if ((iterations == 26) && (currP == 34))
				if ((iterations == 71) && (currP == 99))
					cout << "";

				collectIntersectingSegmentsInNHood(currP, nhood, pStates, neighbors, segments, existsClosedCurve, openCurveCount);

				// remove points which have a closed segment (all MANIFOLD points), which intersects twice with nhood disc
				if (existsClosedCurve)
				{
/*
					bool isConnecting = false;

					if (openCurveCount >= 1)
					{
						for (auto segment:segments)
						{
							int startP = segment.front();

							if ((neighbors[startP][0] == currP) || (neighbors[startP][1] == currP))
								isConnecting = true;
						}
					}

					if (isConnecting)
					{
						// make point a LEAF
						changePointState(currP, PState::LEAF, pStates, pIter, unhandledPointSet);
						affectedPointSet.insert(currP);
					}
					else
					{
						// eliminate current point and overlaps if exist
						for (int side = 0; side < 2; side++)
						{
							int currN = neighbors[currP][side];

							if (currN != -1)
							{
								int nSide = (neighbors[currN][0] == currP) ? 0 : 1;
								removeOverlap(currN, nSide, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
							}
						}

						changePointState(currP, PState::ELIMINATED, pStates, pIter, unhandledPointSet);
					}
*/
				}
				else	// then only try to connect unreferenced points
//				if ((neighbors[currP][0] == -1) && (neighbors[currP][1] == -1))
				{
					// connect branches which have a manifold chain intersecting the nhood disc each
					// if 2 segment which start with LEAF or CONF and end with MANIFOLD
					if (openCurveCount == 2)
					{
						assert(segments.size() == 2);	// TODO: else select the 2 open curves

						// test if each open curve has >=2 points
//						if ((segments.front().size() >= 2) && (segments.back().size() >= 2))

						// test whether sharp condition fulfilled
						if (isSharpPoint(currP, points, segments, pStates, neighbors, nhood))
						{
							// first, remove old overlaps if existing
							if (oldState == PState::CONFORM)
							{
								pStates[currP] = PState::CONFORM;
								removeOverlaps(currP, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
							}

							// mark current point as MANIFOLD (TODO: mark current point as SHARP)
							changePointState(currP, PState::SHARP, pStates, pIter, unhandledPointSet);
							affectedPointSet.insert(currP);

							// mark both LEAF/CONF of the open curves as MANIFOLD and connect them
							int pSide = 0;

							for (auto segment:segments)
							{
								int lastP = segment.front();

								if (pStates[lastP] != PState::MANIFOLD)
								{
									// mark point as MANIFOLD
									changePointState(lastP, PState::MANIFOLD, pStates, pIter, unhandledPointSet);
									affectedPointSet.insert(lastP);

									// connect neighbors
									int nSide = (neighbors[lastP][0] == -1) ? 1 : 0;
									neighbors[lastP][1 - nSide] = currP;
								}

								neighbors[currP][pSide++] = lastP;
							}
						}
					}
					else
					if (openCurveCount >= 2)
					{
						// if >2 -> create T-junction
//						assert(false);	// TODO
					}
				}
			}
//#endif
#endif

			// while affected points remain, remove overlaps affected by them - to ensure consistency
			while (affectedPointSet.size() > 0)
			{
				set<int>::iterator p2Iter = affectedPointSet.begin();

				while (p2Iter != affectedPointSet.end())
				{
					int curr2P = *p2Iter;	// advance to next point in set
					affectedPointSet.erase(p2Iter++);
					set<int> refSet;
					refSet.insert(nhoodRefsMap[curr2P].begin(), nhoodRefsMap[curr2P].end());

					for (auto n:refSet)
						if ((neighbors[n][0] != -1) || (neighbors[n][1] != -1))	// test only points with overlaps
						{
							bool isConsistent = true;

							// DEBUG
							if ((iterations == 1) && (currP == 91) && (n == 91))
								cout << "";

							// verify that SHARP points have both MANIFOLD neighbors
							if ((pStates[n] == PState::SHARP) &&
								((pStates[neighbors[n][0]] != PState::MANIFOLD) || (pStates[neighbors[n][1]] != PState::MANIFOLD)))
								isConsistent = false;

							// verify that neighbor is on same side as nhood (except for SHARP points and SHARP neighbors)
							if (isConsistent && (pStates[n] != PState::SHARP))
								for (side = 0; side < 2; side++)
								{
									int currN = neighbors[n][side];

									if ((currN != -1) && (pStates[currN] != PState::SHARP) &&
										(find(nhood[side][n].begin(), nhood[side][n].end(), currN) == nhood[side][n].end()))
										isConsistent = false;
								}

							// test overlaps of n if they are consistent (except for SHARP points)
							if (isConsistent && (pStates[n] != PState::SHARP))
								for (side = 0; side < 2; side++)
								{
									if ((neighbors[n][side] != -1) && (pStates[neighbors[n][side]] != PState::SHARP))
									{
										set<int> consistentPointSet;

										determineConsistentPointSetInNhood(n, side, nhood, neighbors, pStates, consistentPointSet);

										// for all remaining handled points (MANIFOLD, LEAF or SHARP)
										for (auto currN:nhood[side][n])
											if ((consistentPointSet.find(currN) == consistentPointSet.end()) &&
												((pStates[currN] == PState::MANIFOLD) || (pStates[currN] == PState::LEAF) ||
												(pStates[currN] == PState::SHARP)))
												isConsistent = false;
									}
								}

							if (!isConsistent)
							{
								if (DEBUG)
									cout << "REMOVE OVERLAPS for " << n << endl;

								removeOverlaps(n, neighbors, pStates, pIter, unhandledPointSet, affectedPointSet);
							}
						}
				}
			}

			// keep overlaps and old neighborhood (since neither SHARP nor projects inside boundary edge)
			if ((pStates[currP] == PState::NONCONFORM) && (oldState == PState::CONFORM))
			{
				// restore old neighborhood
				circles[currP] = oldCircle;
				normals[currP] = oldNormal;

				for (i = 0; i < 2; i++)
					nhood[i][currP] = oldNhood[i];

				pStates[currP] = PState::CONFORM;

				if (DEBUG)

					cout << "reverted ";
			}

			// check if next point still in set (in case it was erased during this loop)
			if (unhandledPointSet.find(*pIter) == unhandledPointSet.end())
				pIter = unhandledPointSet.end();

			if (DEBUG)
				cout << endl;

			// ASSERTION
			for (i = 0; i < (int)points.size(); i++)
				checkRuleForPoint(i, pStates, nhood, neighbors);
		}

		for (i = 0; i < (int)points.size(); i++)
		{
			Point center(circles[i].a, circles[i].b);
			Point centerNormal(points[i] - center);
			centerNormal.normalize();

			// determine arc end points in degrees
			int deg[3];	// start/center/end
			Point pNormal[3] = { points[nhood[0][i].back()] - center, centerNormal, points[nhood[1][i].back()] - center };

			for (j = 0; j < 3; j++)
			{
				pNormal[j].normalize();
				deg[j] = (int)((computeDegreesFromVector(pNormal[j])+PI)/PI*180);
			}

			// reverse if not oriented in positive direction
			if (((abs(deg[1] - deg[0] + 360) % 360) > 180) || ((abs(deg[2] - deg[1] + 360) % 360) > 180))
				swap(deg[0], deg[2]);

			arcs[i] = pair<int, int>(deg[0], deg[2]);
		}

		iterations++;
	}

	outputSize = 0;

	for (i = 0; i < (int)points.size(); i++)
	{
		if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
			pClasses[i] = FITTED;	// has one or two neighbors
		else
		if (pStates[i] == PState::OUTLIER)
			pClasses[i] = OUTLIER;
		else
			pClasses[i] = NONFITTED;
	}

	// only display arcs of points in boundary
	for (i = 0; i < (int)points.size(); i++)
		if ((pStates[i] != PState::MANIFOLD) && (pStates[i] != PState::SHARP) && (pStates[i] != PState::LEAF))
			circles[i].r = 0.0;

	// insert eliminated points
	// TODO: which n'hood to assign, or use manifold points' n'hoods? assign into maps when eliminating, then reverse-reassign
	if (DEBUG)
		for (auto entry:eliminatedMap)
			cout << "#" << entry.first << ": referenced from #" << entry.second << endl;

//	insertEliminatedPoints(eliminatedMap, points, neighbors);

	// project points
	projectPoints(points, pStates, nhood, neighbors, normals, circles, projPoints);	// needs original normals

	if (!(mode & MODE_BLEND))
	{
		// do not project points
		for (i = 0; i < (int)points.size(); i++)
			projPoints[i] = points[i];
	}

	// consistently orient directions and normals
	bool isClosed = orientManifoldPoints(points, pStates, nhood, neighbors, normals, circles);

	// dummy initialize for testing
	for (i = 0; i < (int)points.size(); i++)
		denoisedPoints[i] = points[i];

	// denoise points
//	denoisePointsL2(points, pStates, nhood, neighbors, normals, circles, noise, denoisedPoints);
/*
	for (i = 0; i < (int)points.size(); i++)
		projPoints[i] = denoisedPoints[i];
*/
	timer.stop();

	for (i = 0; i < (int)points.size(); i++)
	{
		if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
		{
			outputPoints.push_back(projPoints[i]);
			outputSize++;
		}
	}

	for (i = 0; i < (int)points.size(); i++)
		if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
			squaredFitCount += SQR(1 + nhood[0][i].size() + nhood[1][i].size());

	runtime = timer.getReal();

	if (DEBUG)
		cout << "handled " << handledPointCount << " operations in " << handledFitCount << " fits for " <<
			points.size() << " samples in " << iterations << " iterations, vs. " << squaredFitCount <<
			" output complexity, in " << runtime << "s" << endl;

	for (i = 0; i < (int)points.size(); i++)
	{
		if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::SHARP) || (pStates[i] == PState::LEAF))
		{
			// draw edges
			for (j = 0; j < 2; j++)
			{
				int n = neighbors[i][j];

				if ((n != -1) && ((pStates[n] == PState::MANIFOLD) || (pStates[n] == PState::SHARP) || (pStates[n] == PState::LEAF)))
				{
					pair<int, int> edge = make_pair(i, n);

					if (edge.first > edge.second)
						swap(edge.first, edge.second);

					visEdgeMap[edge] = BIJECTIVE;
				}
			}
		}

		// classify points
		if ((pStates[i] == PState::MANIFOLD) || (pStates[i] == PState::LEAF))
			pClasses[i] = FITTED;
		else
		if (pStates[i] == PState::SHARP)
			pClasses[i] = SHARP;
		else
		if (pStates[i] == PState::OUTLIER)
			pClasses[i] = OUTLIER;
		else
			pClasses[i] = NONFITTED;
	}

	return isClosed;
}

map<pair<int, int>, EdgeEnum> Reconstruct2D::getEdgeMap()
{
	return visEdgeMap;
}

vector<Point> Reconstruct2D::getProjectedPoints()
{
	return projPoints;
}

vector<Point> Reconstruct2D::getOutputPoints()
{
	return outputPoints;
}

vector<Point> Reconstruct2D::getDenoisedPoints()
{
	return denoisedPoints;
}

vector<Point> Reconstruct2D::getNormals()
{
	return normals;
}

vector<Circle> Reconstruct2D::getCircles()
{
	return circles;
}

vector<pair<int, int> > Reconstruct2D::getArcs()
{
	return arcs;
}

vector<PointsEnum> Reconstruct2D::getPointClassification()
{
	return pClasses;
}
