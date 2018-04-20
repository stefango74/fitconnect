//============================================================================
// Author      : Stefan Ohrhallinger
// Version     :
// Copyright   : GPL v3
// Description : Test driver for: Reconstruct a curve from noisy 2D points
//============================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <png++/png.hpp>
#include <GL/glew.h>
#include <GL/glut.h>

#include "Reconstruct2D.h"

using namespace std;

#define INIT_WIDTH 704	// for screen height of 768 not to squash square window
#define INIT_HEIGHT 704
int viewWidth = INIT_WIDTH;
int viewHeight = INIT_HEIGHT;

template<typename T> inline T CUB(T t) { return (SQR(t)*(t)); };

enum class PointState: char { EPSHALF, EPSONE, EPSLARGE, EPSHALF_NM, EPSONE_NM, EPSLARGE_NM };

vector<float> noise;

class TestReconstruct2D
{
private:
	int viewWidth, viewHeight;
	float scale;
	float offset[2], lfsMin, lfsMax;
	bool isClosed;
	int outputSize, iterations, handledFits, handledPoints, squaredFits;
	float runtime;
	vector<Point> points, projPoints, curvePoints, extremaPoints, mPoints, outputPoints, denoisedPoints;
	vector<float> lfsCurvePoints, kCurvePoints;
	vector<PointState> pointStates;
	vector<Point> normals;
	vector<PointsEnum> pClasses;
	map<pair<int, int>, EdgeEnum> edgeMap;
	vector<Circle> circles;
	vector<pair<int, int> > arcs;

	void drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap);
	void drawPoints(vector<Point> &points);
	void drawPoints(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawPointsWithStates(vector<Point> &points, vector<PointState> &pointStates);
	void drawPointsWithClasses(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawCurvePoints(vector<Point> &points, vector<float> &lfsCurvePoints, float minLfs, float maxLfs);
	void drawCircle(Circle &circle);
	void drawCircles(vector<Circle> &circles);
	void drawArc(Circle &c, pair<int, int> arc);
	void drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arc);
	void drawCover(Circle &c, pair<int, int> arc);
	void drawCovers(vector<Circle> &circles, vector<pair<int, int> > &arc);
	void drawLabels(vector<Point> &points);
	void drawNormals(vector<Point> &points, vector<PointsEnum> &pClasses);
	void drawNoiseExtent(vector<Point> &points, vector<float> &noise);
	Point transform(Point);
	float scaleOfPoints(float *translate, vector<Point> &points);

public:
	TestReconstruct2D();
	virtual ~TestReconstruct2D();
	void reshape(int width, int height);
	void draw(int items, float pointSize);
	void draw(string filename, int items, float pointSize);
	void loadPointSet(string name);
	void invertY();
	void loadCurvePointsFromPNG(char *name);
	void generateCurveFromBezierString(string str);
	void generateCircle();
	void generatePointSetCosinesOnCircleWithVaryingAmplitude(float minAmp, float maxAmp);
	void generatePointSetCosinesOnCircleWithVaryingFrequency(float minFreq, float maxFreq);
	void sampleCurveByEps(float maxEps, float error, bool isClosed, float perturb);
	void sampleCurveByReach(float maxReach, float error, bool isClosed, float perturb);
	void generatePointSetBySamplingCurveEpsAndManifold(float error);
	void generateCurveCosinesOnCircleWithVaryingFrequency(int count, float minFreq, float maxFreq, float amplitude);
	bool reconstruct(int mode, int maxIter);
	void reconstruct(vector<float> &noise, int mode, int maxIter);
	void setMaxIter(int p_maxIter);
	void computeScale();
	vector<Point> getPoints();
	vector<Point> getProjPoints();
	map<pair<int, int>, EdgeEnum> getEdges();
	void getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtime);
	void perturbPoints(float stddev);
	void perturbPointsMax(float max);
	void addCornerOutliers();
	void addNoiseToPoints(int multiplier, float maxDistance);
	void determinePointStates(float maxEps);
	float scaleOfCurvePoints();
	void generatePointsFromCircle();
	void generatePointsFrom2Circles(float noise, int seed);
	void generateHoleOutliers();
	void generateNoisySharpCorner(float aspectRatio, float noise, int seed);
	void generateNoisyLeafLine(float noise, int seed);
	void generateNoisyLeafFeature(float noise, int seed);
	void generateNoisySpiral(float noise, int seed);
	void generateVaryingNoisyCircle(float noise, int seed);
	void generateNoisyCircle(float noise, int seed);
	void generateZigZagRectangle();
	void generateZigZagLine();
	void generateZigZagCircle();
	void generateNoisyHole(float noise, float holeSize, int seed);
	void generateNoisyTCrossing(float noise, int seed);
	void generateNoisyFeatures(float noise, int seed);
	void generateOpenCurve();
	void generateTCurve();
	void generateSharpCorner();
	void generateHighNoise();
	void generateHighNoiseClosed();
	void generateOpenCurves();
	void reconstruct_otr(int iterations);
};

/*
 * determine scale of curve points for unit square
 */
float TestReconstruct2D::scaleOfCurvePoints()
{
	float offset[2];

	return scaleOfPoints(offset, curvePoints);
}

void TestReconstruct2D::addCornerOutliers()
{
	points.push_back(Point(-0.9, -2.9));
	points.push_back(Point(4.9, 2.9));
}

vector<Point> TestReconstruct2D::getPoints()
{
	return points;
}

vector<Point> TestReconstruct2D::getProjPoints()
{
	return projPoints;
}

map<pair<int, int>, EdgeEnum> TestReconstruct2D::getEdges()
{
	return edgeMap;
}

void TestReconstruct2D::getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtime)
{
	p_output = outputSize;
	p_iter = iterations;
	p_fit = handledFits;
	p_point = handledPoints;
	p_squared = squaredFits;
	p_runtime = runtime;
}

/*
 * determine scale of points for unit square
 */
float TestReconstruct2D::scaleOfPoints(float *offset, vector<Point> &points)
{
	int i;
	float min[2] = { numeric_limits<float>::max(), numeric_limits<float>::max() };
	float max[2] = { -numeric_limits<float>::max(), -numeric_limits<float>::max() };

	for (auto p:points)
	{
		for (i = 0; i < 2; i++)
		{
			if (p[i] < min[i])
				min[i] = p[i];

			if (p[i] > max[i])
				max[i] = p[i];
		}
	}

	float dim[2] = { max[0] - min[0], max[1] - min[1] };

	for (i = 0; i < 2; i++)
		offset[i] = min[i];

	i = (dim[0] > dim[1]) ? 0 : 1;

	offset[1 - i] -= (dim[i] - dim[1 - i])/2;

	return dim[i];
}

bool TestReconstruct2D::reconstruct(int mode, int maxIter)
{
	Reconstruct2D *instance = new Reconstruct2D(points, noise, mode);
	instance->setMaxIter(maxIter);
	bool isClosed = instance->reconstructNoisy();
	instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtime);
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	outputPoints = instance->getOutputPoints();
	denoisedPoints = instance->getDenoisedPoints();

	delete instance;

	return isClosed;
}

void TestReconstruct2D::reconstruct(vector<float> &p_noise, int mode, int maxIter)
{
	noise = p_noise;
	Reconstruct2D *instance = new Reconstruct2D(points, noise, mode);
	instance->setMaxIter(maxIter);
	instance->reconstructNoisy();
	instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtime);
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	outputPoints = instance->getOutputPoints();
	denoisedPoints = instance->getDenoisedPoints();

	delete instance;
}

void TestReconstruct2D::computeScale()
{
	scale = scaleOfPoints(offset, points);
}

/*
 * create a random float in the unit range [0..1]
 */
double unitrand()
{
	return (double)rand()/RAND_MAX;
}

void TestReconstruct2D::perturbPoints(float range)
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
	{
		float sigma = range*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += sigma*cos(angle);
		points[i][1] += sigma*sin(angle);
	}
}

/*
 * load point set from file with name
 */
void TestReconstruct2D::loadPointSet(string name)
{
	int i;

	points.clear();
	ifstream file(name);

	if (file)
	{
		while (file)
		{
			float x, y;
			file >> x >> y;

			if (file)
				points.push_back(Point(x, y));
		}
	}
	else
	{
		cerr << "ERROR: input file " << name << " could not be read." << endl;
		exit(2);
	}

	if (points.size() < 3)
	{
		cerr << "ERROR: input file " << name << " contains less than 3 points." << endl;
		exit(3);
	}

	for (i = 0; i < (int)points.size(); i++)
	{
		pClasses.push_back(CONFORMING);
		pointStates.push_back(PointState::EPSLARGE_NM);
	}

	computeScale();
}

/*
 * invert y-coordinate
 */
void TestReconstruct2D::invertY()
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
		points[i][1] = -points[i][1];
}

/*
 * converts a string two comma-separated floats into a Point
 */
Point str2point(string str)
{

	int ofs = str.find(",");
//	float t = stof(str);

//	string test = str.substr(0, ofs);
//	string test2 = str.substr(ofs + 1, str.length() - ofs - 1);
//	float t = stof(test);

	return Point(stof(str.substr(0, ofs)), stof(str.substr(ofs + 1, str.length() - ofs - 1)));
}

/*
 * parses SVG curveto string and extracts bezier control points
 */
void parseSVGCurveToString(string str, vector<Point> &bezierVec)
{
	int i;
	Point p[2];
	assert((str[0] == 'M') || (str[0] == 'm'));
	string::size_type ofs = str.find(' ', 2);
	string pStr = str.substr(2, ofs - 1);
	p[0] = str2point(pStr);
	char cc = str[ofs + 1];
	bool isRelative = (cc == 'c');
	assert(isRelative || (cc == 'C'));
	bool isEnd = false;
	ofs += 2;

	do
	{
		bezierVec.push_back(p[0]);

		for (i = 0; i < 3; i++)
		{
			int prevOfs = ofs;
			ofs = str.find(' ', ofs + 1);

			if (ofs == string::npos)
			{
				isEnd = true;
				ofs = str.length();
			}

			pStr.clear();
			pStr = str.substr(prevOfs, ofs - prevOfs);
			p[1] = str2point(pStr);

			if (isRelative)
				p[1] = p[1] + p[0];

			bezierVec.push_back(p[1]);
		}

		if (!isEnd)
			p[0] = p[1];

	} while (!isEnd);
}

/*
 * generate curve point set from SVG curveto string
 */
void TestReconstruct2D::generateCurveFromBezierString(string str)
{
	int i, j;
	const int SAMPLE_COUNT = 300;
	string bezierStr(str);
	vector<Point> b;
	parseSVGCurveToString(bezierStr, b);

	// iterate all cubic bezier curves
	for (j = 0; j < (int)b.size()/4; j++)
	{
		// sample the cubic bezier curve by evaluating with parameter t [0..1]
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			float t = (float)i/SAMPLE_COUNT;
			Point p = b[j*4]*CUB(1 - t) + b[j*4 + 1]*3*t*SQR(1 - t) + b[j*4 + 2]*3*SQR(t)*(1 - t) + b[j*4 + 3]*CUB(t);
			curvePoints.push_back(p);
		}
	}
}

/*
 * generate point set from sharp corner, perturb points by noise
 */
void TestReconstruct2D::generateNoisySharpCorner(float aspectRatio, float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float x = (float)i/SAMPLE_COUNT;
		float y = (float)i/SAMPLE_COUNT*aspectRatio;
		points.push_back(Point(x - 1.0, y));
		points.push_back(Point(x - 1.0, -y));
		points.push_back(Point(x, aspectRatio - y));
		points.push_back(Point(x, -aspectRatio + y));
	}

	perturbPoints(noise);
}

/*
 * generate point set from line, perturb points by noise
 */
void TestReconstruct2D::generateNoisyLeafLine(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));

	perturbPoints(noise);
}

/*
 * generate point set from line with feature, perturb points by noise
 */
void TestReconstruct2D::generateNoisyLeafFeature(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25, FEATURE_COUNT = 5;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));

	for (i = 0; i < FEATURE_COUNT; i++)
		points.push_back(Point(0.0, (float)i/SAMPLE_COUNT));


	perturbPoints(noise);
}

/*
 * generate point set from spiral, perturb points by noise
 * distance between boundary = 1
 * density by angle, therefore features denser sampled in interior
 */
void TestReconstruct2D::generateNoisySpiral(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;
	const int rotations = 2;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float i2 = 1.0 - SQR(1.0 - (float)i/(SAMPLE_COUNT*1.2));	// same range as i/SAMPLE_COUNT but growing superlinearly
		float angle = rotations*2*PI*i2;
		float x = (0.0 + (float)i2*rotations)*cos(angle);
		float y = (0.0 + (float)i2*rotations)*sin(angle);
		points.push_back(Point(x, y));
	}

	perturbPoints(noise);
}

/*
 * generate point set from circle, perturb points by varying noise
 */
void TestReconstruct2D::generateVaryingNoisyCircle(float noiseExtent, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;
	const int CLUSTER_COUNT = 2;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);
		points.push_back(Point(x, y));
	}

	// perturb with varying noise up to stddev 'noise'
	noise.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
	{
		float level = abs(1.0 - fmod(CLUSTER_COUNT*2*(float)i/points.size(), 2.0));
		float sigma = noiseExtent*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += level*sigma*cos(angle);
		points[i][1] += level*sigma*sin(angle);
		noise[i] = level*noiseExtent*scale;
	}
}

/*
 * generate point set from circle, perturb points by fixed noise
 */
void TestReconstruct2D::generateNoisyCircle(float noiseExtent, int seed)
{
	int i;
	const int SAMPLE_COUNT = 100;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);
		points.push_back(Point(x, y));
	}

	// perturb with varying noise up to stddev 'noise'
	noise.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
	{
		float sigma = noiseExtent*scale*unitrand();
		float angle = 2*PI*unitrand();
		points[i][0] += sigma*cos(angle);
		points[i][1] += sigma*sin(angle);
		noise[i] = noiseExtent*scale;
	}
}

void TestReconstruct2D::generateZigZagRectangle()
{
	int i;
	const int SIDE_COUNT = 4;
	const int SAMPLE_COUNT = 4*SIDE_COUNT;	// 4 sides
	const float STEP = 1.0;

	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SIDE_COUNT; i++)
	{
		int odd = i & 1;
		points[i] = Point(STEP/2 + i*STEP, STEP/2 - odd);
		points[SIDE_COUNT + i] = Point(STEP/2 + SIDE_COUNT*STEP + odd, STEP/2 + i*STEP);
		points[2*SIDE_COUNT + i] = Point(STEP/2 + (SIDE_COUNT - i)*STEP, STEP/2 + SIDE_COUNT*STEP + odd);
		points[3*SIDE_COUNT + i] = Point(STEP/2 - odd, STEP/2 + (SIDE_COUNT - i)*STEP);
	}

	for (i = 0; i < SAMPLE_COUNT; i++)
		noise[i] = STEP/2;
}

#ifdef OLD
void TestReconstruct2D::generateZigZagCircle()
{
	int i;
	const int SAMPLE_COUNT = 100;
	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		points[i] = Point(sin(i*2.0*PI/SAMPLE_COUNT), cos(i*2.0*PI/SAMPLE_COUNT));
		noise[i] = 0.1;
	}
}
#endif

void TestReconstruct2D::generateZigZagLine()
{
	int i;
	const int SAMPLE_COUNT = 10;
	const float STEP = 1.0;

	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
		points[i] = Point(STEP/2 + i*STEP, STEP/2 - (i & 1));

	for (i = 0; i < SAMPLE_COUNT; i++)
		noise[i] = STEP/2;
}

/*
 * generate point set from circle with non-random zig zag perturbation
 */
void TestReconstruct2D::generateZigZagCircle()
{
	int i;
	const int SAMPLE_COUNT = 40;
	const float NOISE = 0.1;
	points.resize(SAMPLE_COUNT);
	noise.resize(SAMPLE_COUNT);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float angle = 2*PI/SAMPLE_COUNT*i;
		float x = cos(angle);
		float y = sin(angle);

		if (i & 1)
			points[i] = Point(x, y)*(1.0 + NOISE);
		else
			points[i] = Point(x, y)*(1.0 - NOISE);

		noise[i] = NOISE;
	}
}

/*
 * generate point set from spiral, perturb points by noise
 */
void TestReconstruct2D::generateNoisyTCrossing(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 25;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		if ((i & 1) == 0)
			points.push_back(Point((float)i/SAMPLE_COUNT, 0.0));
		else
			points.push_back(Point(0.5, (float)i/SAMPLE_COUNT));
	}

	perturbPoints(noise);
}

/*
 * generate sine curve with frequency FREQ and increasing amplitude, perturb points by noise
 */
void TestReconstruct2D::generateNoisyFeatures(float noise, int seed)
{
	int i;
	const int SAMPLE_COUNT = 50;
	const int FREQ = 5;
	const float AMP = 0.1;

	srand(seed);

	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float x = ((float)i)/SAMPLE_COUNT*2 - 1.0;
		float y = AMP*i*sin((float)i/SAMPLE_COUNT*2*PI*FREQ)/SAMPLE_COUNT;
		points.push_back(Point(x, y));
	}

	for (i = 0; i < (int)points.size(); i++)
		points[i][1] += noise*(2*(double)rand()/RAND_MAX - 1.0);
}

/*
 * generate point set from 2 circles
 */
void TestReconstruct2D::generatePointsFrom2Circles(float noise, int seed)
{
	int i, j;
	const int SAMPLE_COUNT = 10;

	srand(seed);

	for (j = 0; j < 2; j++)
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			double t = (double)i/SAMPLE_COUNT*2*PI;
			Point p(4*j + sin(t), cos(t));
			points.push_back(p);
		}

	// add outliers
	points.push_back(Point(2.0, -1.1));
	points.push_back(Point(2.0, 1.1));

	computeScale();
	perturbPoints(noise);
}

/*
 * generate point set from 2 circles inside each other, with outliers
 */
void TestReconstruct2D::generateHoleOutliers()
{
	int i, j;
	const int SAMPLE_COUNT = 10;

	for (j = 0; j < 2; j++)
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			double t = (double)i/SAMPLE_COUNT*2*PI;
			double size = (0.6 + j);
			Point p(size*sin(t), size*cos(t));
			points.push_back(p);
		}

	// add outliers
	points.push_back(Point(-1.1, -2.0));
	points.push_back(Point(1.1, -2.0));

	computeScale();
}

/*
 * generate closed high noise point set
 */
void TestReconstruct2D::generateHighNoiseClosed()
{
	float coords[34][2] = {
			{ 0.0, 0.0 },
			{ 1.0, 0.1 },
			{ 2.0, -0.2 },
			{ 8.0, 0.2 },
			{ 9.0, -0.3 },
			{ 10.0, -0.1 },
			{ 0.0, 2.0 },
			{ 0.01, 4.0 },
			{ 0.0, 7.0 },
			{ 4.0, 7.01 },
			{ 7.0, 7.01 },
			{ 10.0, 7.0 },
			{ 10.01, 4.0 },
			{ 10.0, 2.0 },
		};

	srand(0);

	for (int i = 0; i < 20; i++)
	{
		coords[14 + i][0] = 5.0 + 2.0*(2.0*(float)rand()/RAND_MAX - 1.0);
		coords[14 + i][1] = 0.0 + 2.0*(2.0*(float)rand()/RAND_MAX - 1.0);
	}

	for (auto coord:coords)
		points.push_back(Point(coord[0], coord[1]));

	pointStates.resize(points.size());
}

/*
 * generate two open curves (semicircles)
 */
void TestReconstruct2D::generateOpenCurves()
{
	int i, j;
	const int SAMPLE_COUNT = 10;

	for (j = 0; j < 2; j++)
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			double t = (double)i/SAMPLE_COUNT*PI + (1.5 + j)*PI;
			Point p(sin(t), 4*j + cos(t));
			points.push_back(p);
		}
}

/*
 * return radius for circle through point p with normalized normal n and point q
 * q can be mirrored on the line through n, therefore the radius is the circumradius of the triangle pqq'
 */
float radiusForCircleThrough2PointsandNormal(Point p, Point n, Point q)
{
	float a = p.distance(q);
	Point n2(-n[1], n[0]);
	Point pq = p - q;
	float dist = abs(pq*n2);
	float b = 2*dist;

	if (b == 0)
		return 0.5*sqrt(pq.squared_length());	// distance pq = diameter

	float e = (2*a + b)*SQR(b)*(2*a - b);	// r=abc/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) -> isosceles a=b

	if (e <= 0.0)
		return numeric_limits<float>::max();	// triangle points are collinear, infinite radius

	float d = sqrt(e);
	return abs(SQR(a)*b/d);	// circumradius of triangle adapted to isosceles version
}

/*
 * return radius for circle through point p with normalized normal n and point q
 */
float radiusForCircleThrough2PointsandNormal2(Point p, Point n, Point q)
{
	// circle center c=p+t*n, |pc|=|qc|
	double px = p[0];
	double py = p[1];
	double qx = q[0];
	double qy = q[1];
	double nx = n[0];
	double ny = n[1];
	double vx = qx - px;
	double vy = qy - py;

	return abs((SQR(vx) + SQR(vy))/(2*(vx*nx + vy*ny)));
}

/*
 * compute LFS values for curve points
 */
void computeLFSForCurve(vector<Point> &curvePoints, vector<float> &lfsCurvePoints, vector<Point> &mPoints, bool isClosed)
{
	int i, j;

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		float minR = numeric_limits<float>::max();

		for (j = 0; j < (int)curvePoints.size(); j++)
			if (i != j)
			{
				// at point p, determine radius r for maximum empty circle with center through normal n (one neighbor q on its boundary)
				Point curr2P = curvePoints[j];
				float radius = radiusForCircleThrough2PointsandNormal(currP, normal, curr2P);
				radius = radiusForCircleThrough2PointsandNormal2(currP, normal, curr2P);

				if (radius < minR)
				{
					minR = radius;
					float direction = (normal*(curr2P - currP) < 0) ? -1.0 : 1.0;
					mPoints[i] = currP + normal*(radius*direction);
				}
			}
	}

	// compute lfs from nearest medial axis points
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points;

	ann_points = annAllocPts(mPoints.size(), 2);

	for(i = 0; i < (int)mPoints.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = mPoints[i][0];
		p[1] = mPoints[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, mPoints.size(), 2);
	ANNpointArray search_point = annAllocPts(1, 2);
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray distances = new ANNdist[1];

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		// get nearest neighbor in medial axis
		search_point[0][0] = curvePoints[i][0];
		search_point[0][1] = curvePoints[i][1];
		kdTree->annkSearch(search_point[0], 1, nnIdx, distances);
		lfsCurvePoints[i] = sqrt(distances[0]);
	}

	delete nnIdx;
	delete distances;
	annDeallocPts(ann_points);
}

/*
 * return distance of p0 from line p1-p2
 */
float distancePFromLine(Point p0, Point p1, Point p2)
{
	Point normal(p2 - p1);
	swap(normal[0], normal[1]);
	normal[0] = -normal[0];
	normal.normalize();
	return abs(normal*(p0 - p1));
}

/*
 * generate point set by sampling with condition < epsMax, max error to original curve and optionally manifold condition
 */
void TestReconstruct2D::sampleCurveByEps(float maxEps, float error, bool isClosed, float perturb)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
//	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	i = 0;

	do
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		int prevI = i;
		prevP = curvePoints[prevI];
		Point newP = prevP;

		if (perturb != 0.0)
		{
			// perturb p by up to perturb*lfs(p) along its normal in any direction
			float random = (float)rand()/RAND_MAX*2 - 1.0;
			newP = newP + normal*(random*perturb*lfsCurvePoints[i]);
		}

		points.push_back(newP);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];

			j = prevI + 1;
			isConforming = true;

			// test eps and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for epsilon condition (a sample within dist/lfs < maxEps)
					float lfs = lfsCurvePoints[j];

					float dist = currP.distance(prevP);

					if (dist/lfs < maxEps)
						isConforming = true;
					else
					{
						dist = currP.distance(p);
						isConforming = (dist/lfs < maxEps);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

//	determinePointStates(maxEps);
	computeScale();
}

/*
 * generate point set by sampling with condition < maxReach, max error to original curve
 */
void TestReconstruct2D::sampleCurveByReach(float maxReach, float error, bool isClosed, float perturb)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
//	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	i = 0;

	do
	{
		// compute normal
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];

		int prevI = i;
		prevP = curvePoints[prevI];
		Point newP = prevP;

		if (perturb != 0.0)
		{
			// perturb p by up to perturb*lfs(p) along its normal in any direction
			float random = (float)rand()/RAND_MAX*2 - 1.0;
			newP = newP + normal*(random*perturb*lfsCurvePoints[i]);
		}

		points.push_back(newP);
		noise.push_back(perturb*lfsCurvePoints[i]);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];
			float reach = numeric_limits<float>::max();

			for (j = prevI; j <= i; j++)
			{
				float lfs = lfsCurvePoints[j];

				if (lfs < reach)
					reach = lfs;
			}

			j = prevI + 1;
			isConforming = true;

			// test reach and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for reach condition (a sample within dist/reach < maxReach)
					// use prev/next points so that discrete interval of curve points does not impact
					float dist = curvePoints[j + 1].distance(prevP);

					if (dist/reach < maxReach)
						isConforming = true;
					else
					{
						dist = curvePoints[j - 1].distance(p);
						isConforming = (dist/reach < maxReach);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

	if (!isClosed)
		points.push_back(curvePoints[curvePoints.size() - 1]);

//	determinePointStates(maxReach);
	computeScale();
}

TestReconstruct2D::TestReconstruct2D()
{
	viewWidth = 1;
	viewHeight = 1;
	scale = 1.0;
	offset[0] = 0.0;
	offset[1] = 0.0;
	lfsMin = 0.0;
	lfsMax = 0.0;
}

TestReconstruct2D::~TestReconstruct2D()
{
}

void TestReconstruct2D::reshape(int width, int height)
{
	this->viewWidth = width;
	this->viewHeight = height;
}

/*
 * transform point into unit square
 */
Point TestReconstruct2D::transform(Point p)
{
	int i;
	Point tp;

	for (i = 0; i < 2; i++)
		tp[i] = 0.05 + (p[i] - offset[i])/scale*0.9;	// add 5% border

	return tp;
}

void TestReconstruct2D::drawPoints(vector<Point> &points)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawPoints(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
		if (pClasses[i] == FITTED)
		{
			Point tp = transform(points[i]);
			glVertex2f(tp[0], tp[1]);
		}

	glEnd();
}

void TestReconstruct2D::drawPointsWithClasses(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		if (pClasses[i] == FITTED)
			glColor3f(0.0, 0.0, 0.0);	// black
		else
			glColor3f(0.5, 0.5, 0.5);	// gray

		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap)
{
	int i;

	glBegin(GL_LINES);

	for (auto edgeItem:edgeMap)
	{
		for (i = 0; i < 2; i++)
		{
			int index = ((i == 0) ? edgeItem.first.first : edgeItem.first.second);
			Point tp = transform(points[index]);
			glVertex2f(tp[0], tp[1]);
		}
	}

	glEnd();
}

void TestReconstruct2D::drawArc(Circle &c, pair<int, int> arc)
{
	int i;

	glBegin(GL_LINES);
	Point oldTP;

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	for (i = arc.first; i < arc.second; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point tp = transform(p);

		if (i != arc.first)
		{
			glVertex2f(oldTP[0], oldTP[1]);
			glVertex2f(tp[0], tp[1]);
		}

		oldTP = tp;
	}

	glEnd();
}

void TestReconstruct2D::drawCircle(Circle &c)
{
	pair<int, int> arc(0, 360);
	drawArc(c, arc);
}

void TestReconstruct2D::drawCover(Circle &c, pair<int, int> arc)
{
	int i, j;

	glBegin(GL_QUADS);
	Point oldTP[2];

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	for (i = arc.first; i < arc.second; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point n(cos(degInRad)*c.variance, sin(degInRad)*c.variance);
		Point pp[2] = { p + n, p - n };
		Point tp[2];

		for (j = 0; j < 2; j++)
			tp[j] = transform(pp[j]);

		if (i != arc.first)
		{
			glVertex2f(oldTP[0][0], oldTP[0][1]);
			glVertex2f(tp[0][0], tp[0][1]);
			glVertex2f(tp[1][0], tp[1][1]);
			glVertex2f(oldTP[1][0], oldTP[1][1]);
		}

		for (j = 0; j < 2; j++)
			oldTP[j] = tp[j];
	}

	glEnd();
}

void TestReconstruct2D::drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arcs)
{
	int i;

	if (arcs.size() > 0)
		for (i = 0; i < (int)circles.size(); i++)
			drawArc(circles[i], arcs[i]);
}

void TestReconstruct2D::drawCovers(vector<Circle> &circles, vector<pair<int, int> > &arcs)
{
	int i;

	if (arcs.size() > 0)
		for (i = 0; i < (int)circles.size(); i++)
			drawCover(circles[i], arcs[i]);
}

void drawNormal(Point tp, Point n)
{
	glVertex2f(tp[0], tp[1]);
	Point endP = tp + n*0.05;
	glVertex2f(endP[0], endP[1]);

	// draw arrow hooks
	Point hookVec;
	float cs = cos(0.75*PI);
	float sn = sin(0.75*PI);
	hookVec[0] = 0.01*(cs*n[0] - sn*n[1]);
	hookVec[1] = 0.01*(sn*n[0] + cs*n[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] + hookVec[0], endP[1] + hookVec[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] - hookVec[1], endP[1] + hookVec[0]);
}

void TestReconstruct2D::drawNormals(vector<Point> &points, vector<PointsEnum> &pClasses)
{
	int i;

	glBegin(GL_LINES);

	for (i = 0; i < (int)points.size(); i++)
		if (pClasses[i] == FITTED)
		{
			// draw arrow line
			Point tp = transform(points[i]);
			drawNormal(tp, normals[i]);
		}

	glEnd();
}

void TestReconstruct2D::drawLabels(vector<Point> &points)
{
	int i, j;

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glRasterPos2f(tp[0] + 0.01, tp[1] - 0.00);
		char str[5];
		sprintf(str, "%d", i);

		for (j = 0; j < (int)strlen(str); j++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, str[j]);
	}
}

void TestReconstruct2D::drawNoiseExtent(vector<Point> &points, vector<float> &noise)
{
	int i;

	for (i = 0; i < (int)points.size(); i++)
	{
		Circle circle(points[i][0], points[i][1], noise[i]);
		drawCircle(circle);
	}
}

const int DRAW_POINT = 1;
const int DRAW_PROJPOINT = 2;
const int DRAW_CURVEPOINT = 4;
const int DRAW_MEDIALPOINT = 8;
const int DRAW_NORMAL = 16;
const int DRAW_LABEL = 32;
const int DRAW_EDGE = 64;
const int DRAW_ARC = 128;
const int DRAW_COVER = 256;
const int DRAW_POINTCLASS = 512;
const int DRAW_DENOISEDPOINT = 1024;
const int DRAW_NOISEEXTENT = 2048;

void TestReconstruct2D::draw(int items, float pointSize)
{
	// clear screen to white
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	if (items & DRAW_COVER)
	{
		// draw covers
		glEnable (GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4f(0.0, 1.0, 0.0, 0.15);	// green
		glLineWidth(1.0);
		drawCovers(circles, arcs);
		glDisable(GL_BLEND);
	}

	if (items & DRAW_ARC)
	{
		// draw arcs
		glColor3f(0.0, 0.0, 1.0);	// blue
		glLineWidth(1.0);
		drawArcs(circles, arcs);
	}

	if (items & DRAW_EDGE)
	{
		// draw lines
		glColor3f(1.0, 0.0, 0.0);	// red
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
		drawEdges(projPoints, edgeMap);
	}

	if (items & DRAW_POINT)
	{
		// draw points
		glColor3f(0.0, 0.0, 0.0);	// black
		glPointSize(pointSize);
		glEnable(GL_POINT_SMOOTH);
		drawPoints(points);
	}

	if (items & DRAW_POINTCLASS)
	{
		// draw points
		glPointSize(pointSize);
		glEnable(GL_POINT_SMOOTH);
		drawPointsWithClasses(points, pClasses);
	}

	if (items & DRAW_PROJPOINT)
	{
		// draw projected points
		glPointSize(pointSize);
		glColor3f(0.0, 0.0, 1.0);	// blue
		drawPoints(projPoints, pClasses);
	}

	if (items & DRAW_DENOISEDPOINT)
	{
		// draw projected points
		glPointSize(pointSize);
		glColor3f(0.0, 0.0, 1.0);	// blue
		drawPoints(denoisedPoints, pClasses);
	}

	if (items & DRAW_NORMAL)
	{
		// draw normals
		glColor3f(0.0, 0.0, 0.0);	// black
		glLineWidth(1.0);
		drawNormals(projPoints, pClasses);
	}

	if (items & DRAW_LABEL)
	{
		// draw labels
		glColor3f(0.0, 0.0, 0.0);	// black
		drawLabels(points);
	}

	if (items & DRAW_NOISEEXTENT)
	{
		glColor3f(0.5, 0.5, 0.5);	// grey
		drawNoiseExtent(points, noise);
	}
}

/*
 * write screen to PNG file
 */
void writePNG(int width, int height, string fileName)
{
    int x, y, npixels = width*height;
	GLfloat* pixels = new GLfloat[npixels*3];
	glReadPixels(0.0, 0.0, width, height, GL_RGB, GL_FLOAT, pixels);
	png::image<png::rgb_pixel> image(width, height);
	double R, G, B;

	for (y = height - 1; y >= 0; y--)
		for (x = 0; x < width; x++)
		{
			R = *pixels++;
			G = *pixels++;
			B = *pixels++;
			image[y][x] = png::rgb_pixel(255*R, 255*G, 255*B);     // set pixel to position (x, y)
		}

	image.write(fileName);
}

/*
 * draw and write to PNG file
 */
void TestReconstruct2D::draw(string filename, int items, float pointSize)
{
	draw(items, pointSize);
	writePNG(viewWidth, viewHeight, filename);
}

static TestReconstruct2D *instance = NULL;

///////////////////////////// GL FUNCTIONS /////////////////////////////////

// GLUT idle function
void idle()
{
}

// GLUT display function
void display(string filename, int items, float pointSize)
{
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);
	glViewport(0, 0, viewWidth, viewHeight);
	instance->draw(filename, items, pointSize);

	// restore the stored viewport dimensions
	glViewport(vp[0], vp[1], vp[2], vp[3]);
	glutSwapBuffers();
}

// GLUT reshape function
void reshape(int newWidth, int newHeight)
{
	if (newWidth == 0)
		newWidth = 1;

	if (newHeight == 0)
		newHeight = 1;

	viewWidth = newWidth;
	viewHeight = newHeight;

	glViewport(0, 0, viewWidth, viewHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// update application
	instance->reshape(viewWidth, viewHeight);
}

void do_realdata(TestReconstruct2D *instance, string pointsetname)
{
	instance->loadPointSet("./data/" + pointsetname + ".txt");
//	instance->loadPointSet("/home/stef/web/amit2d/" + pointsetname + ".txt");
	instance->invertY();
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
}

void fig_realdata(bool edge, string pointsetname, string filename, float pointsize)
{
	instance = new TestReconstruct2D();
	do_realdata(instance, pointsetname);
	reshape(viewWidth, viewHeight);
	display(filename, DRAW_POINT | (edge ? DRAW_EDGE : 0), pointsize);
}

void fig_realdata_a()
{
	fig_realdata(false, "Keyboard0", "fig_realdata_a.png", 1.0);
}

void fig_realdata_b()
{
	fig_realdata(false, "Monitor0", "fig_realdata_b.png", 1.0);
}

void fig_realdata_c()
{
	fig_realdata(false, "Cup0", "fig_realdata_c.png", 1.0);
}

void fig_realdata_d()
{
	fig_realdata(true, "Keyboard0", "fig_realdata_d.png", 1.0);
}

void fig_realdata_e()
{
	fig_realdata(true, "Monitor0", "fig_realdata_e.png", 1.0);
}

void fig_realdata_f()
{
	fig_realdata(true, "Cup0", "fig_realdata_f.png", 1.0);
}

void fig_mouse_a()
{
	fig_realdata(false, "Mouse0", "fig_mouse_a.png", 5.0);
}

void fig_mouse_b()
{
	fig_realdata(true, "Mouse0", "fig_mouse_b.png", 5.0);
}

void do_realdata_rect_otr(TestReconstruct2D *instance, int iterations)
{
	instance->generateHighNoiseClosed();
	instance->computeScale();
	instance->reconstruct_otr(iterations);
}

void fig_realdata_rect_otr(bool edge, string filename, float pointsize, int iterations)
{
	instance = new TestReconstruct2D();
	do_realdata_rect_otr(instance, iterations);
	reshape(viewWidth, viewHeight);
	display(filename, DRAW_POINT | (edge ? DRAW_EDGE : 0), pointsize);
}

void do_realdata_bunny_otr(TestReconstruct2D *instance, int iterations)
{
	string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";
	srand(1);
	instance->generateCurveFromBezierString(bunny);
	instance->sampleCurveByReach(0.43, 0.0, true, 0.0);	// rho=0.43 -> eps>0.3
	instance->computeScale();
	instance->reconstruct_otr(iterations);
}

void fig_realdata_bunny_otr(bool edge, string filename, float pointsize, int iterations)
{
	instance = new TestReconstruct2D();
	do_realdata_bunny_otr(instance, iterations);
	reshape(viewWidth, viewHeight);
	display(filename, DRAW_POINT | (edge ? DRAW_EDGE : 0), pointsize);
}

void fig_rect_otr_a()
{
	fig_realdata_rect_otr(true, "fig_rect_otr_a.png", 5.0, 25);
}

void fig_rect_otr_b()
{
	fig_realdata_rect_otr(true, "fig_rect_otr_b.png", 5.0, 30);
}

void fig_rect_otr_c()
{
	fig_realdata_rect_otr(true, "fig_rect_otr_c.png", 5.0, 35);
}

void fig_bunny_otr_a()
{
	fig_realdata_bunny_otr(true, "fig_bunny_otr_a.png", 5.0, 90);
}

void fig_bunny_otr_b()
{
	fig_realdata_bunny_otr(true, "fig_bunny_otr_b.png", 5.0, 100);
}

void fig_bunny_otr_c()
{
	fig_realdata_bunny_otr(true, "fig_bunny_otr_c.png", 5.0, 110);
}

void fig_boundary_a()
{
	instance = new TestReconstruct2D();
	instance->generateNoisyLeafFeature(0.05, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_boundary_a.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_boundary_b()
{
	instance = new TestReconstruct2D();
	instance->generateNoisyLeafLine(0.1, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_boundary_b.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_blend_a()
{
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(0.5, 1);
	instance->computeScale();
	instance->reconstruct(0, -1);
	reshape(viewWidth, viewHeight);
	display("fig_blend_a.png", DRAW_POINT, 10.0);
}

void fig_blend_b()
{
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(0.5, 1);
	instance->computeScale();
	instance->reconstruct(0, -1);
	reshape(viewWidth, viewHeight);
	display("fig_blend_b.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void do_blend_c()
{
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(0.5, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
}

void fig_blend_c()
{
	do_blend_c();
	display("fig_blend_c.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_iter(int iter, string filename)
{
	instance = new TestReconstruct2D();
	instance->generateHighNoiseClosed();
	instance->computeScale();
	instance->reconstruct(0, iter);
	reshape(viewWidth, viewHeight);
	display(filename, DRAW_POINTCLASS | DRAW_EDGE | DRAW_ARC | DRAW_COVER, 20.0);
}

void fig_iter_a()
{
	fig_iter(0, "fig_iter_a.png");
}

void fig_iter_b()
{
	fig_iter(1, "fig_iter_b.png");
}

void fig_iter_c()
{
	fig_iter(3, "fig_iter_c.png");
}

void fig_iter_d()
{
	fig_iter(15, "fig_iter_d.png");
}

void fig_iter_e()
{
	fig_iter(22, "fig_iter_e.png");
}

void do_basic_a()
{
	instance = new TestReconstruct2D();
//	instance->generatePointsFrom2Circles(0.0, 1);
	instance->generateHoleOutliers();
	instance->computeScale();
//	instance->reconstruct(MODE_BLEND, -1);
	instance->reconstruct(0, -1);
	reshape(viewWidth, viewHeight);
}

void fig_basic_a()
{
	do_basic_a();
	display("fig_basic_a.png", DRAW_POINT | DRAW_EDGE, 20.0);
}

void fig_basic_b()
{
	instance = new TestReconstruct2D();
	instance->generateHighNoiseClosed();
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_basic_b.png", DRAW_POINT | DRAW_EDGE, 20.0);
}

void do_basic_c()
{
	instance = new TestReconstruct2D();
	instance->generateNoisySpiral(0.25, 1);	// noise 0.25 is 50% of distance between rotated boundaries
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
}

void fig_basic_c()
{
	do_basic_c();
	display("fig_basic_c.png", DRAW_POINT | DRAW_EDGE, 20.0);
}

void do_basic_d()
{
	instance = new TestReconstruct2D();
	instance->generateOpenCurves();
	instance->computeScale();
	instance->reconstruct(0, -1);
	reshape(viewWidth, viewHeight);
}

void fig_basic_d()
{
	do_basic_d();
	display("fig_basic_d.png", DRAW_POINT | DRAW_EDGE, 20.0);
}

void do_circle(TestReconstruct2D *instance, float noise)
{
	instance->generateVaryingNoisyCircle(noise, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
}

void fig_circle(float noise, string subfig, bool extents)
{
	instance = new TestReconstruct2D();
	do_circle(instance, noise);
	reshape(viewWidth, viewHeight);
	display("fig_circle_" + subfig + ".png", DRAW_POINT | DRAW_EDGE | (extents ? (DRAW_ARC | DRAW_COVER) : 0), 10.0);
}

void fig_circle_a()
{
	fig_circle(0.1, "a", false);
}

void fig_circle_b()
{
	fig_circle(0.25, "b", false);
}

void fig_circle_c()
{
	fig_circle(0.5, "c", false);
}

void fig_circle_d()
{
	fig_circle(0.75, "d", false);
}

void fig_circle_e()
{
	fig_circle(1.0, "e", false);
}

void fig_circle_f()
{
	fig_circle(1.25, "f", false);
}

void fig_circle_extent_a()
{
	fig_circle(0.1, "a2", true);
}

void fig_circle_extent_b()
{
	fig_circle(0.5, "b2", true);
}

void fig_circle_extent_c()
{
	fig_circle(1.0, "c2", true);
}

void do_lfs(TestReconstruct2D *instance, float noise)
{
	string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";
	srand(1);
	instance->generateCurveFromBezierString(bunny);
	instance->sampleCurveByReach(0.43, 0.0, true, noise);	// rho=0.43 -> eps>0.3
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
}

void fig_lfs(float noise, string subfig)
{
	instance = new TestReconstruct2D();
	do_lfs(instance, noise);
	reshape(viewWidth, viewHeight);
	display("fig_lfs_" + subfig + ".png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_lfs_a()
{
	fig_lfs(0.0, "a");
}

void fig_lfs_b()
{
	fig_lfs(0.1, "b");
}

void fig_lfs_c()
{
	fig_lfs(1.0/3.0, "c");
}

void fig_lfs_d()
{
	fig_lfs(0.5, "d");
}

void fig_limit_a()
{
	instance = new TestReconstruct2D();
	instance->generateNoisySharpCorner(1.0/3, 0.05, 1);	// test cases with seeds 1-30
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_limit_a.png", DRAW_POINT | DRAW_EDGE, 5.0);
}

void fig_limit_b()
{
	instance = new TestReconstruct2D();
	instance->generateNoisyTCrossing(0.05, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_limit_b.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_feature(float noise, string subfig)
{
	instance = new TestReconstruct2D();
	instance->generateNoisyFeatures(noise, 1);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display("fig_feature_" + subfig + ".png", DRAW_POINT | DRAW_EDGE, 5.0);
}

void fig_feature_a()
{
	fig_feature(0.0, "a");
}

void fig_feature_b()
{
	fig_feature(0.05, "b");
}

void fig_teaser_a()
{
	do_blend_c();
	display("fig_teaser_a.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_teaser_b()
{
	do_basic_a();
	display("fig_teaser_b.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void fig_teaser_c()
{
	fig_realdata(true, "Monitor0", "fig_teaser_c.png", 1.0);
}

void fig_teaser_d()
{
	do_basic_c();
	display("fig_teaser_d.png", DRAW_POINT | DRAW_EDGE, 10.0);
}

void tab_circle_reconstruct(stringstream &sout, float noise)
{
	int i, j;
	instance = new TestReconstruct2D();
	instance->generateVaryingNoisyCircle(noise, 1);

	// compute mean error of noisy points
	vector<Point> points = instance->getPoints();
	float sumDist = 0.0;
	float sqrSumDist = 0.0;
	float inputMaxErr = 0.0;

	for (i = 0; i < (int)points.size(); i++)
	{
		float dist = abs(sqrt(SQR(points[i][0]) + SQR(points[i][1])) - 1.0);
		sumDist += dist;
		sqrSumDist += SQR(dist);

		if (dist > inputMaxErr)
			inputMaxErr = dist;
	}

	float inputMeanErr = sumDist/points.size();
	float inputRMSErr = sqrt(sqrSumDist/points.size());

	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);

	// compute max and root mean square error of reconstructed polygon (sampled at edges with constant density)
	points = instance->getProjPoints();
	map<pair<int, int>, EdgeEnum> edgeMap = instance->getEdges();
	sumDist = 0.0;
	sqrSumDist = 0.0;
	float outputMaxErr = 0.0;
	int count = 0;

	for (auto edgeItem:edgeMap)
	{
		Point p0 = points[edgeItem.first.first], p1 = points[edgeItem.first.second];
		Point v = p1 - p0;
		float edgeLen = sqrt(v.squared_length());

		for (j = 0; j < (int)(edgeLen*100.0); j++)
		{
			float len = 0.01/edgeLen;
			Point p = p0 + v*j*len;
			float dist = abs(sqrt(SQR(p[0]) + SQR(p[1])) - 1.0);
			sumDist += dist;
			sqrSumDist += SQR(dist);

			if (dist > outputMaxErr)
				outputMaxErr = dist;

			count++;
		}
	}

	float outputMeanErr = sumDist/count;
	float outputRMSErr = sqrt(sqrSumDist/count);

	sout << noise << "\t| " << inputMaxErr << "\t| " << inputMeanErr << "\t| " << inputRMSErr << "\t| " << outputMaxErr << "\t| " << outputMeanErr << "\t| " << outputRMSErr << endl;
}

void tab_circle()
{
	int i;
	stringstream sout;
	sout.setf(ios::fixed, ios::floatfield);
	sout.precision(3);
	float noise[] = { 0.1, 0.25, 0.5, 0.75, 1.0 };
	sout << "Table 2:" << endl;
	sout << "Noise \t| maxI \t| meanI\t| rmsI \t| maxO \t| meanO\t| rmsO" << endl;
	sout << "=======================================================" << endl;

	for (i = 0; i < 5; i++)
		tab_circle_reconstruct(sout, noise[i]);

	// output table
	cout << sout.str();
}

void printUsage(char *argv[])
{
	cout << "Usage: " << argv[0] << " fig1a|fig1b|fig1c|fig1d|fig4a|fig4b|fig4c|" <<
		"fig4d|fig4e|fig5a|fig5b|fig5c|fig6a|fig6b|fig6c|fig6d|fig6e|fig6f|fig6g|fig6h|fig7a|fig7b|fig7c|" <<
		"fig8a|fig8b|fig8c|fig8d|fig8e|fig9a|fig9b|fig9c|fig9d|"
		"fig10a|fig10b|fig10c|fig11a|fig11b|fig11c|fig11d|fig11e|fig11f|fig12a|fig12b|fig12c|fig12d|fig13a|fig13b|fig13c|fig14a|fig14b|fig15a|fig15b|"
		"tab1|tab2" << endl;
}

/*
 * construct line touching the top or bottom of both discs (p0, r0) and (p1, r1) as point tp + normalized vector tv
 */
bool constructLine(Point p0, Point p1, float r0, float r1, Point &tp, Point &tv, bool isTop0, bool isTop1)
{
	float centerDist = p0.distance(p1);

	// ensure r0 >= r1
	if (r0 < r1)
	{
		swap(p0, p1);
		swap(r0, r1);
		swap(isTop0, isTop1);
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
		float root = sqrt(denom - SQR(r0));
		Point tx0(p0[0] + (SQR(r0)*dx + r0*dy*root)/denom, p0[1] + (SQR(r0)*dy - r0*dx*root)/denom);
		Point tx1(p0[0] + (SQR(r0)*dx - r0*dy*root)/denom, p0[1] + (SQR(r0)*dy + r0*dx*root)/denom);
		Point tx = (isTop0 ^ (tx0[1] < tx1[1])) ? tx0 : tx1;
		tv = tx - tp;
	}

	tv.normalize();

	return true;
}

/*
 * linear least squares fit line to points
 */
void determineLinearLSFit(vector<Point> &points, vector<int> &indices, float &a, float &b)
{
	int i, count = indices.size();

	float meanx = 0.0, meany = 0.0, ssxx = 0.0, ssxy = 0.0;

	for (i = 0; i < count; i++)
	{
		meanx += points[indices[i]][0];
		meany += points[indices[i]][1];
	}

	meanx /= count;
	meany /= count;

	for (i = 0; i < count; i++)
	{
		ssxx += SQR(points[indices[i]][0] - meanx);
		ssxy += (points[indices[i]][0] - meanx)*(points[indices[i]][1] - meany);
	}

	b = ssxy/ssxx;
	a = meany - b*meanx;
}

void handleKeypress(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 27:
			exit(0);
	}
}

/*
 * intersects lines p0+s*v0 and p1+t*v1 and returns false if parallel, else intersection point x
 */
bool intersectLinesTest(Point p0, Point v0, Point p1, Point v1, Point &x)
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

void do_fig_comparison(string filename, string figname, float pointSize)
{
	instance = new TestReconstruct2D();
//	instance->loadPointSet("/home/stef/web/amit2d/" + filename + ".txt");
	instance->loadPointSet("./data/" + filename + ".txt");
	instance->invertY();
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
	reshape(viewWidth, viewHeight);
	display(figname + ".png", DRAW_POINT | DRAW_EDGE, pointSize);
}

void fig_comparison_a()
{
	do_fig_comparison("apple_2percent_noise", "fig_comparison_a", 7.0);
}

void fig_comparison_b()
{
	do_fig_comparison("butterfly_2percent_noise", "fig_comparison_b", 7.0);
}

void fig_comparison_c()
{
	do_fig_comparison("crab_2percent_noise", "fig_comparison_c", 7.0);
}

void fig_comparison_d()
{
	do_fig_comparison("dolphin_2percent_noise", "fig_comparison_d", 7.0);
}

void fig_comparison2_a()
{
	do_fig_comparison("fish", "fig_comparison2_a", 5.0);
}

void fig_comparison2_b()
{
	do_fig_comparison("bot", "fig_comparison2_b", 5.0);
}

void do_highnoise()
{
	string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";
	instance->generateCurveFromBezierString(bunny); instance->sampleCurveByReach(0.01, 0.0, true, 1.0/3);
	instance->computeScale();
	instance->reconstruct(MODE_BLEND, -1);
}

void fig_highnoise()
{
	srand(1);
	instance = new TestReconstruct2D();
	do_highnoise();
	reshape(viewWidth, viewHeight);
	display("fig_highnoise.png", DRAW_POINT | DRAW_EDGE, 5);
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_2                                          Point2;
typedef K::Segment_2                                        Segment;
typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr_2;

void TestReconstruct2D::reconstruct_otr(int iterations)
{
	int i;
	vector<Point2> _points;

	for (auto p:points)
		_points.push_back(Point2(p.x(), p.y()));

	Otr_2 otr2(_points);
	otr2.run(iterations);

	vector<Point2> isolated_points;
	vector<Segment> segments;
	otr2.list_output(back_inserter(isolated_points), back_inserter(segments));
	vector<Point2>::iterator pit;

	points.clear();
	projPoints.clear();

	// add segments to our structure
	for (auto p:_points)
		points.push_back(Point(p[0], p[1]));

	// collect new points into projPoints
	for (auto p:isolated_points)
		projPoints.push_back(Point(p[0], p[1]));

	map<Point2, int> pointMap;

	for (auto s:segments)
		for (i = 0; i < 2; i++)
			pointMap[s[i]] = -1;

	for (auto entry:pointMap)
	{
		Point2 p = entry.first;
		pointMap[p] = projPoints.size();
		projPoints.push_back(Point(p[0], p[1]));
	}

	for (auto s:segments)
	{
		int p0 = pointMap[s[0]];
		int p1 = pointMap[s[1]];
		edgeMap[pair<int, int>(p0, p1)] = BIJECTIVE;
	}

	pClasses.resize(points.size());
}

void tab_runtime()
{
	int i, outputSize, iterations, handledFits, handledPoints, squaredFits;
	float runtime;
	string name[] = { "Circle 0.1r", "Circle 0.25r", "Circle 0.5r", "Circle 0.75r", "Circle 1r",
			"0\t", "0.1lfs\t", "1/3lfs\t", "0.5lfs\t", "Keyboard", "Monitor\t", "Cup\t", "Mouse\t",
			"Apple\t", "Butterfly\t", "Crab\t", "Dolphin\t", "Fish\t", "Bottle\t", "Bunny\t"
	};
	stringstream sout;
	sout.setf(ios::fixed, ios::floatfield);
	sout.precision(3);

	sout << "Table 1:" << endl;
	sout << "Name\t| Input\t| Output\t| Iterations\t| Operations\t| Complexity\t| Runtime" << endl;
	sout << "===================================================================================" << endl;

	for (i = 0; i < 20; i++)
	{
		instance = new TestReconstruct2D();

		if (i == 0)
			do_circle(instance, 0.1);
		else
		if (i == 1)
			do_circle(instance, 0.25);
		else
		if (i == 2)
			do_circle(instance, 0.5);
		else
		if (i == 3)
			do_circle(instance, 0.75);
		else
		if (i == 4)
			do_circle(instance, 1.0);
		else
		if (i == 5)
			do_lfs(instance, 0.0);
		else
		if (i == 6)
			do_lfs(instance, 0.1);
		else
		if (i == 7)
			do_lfs(instance, 1.0/3.0);
		else
		if (i == 8)
			do_lfs(instance, 0.5);
		else
		if (i == 9)
			do_realdata(instance, "Keyboard0");
		else
		if (i == 10)
			do_realdata(instance, "Monitor0");
		else
		if (i == 11)
			do_realdata(instance, "Cup0");
		else
		if (i == 12)
			do_realdata(instance, "Mouse0");
		else
		if (i == 13)
			do_realdata(instance, "apple_2percent_noise");
		else
		if (i == 14)
			do_realdata(instance, "butterfly_2percent_noise");
		else
		if (i == 15)
			do_realdata(instance, "crab_2percent_noise");
		else
		if (i == 16)
			do_realdata(instance, "dolphin_2percent_noise");
		else
		if (i == 17)
			do_realdata(instance, "fish");
		else
		if (i == 18)
			do_realdata(instance, "bot");
		else
		if (i == 19)
			do_highnoise();

		instance->getData(outputSize, iterations, handledFits, handledPoints, squaredFits, runtime);
		vector<Point> points = instance->getPoints();
		sout << name[i] << "\t" << points.size() << "\t" << outputSize << "\t|" << iterations << "\t| " << handledPoints << "\t| " << squaredFits << "\t| " << runtime << endl;
		delete instance;
	}

	// output table
	cout << sout.str();
}

void fig_all()
{
	fig_teaser_a();
	fig_teaser_b();
	fig_teaser_c();
	fig_teaser_d();
	fig_iter_a();
	fig_iter_b();
	fig_iter_c();
	fig_iter_d();
	fig_iter_e();
	fig_blend_a();
	fig_blend_b();
	fig_blend_c();
	fig_realdata_a();
	fig_realdata_b();
	fig_realdata_c();
	fig_realdata_d();
	fig_realdata_e();
	fig_realdata_f();
	fig_basic_a();
	fig_basic_b();
	fig_basic_c();
	fig_basic_d();
	fig_circle_a();
	fig_circle_b();
	fig_circle_c();
	fig_circle_d();
	fig_circle_e();
	fig_comparison_a();
	fig_comparison_b();
	fig_comparison_c();
	fig_comparison_d();
	fig_highnoise();
	fig_rect_otr_a();
	fig_rect_otr_b();
	fig_rect_otr_c();
	fig_bunny_otr_a();
	fig_bunny_otr_b();
	fig_bunny_otr_c();
	fig_lfs_a();
	fig_lfs_b();
	fig_lfs_c();
	fig_lfs_d();
	fig_circle_extent_a();
	fig_circle_extent_b();
	fig_circle_extent_c();
	fig_feature_a();
	fig_feature_b();
	fig_mouse_a();
	fig_mouse_b();
	fig_limit_a();
	fig_limit_b();
	tab_runtime();
	tab_circle();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowSize(viewWidth, viewHeight);
	glutCreateWindow("Reconstruct from 2D Points");

	if (argc <= 1)
		printUsage(argv);
	else
	{
		string a = argv[1];

		if (!a.compare("fig1a"))
			fig_teaser_a();
		else
		if (!a.compare("fig1b"))
			fig_teaser_b();
		else
		if (!a.compare("fig1c"))
			fig_teaser_c();
		else
		if (!a.compare("fig1d"))
			fig_teaser_d();
		else
		if (!a.compare("fig4a"))
			fig_iter_a();
		else
		if (!a.compare("fig4b"))
			fig_iter_b();
		else
		if (!a.compare("fig4c"))
			fig_iter_c();
		else
		if (!a.compare("fig4d"))
			fig_iter_d();
		else
		if (!a.compare("fig4e"))
			fig_iter_e();
		else
		if (!a.compare("fig5a"))
			fig_blend_a();
		else
		if (!a.compare("fig5b"))
			fig_blend_b();
		else
		if (!a.compare("fig5c"))
			fig_blend_c();
		else
		if (!a.compare("fig6a"))
			fig_realdata_a();
		else
		if (!a.compare("fig6b"))
			fig_realdata_b();
		else
		if (!a.compare("fig6c"))
			fig_realdata_c();
		else
		if (!a.compare("fig6d"))
			fig_realdata_d();
		else
		if (!a.compare("fig6e"))
			fig_realdata_e();
		else
		if (!a.compare("fig6f"))
			fig_realdata_f();
		else
		if (!a.compare("fig6g"))
			fig_mouse_a();
		else
		if (!a.compare("fig6h"))
			fig_mouse_b();
		else
		if (!a.compare("fig7a"))
			fig_basic_a();
		else
		if (!a.compare("fig7b"))
			fig_basic_b();
		else
		if (!a.compare("fig7c"))
			fig_basic_c();
		else
		if (!a.compare("fig7d"))
			fig_basic_d();
		else
		if (!a.compare("fig8a"))
			fig_circle_a();
		else
		if (!a.compare("fig8b"))
			fig_circle_b();
		else
		if (!a.compare("fig8c"))
			fig_circle_c();
		else
		if (!a.compare("fig8d"))
			fig_circle_d();
		else
		if (!a.compare("fig8e"))
			fig_circle_e();
		else
		if (!a.compare("fig9a"))
			fig_comparison_a();
		else
		if (!a.compare("fig9b"))
			fig_comparison_b();
		else
		if (!a.compare("fig9c"))
			fig_comparison_c();
		else
		if (!a.compare("fig9d"))
			fig_comparison_d();
		else
		if (!a.compare("fig10a"))
			fig_comparison2_a();
		else
		if (!a.compare("fig10b"))
			fig_comparison2_b();
		else
		if (!a.compare("fig10c"))
			fig_highnoise();
		else
		if (!a.compare("fig11a"))
			fig_rect_otr_a();
		else
		if (!a.compare("fig11b"))
			fig_rect_otr_b();
		else
		if (!a.compare("fig11c"))
			fig_rect_otr_c();
		else
		if (!a.compare("fig11d"))
			fig_bunny_otr_a();
		else
		if (!a.compare("fig11e"))
			fig_bunny_otr_b();
		else
		if (!a.compare("fig11f"))
			fig_bunny_otr_c();
		else
		if (!a.compare("fig12a"))
			fig_lfs_a();
		else
		if (!a.compare("fig12b"))
			fig_lfs_b();
		else
		if (!a.compare("fig12c"))
			fig_lfs_c();
		else
		if (!a.compare("fig12d"))
			fig_lfs_d();
		else
		if (!a.compare("fig13a"))
			fig_circle_extent_a();
		else
		if (!a.compare("fig13b"))
			fig_circle_extent_b();
		else
		if (!a.compare("fig13c"))
			fig_circle_extent_c();
		else
		if (!a.compare("fig14a"))
			fig_feature_a();
		else
		if (!a.compare("fig14b"))
			fig_feature_b();
		else
		if (!a.compare("fig15a"))
			fig_limit_a();
		else
		if (!a.compare("fig15b"))
			fig_limit_b();
		else
		if (!a.compare("tab1"))
			tab_runtime();
		else
		if (!a.compare("tab2"))
			tab_circle();
		else
		{
			cout << "No valid parameter chosen" << endl;
			printUsage(argv);
		}
	}

	return 0;
}
