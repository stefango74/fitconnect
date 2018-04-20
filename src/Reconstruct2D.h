//============================================================================
// Author      : Stefan Ohrhallinger
// Version     :
// Copyright   : GPL v3
// Description : Reconstruct a curve from noisy 2D points
//============================================================================

#ifndef RECONSTRUCT2D_H_
#define RECONSTRUCT2D_H_

#include <set>
#include <stack>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <ANN/ANN.h>
#include "Point.h"
#include "circle.h"
#include "PrecTimer.h"

#define PI 3.1415926

using namespace std;

typedef enum { CONFORMING, NONCONFORMING, NOISY, OUTLIER, SHARP, FITTED, NONFITTED } PointsEnum;
typedef enum { UNIJECTIVE, BIJECTIVE } EdgeEnum;

const int MODE_BLEND = 1;

class Reconstruct2D {
private:
	bool isClosed;
	int mode, maxIter;
	int outputSize, iterations, handledFitCount, handledPointCount, squaredFitCount;
	float runtime;
	vector<Point> points, projPoints, normals, outputPoints, denoisedPoints;
	vector<float> noise;
	vector<PointsEnum> pClasses;
	map<pair<int, int>, EdgeEnum> visEdgeMap;
	vector<Circle> circles;
	vector<pair<int, int> > arcs;
	void determineNeighbors(ANNkd_tree *kdTree, ANNpointArray ann_points, int kMax,
			vector<pair<int, int> > &manifoldNeighbors);

public:
	Reconstruct2D(const vector<Point> &p_points, int mode);
	Reconstruct2D(const vector<Point> &p_points, const vector<float> & p_noise, int mode);
	bool reconstructNoisy();
	void setMaxIter(int p_maxIter);
	map<pair<int, int>, EdgeEnum> getEdgeMap();
	vector<Point> getProjectedPoints();
	vector<Point> getOutputPoints();
	vector<Point> getDenoisedPoints();
	vector<Circle> getCircles();
	vector<pair<int, int> > getArcs();
	vector<PointsEnum> getPointClassification();
	vector<Point> getNormals();
	void getData(int &p_output, int &p_iter, int &p_fit, int &p_point, int &p_squared, float &p_runtime);
};

typedef struct Node *NodePtr;

struct Node
{
	int index;
	vector<pair<NodePtr, int> > neighbors;
};

#endif /* RECONSTRUCT2D_H_ */
