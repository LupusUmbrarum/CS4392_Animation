#pragma once
#include <vector>
#include <iostream>

 #define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#define PI 3.1415926535897

struct TableEntry;

// I was able to get the curve working by following this tutorial:
// https://www.habrador.com/tutorials/interpolation/1-catmull-rom-splines/
class Curve
{
public:
	// ttotal is the total amount of time it should take to traverse the loop,
	// tcurrent is the current time that we're at, and tinc is the value that 
	// tcurrent will be incremented by
	float t1 = 0.5f, t2 = 0.5f, ttotal = 10.0f, tcurrent = 0.0f, tinc = 0.0f;
	// adjt2 is the adjusted time for t2, it equals ttotal - t2, giving the actual time position of t2
	float adjt2;
	float tau = 0.5; // Coefficient for catmull-rom spline
	int num_points_per_segment = 200;
	int numPoints = 0;
	int tableSize = 0;
	int aircraftPointIndex = 0;
	bool aircraftIsMoving = false;
	bool calcCurve = false;

	Curve();
	~Curve();
	
	void init(bool calc = false);
	void calculate_curve();
	float getNextDistance(float distance);
	float getTotalArcLength();
	int getNextIndex(float distance);

	std::vector<glm::vec3> control_points_pos;
	std::vector<glm::vec3> curvePoints;
	std::vector<glm::quat> control_points_quaternion;
	std::vector<float> distanceBetweenPoints;
	std::vector<TableEntry> table;
	std::vector<glm::quat> quats;
	
	void createSplineCurve(int pos);
	void calculateDistances();
	void createTable();
	void setupOtherValues();
	int clampForSpline(int pos);
	float distance(glm::vec3 a, glm::vec3 b);
	float getVectorMagnitude(glm::vec3 a);
	glm::vec3 crossProduct(glm::vec3 a, glm::vec3 b);
	float angleBetweenPoints(glm::vec3 a, glm::vec3 b);
	glm::vec3 calcPoint(float t, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
	glm::mat4 quatToMat4(glm::quat q);
	void calcAllQuats();
	glm::quat slerp(glm::quat a, glm::quat b, float t);
private:
};

struct TableEntry
{
	int segment, pointIndex;
	// length is sum, distance is individual
	float u, length, distance;
	glm::vec3 point;

	TableEntry()
	{
		segment = 0;
		pointIndex = 0;
		u = 0.0f;
		length = 0.0f;
		distance = 0.0f;
		point = glm::vec3(0, 0, 0);
	}

	TableEntry(int segment, int pointIndex, float u, float length, float distance, glm::vec3 point)
	{
		this->segment = segment;
		this->pointIndex = pointIndex;
		this->u = u;
		this->length = length;
		this->distance = distance;
		this->point = point;
	}

	void print()
	{
		std::cout << "Segment: " << segment << ", Index: " << pointIndex << ", U: ";
		std::cout << u << ", Length: " << length << ", Distance " << distance << ", Point: " << point.x << "," << point.y << "," << point.z << std::endl;
	}
};