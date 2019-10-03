#pragma once
#include <vector>
#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#define PI 3.1415926535897

struct TableEntry;

class Curve
{
public:
	// ttotal is the total amount of time it should take to traverse the loop,
	// tcurrent is the current time that we're at, and tinc is the value that 
	// tcurrent will be incremented by
	float t1 = 0.5f, t2 = 0.5f, ttotal = 5.0f, tcurrent = 0.0f, tinc = 0.0f;
	// adjt2 is the adjusted time for t2, it equals ttotal - t2, giving the actual time position of t2
	float adjt2;
	bool aircraftIsMoving = false;
	bool calcCurve = false;
	const double alpha = 0.5f;

	int numPoints = 0;
	int tableSize = 0;
	int aircraftPointIndex = 0;

	int getNextPoint(float distance);
	void setupOtherValues();

	// I was able to get it working by following this tutorial:
	// https://www.habrador.com/tutorials/interpolation/1-catmull-rom-splines/
	const glm::mat4 M = glm::mat4
	(
		-.5, 1.5, -1.5, .5,
		1, -2.5, 2, -.5,
		-.5, 0, .5, 0,
		0, 1, 0, 0
	);

	Curve();
	~Curve();

	void init(bool calc = false);
	void calculate_curve(bool calc = false);
	void tick(float deltaTime);

	float tau = 0.5; // Coefficient for catmull-rom spline
	int num_points_per_segment = 200;

	std::vector<glm::vec3> control_points_pos;
	std::vector<glm::vec3> curve_points_pos;
	std::vector<float> distanceBetweenPoints;
	std::vector<TableEntry> table;

	void calculateDistances();
	void createTable();
	float ease(float t);

private:
	float distance(glm::vec3 a, glm::vec3 b);
	float distance(int* a, int* b);

	void createSplineCurve(int pos);
	glm::vec3 calcPoint(float t, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
	int clampForSpline(int pos);
};

struct TableEntry
{
	int segment, pointIndex;
	float u, length;
	glm::vec3 point;

	TableEntry()
	{
		segment = 0;
		pointIndex = 0;
		u = 0.0f;
		length = 0.0f;
		point = glm::vec3(0, 0, 0);
	}

	TableEntry(int segment, int pointIndex, float u, float length, glm::vec3 point)
	{
		this->segment = segment;
		this->pointIndex = pointIndex;
		this->u = u;
		this->length = length;
		this->point = point;
	}

	void print()
	{
		std::cout << "Segment: " << segment << ", Index: " << pointIndex << ", U: ";
		std::cout << u << ", Length: " << length << ", Point: " << point.x << "," << point.y << "," << point.z << std::endl;
	}
};