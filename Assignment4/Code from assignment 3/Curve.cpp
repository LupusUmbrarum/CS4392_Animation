#include "Curve.h"

Curve::Curve()
{
}

Curve::~Curve()
{
}

void Curve::init(bool calc)
{
	this->control_points_pos = {
		{ 0.0, 8.5, -2.0 },
		{ -3.0, 11.0, 2.3 },
		{ -6.0, 8.5, -2.5 },
		{ -4.0, 5.5, 2.8 },
		{ 1.0, 2.0, -4.0 },
		{ 4.0, 2.0, 3.0 },
		{ 7.0, 8.0, -2.0 },
		{ 3.0, 10.0, 3.7 }
	};

	calculate_curve(calc);
}

void Curve::calculate_curve(bool calc)
{
	calcCurve = calc;

	// without this, I'd just keep adding to the list
	curve_points_pos.clear();

	if (calc)
	{
		// create the spline curve at each of the points
		for (int i = 0; i < 8; i++)
		{
			createSplineCurve(i);
		}
	}
	else
	{
		// create the straight lines between the points
		curve_points_pos = {
		{ 0.0, 8.5, -2.0 },
		{ -3.0, 11.0, 2.3 },
		{ -6.0, 8.5, -2.5 },
		{ -4.0, 5.5, 2.8 },
		{ 1.0, 2.0, -4.0 },
		{ 4.0, 2.0, 3.0 },
		{ 7.0, 8.0, -2.0 },
		{ 3.0, 10.0, 3.7 }
		};
	}

	// save the current number of points
	numPoints = curve_points_pos.size();

	// calculate the distances between points
	calculateDistances();

	// create the table containing all of the data
	createTable();

	setupOtherValues();
}

void Curve::createSplineCurve(int index)
{
	glm::vec3 p0, p1, p2, p3, prevPoint, newPoint;

	p0 = control_points_pos[clampForSpline(index - 1)];
	p1 = control_points_pos[index];
	p2 = control_points_pos[clampForSpline(index + 1)];
	p3 = control_points_pos[clampForSpline(index + 2)];

	prevPoint = p1;
	float t;

	for (int i = 1; i <= 200; i++)
	{
		// .005 being the resolution of the curve
		t = i * .005;

		// get the new point
		newPoint = calcPoint(t, p0, p1, p2, p3);

		// add it to the list of points
		curve_points_pos.push_back(newPoint);

		// set the previous point to the latest point
		prevPoint = newPoint;
	}
}

float Curve::distance(glm::vec3 a, glm::vec3 b)
{
	float x1, y1, z1, x2, y2, z2;
	x1 = a.x;
	y1 = a.y;
	z1 = a.z;
	x2 = b.x;
	y2 = b.y;
	z2 = b.z;

	float dx, dy, dz;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;

	return sqrt(dx * dx + dy * dy + dz * dz);
}

glm::vec3 Curve::calcPoint(float t, glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 u1, u2, u3, u4;

	u1 = 2.0f * p1;
	u2 = p2 - p0;
	u3 = 2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3;
	u4 = -p0 + 3.0f * p1 - 3.0f * p2 + p3;

	float t2 = t * t;
	float t3 = t2 * t;

	u2 *= t;
	u3 *= t2;
	u4 *= t3;

	return 0.5f * (u1 + u2 + u3 + u4);
}

float Curve::distance(int* a, int* b)
{
	glm::vec3 va, vb;
	va = glm::vec3(a[0], a[1], a[2]);
	vb = glm::vec3(b[0], b[1], b[2]);

	return distance(va, vb);
}

int Curve::clampForSpline(int x)
{
	if (x < 0)
	{
		x = control_points_pos.size() - 1;
	}

	if (x > control_points_pos.size())
	{
		x = 1;
	}
	else if (x > control_points_pos.size() - 1)
	{
		x = 0;
	}

	return x;
}

void Curve::calculateDistances()
{
	distanceBetweenPoints.clear();
	float prevDist = 0.0f;
	float curDist = 0.0f;

	for (int i = 0; i < numPoints - 1; i++)
	{
		curDist = distance(curve_points_pos[i], curve_points_pos[i + 1]);
		distanceBetweenPoints.push_back(curDist);
		prevDist = curDist;
		//std::cout << distanceBetweenPoints[i] << std::endl;
	}

	// get distance between first point and last point
	distanceBetweenPoints.push_back(distance(curve_points_pos[0], curve_points_pos[numPoints - 1]));
}

void Curve::tick(float deltaTime)
{

}

void Curve::createTable()
{
	table.clear();

	float uinc = 1.0f / (float)numPoints;

	table.push_back(TableEntry
	(
		0,
		0,
		0,
		0,
		0,
		curve_points_pos[0]
	));

	float prevTotalDistance = 0.0f;

	for (int i = 1; i < numPoints; i++)
	{
		int segment = i / (numPoints / control_points_pos.size());

		table.push_back(TableEntry
		(
			segment, // segment
			i, // index
			uinc * (float)i, // u
			distanceBetweenPoints[i] + prevTotalDistance, // length
			distanceBetweenPoints[i], // distance
			curve_points_pos[i]// glm::vec3 point
		));

		prevTotalDistance += distanceBetweenPoints[i];

		//std::cout << distanceBetweenPoints[i] << ", " << distanceBetweenPoints[i - 1] << std::endl;
		//table[i].print();
	}

	// get the last point. the one between the last index and the first index
	table.push_back(TableEntry
	(
		table.size() - 1 / (numPoints / control_points_pos.size()),
		table.size(), 
		uinc * (float)numPoints, 
		distanceBetweenPoints[distanceBetweenPoints.size() - 1] + prevTotalDistance,
		distanceBetweenPoints[distanceBetweenPoints.size() - 1],
		curve_points_pos[0]
	));

	tableSize = table.size();

	/*
	for (int i = 0; i < tableSize; i++)
	{
		table[i].print();
	}
	*/

	std::cout << "\n\n" << std::endl;

	//table[table.size() - 1].print();
}

float Curve::getNextDistance(float distance)
{
	for (int i = 0; i < tableSize; i++)
	{
		if (table[i].length > distance)
		{
			return table[i].length;
		}
	}

	return 0;
}

int Curve::getNextIndex(float distance)
{
	for (int i = 0; i < tableSize; i++)
	{
		if (table[i].length > distance)
		{
			return i;
		}
	}

	return 0;
}

float Curve::ease(float t)
{
	if (t < t1)
	{

	}
	
	if (t1 < t && t <= t2)
	{

	}

	// else, t must be greater than t2

	return (float)(glm::sin(t * PI - PI / 2) + 1) / 2;
}

// this function is more like a local reset
void Curve::setupOtherValues()
{
	adjt2 = ttotal - adjt2;
	aircraftPointIndex = 0;

	tcurrent = 0.0f;

	tinc = 1.0f / (float)numPoints;
}

float Curve::getTotalArcLength()
{
	return table[table.size() - 1].length;
}