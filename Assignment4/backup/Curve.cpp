#include "Curve.h"

Curve::Curve()
{
}

Curve::~Curve()
{
}

void Curve::init(bool calc)
{
	calcCurve = calc;

	this->control_points_pos =
	{
		{ 0.0, 8.5, -2.0 },
		{ -3.0, 11.0, 2.3 },
		{ -6.0, 8.5, -2.5 },
		{ -4.0, 5.5, 2.8 },
		{ 1.0, 2.0, -4.0 },
		{ 4.0, 2.0, 3.0 },
		{ 7.0, 8.0, -2.0 },
		{ 3.0, 10.0, 3.7 }
	};

	this->control_points_quaternion = 
	{
		{0.13964   , 0.0481732 , 0.831429 , 0.541043 , },
		{0.0509038 , -0.033869 , -0.579695, 0.811295 , },
		{-0.502889 , -0.366766 , 0.493961 , 0.592445 , },
		{-0.636    , 0.667177  , -0.175206, 0.198922 , },
		{0.693492  , 0.688833  , -0.152595, -0.108237, },
		{0.752155  , -0.519591 , -0.316988, 0.168866 , },
		{0.542054  , 0.382705  , 0.378416 , 0.646269 , },
		{0.00417342, -0.0208652, -0.584026, 0.810619   }
	};

	calculate_curve();
}

void Curve::calculate_curve()
{
	curvePoints.clear();

	if (calcCurve)
	{
		// create the spline curve at each of the points
		for (int i = 0; i < 8; i++)
		{
			createSplineCurve(i);
		}
	}
	else
	{
		curvePoints = {
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

	calcAllQuats();

	// save the current number of points
	numPoints = curvePoints.size();

	// calculate the distances between points
	calculateDistances();

	// create the table containing all of the data
	createTable();

	setupOtherValues();
}

float Curve::angleBetweenPoints(glm::vec3 a, glm::vec3 b)
{
	return acos(
		(a.x * b.x + a.y * b.y + a.z * b.z) 
		/
		sqrt((a.x * a.x + a.y * a.y + a.z * a.z) * (b.x * b.x + b.y * b.y + b.z * b.z)));
}

glm::vec3 Curve::crossProduct(glm::vec3 a, glm::vec3 b)
{
	float cx, cy, cz;
	cx = a.y * b.z - a.z * b.y;
	cy = a.z * b.x - a.x * b.z;
	cz = a.x * b.y - a.y * b.x;
	return glm::vec3(cx, cy, cz);
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
		curvePoints.push_back(newPoint);

		// set the previous point to the latest point
		prevPoint = newPoint;
	}
}

void Curve::calculateDistances()
{
	distanceBetweenPoints.clear();
	float prevDist = 0.0f;
	float curDist = 0.0f;

	for (int i = 0; i < numPoints - 1; i++)
	{
		curDist = distance(curvePoints[i], curvePoints[i + 1]);
		distanceBetweenPoints.push_back(curDist);
		prevDist = curDist;
		//std::cout << distanceBetweenPoints[i] << std::endl;
	}

	// get distance between first point and last point
	distanceBetweenPoints.push_back(distance(curvePoints[0], curvePoints[numPoints - 1]));
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
		curvePoints[0]
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
			curvePoints[i]// glm::vec3 point
		));

		prevTotalDistance += distanceBetweenPoints[i];
	}

	// get the last point. the one between the last index and the first index
	table.push_back(TableEntry
	(
		table.size() - 1 / (numPoints / control_points_pos.size()),
		table.size(),
		uinc * (float)numPoints,
		distanceBetweenPoints[distanceBetweenPoints.size() - 1] + prevTotalDistance,
		distanceBetweenPoints[distanceBetweenPoints.size() - 1],
		curvePoints[0]
	));

	tableSize = table.size();
}

void Curve::setupOtherValues()
{
	adjt2 = ttotal - adjt2;
	aircraftPointIndex = 0;

	tcurrent = 0.0f;

	tinc = 1.0f / (float)numPoints;
}

int Curve::clampForSpline(int pos)
{
	if (pos < 0)
	{
		pos = control_points_pos.size() - 1;
	}

	if (pos > control_points_pos.size())
	{
		pos = 1;
	}
	else if (pos > control_points_pos.size() - 1)
	{
		pos = 0;
	}

	return pos;
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

float Curve::getTotalArcLength()
{
	return table[table.size() - 1].length;
}

float Curve::getVectorMagnitude(glm::vec3 a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
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

glm::mat4 Curve::quatToMat4(glm::quat q)
{
	float x, y, z, s;
	x = q.x;
	y = q.y;
	z = q.z;
	s = q.w;

	float x1, x2, x3, x4,
		y1, y2, y3, y4,
		z1, z2, z3, z4,
		w1, w2, w3, w4;

	x1 = 1 - (2 * (q.y * q.y)) - (2 * (q.z * q.z));
	x2 = 2 * x * y + 2 * s * z;
	x3 = 2 * x * z - 2 * s * y;
	x4 = 0;

	y1 = 2 * x * y - 2 * s * z;
	y2 = 1 - 2 * (x * x) - 2 * (z * z);
	y3 = 2 * y * z + 2 * s * x;
	y4 = 0;

	z1 = 2 * x * z + 2 * s * y;
	z2 = 2 * y * z - 2 * s * x;
	z3 = 1 - 2 * (x * x) - 2 * (y * y);
	z4 = 0;

	w1 = 0;
	w2 = 0;
	w3 = 0;
	w4 = 1;

	float arr[] = { x1, x2, x3, x4, y1, y2, y3, y4,
					z1, z2, z3, z4, w1, w2, w3, w4 };

	return glm::transpose(glm::make_mat4(arr));
}

void Curve::calcAllQuats()
{
	quats.clear();

	for (int i = 0; i < control_points_quaternion.size(); i++)
	{
		float x, y, z, w;
		glm::quat a, b;

		a = control_points_quaternion[i];

		if (i + 1 == control_points_quaternion.size())
		{
			b = control_points_quaternion[0];
		}
		else
		{
			b = control_points_quaternion[i + 1];
		}

		quats.push_back(a);

		for (int j = 1; j < 200; j++)
		{
			quats.push_back(slerp(a, b, (float)j / 200.f));
		}
	}
}

// I couldn't get slerp to work until I found this algorithm at
// https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
glm::quat Curve::slerp(glm::quat a, glm::quat b, float t)
{
	glm::quat q = glm::quat();

	float cosHalfTheta = 
		a.w * b.w +
		a.x * b.x +
		a.y * b.y +
		a.z * b.z;

	if (abs(cosHalfTheta) >= 1.f)
	{
		return a;
	}

	float halfTheta = acos(cosHalfTheta);
	float sinHalfTheta = sqrt(1.f - (cosHalfTheta * cosHalfTheta));

	if (fabs(sinHalfTheta) < .001f)
	{
		return a * .5f + b * .5f;
	}

	float ratioA = sin((1.f - t) * halfTheta) / sinHalfTheta;
	float ratioB = sin(t * halfTheta / sinHalfTheta);

	return a * ratioA + b * ratioB;
}