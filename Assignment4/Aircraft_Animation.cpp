#include "Aircraft_Animation.h"


Aircraft_Animation::Aircraft_Animation()
{
	this->m_model_mat = glm::mat4(1.0f);
}


Aircraft_Animation::~Aircraft_Animation()
{
}

void Aircraft_Animation::init()
{
	reset();
}

void Aircraft_Animation::init(Curve* animation_curve)
{
	m_curve = animation_curve;
	reset();
}

void Aircraft_Animation::update(float delta_time)
{
	if (!m_curve->aircraftIsMoving)
	{
		return;
	}

	checkAcceleration();

	// this sets the real time to logical time
	tcurrent += (delta_time)*tinc;

	distanceCovered = m_curve->getNextDistance(distanceCovered + getDistance());

	calcV0();
	
	// overwriting previous command
	pointIndex = tcurrent;

	if (pointIndex >= m_curve->numPoints)
	{
		pointIndex = 0;
		tcurrent = 0;
	}

	move(pointIndex);
}

void Aircraft_Animation::move(int pointIndex)
{
	glm::quat m1, m2;

	m1 = m_curve->quats[pointIndex / 200];

	// if the next point is out of bounds, get first point
	if (pointIndex + 1 > m_curve->numPoints)
	{
		m2 = m_curve->quats[0];
	}
	else
	{
		m2 = m_curve->quats[(pointIndex / 200) + 1];
	}

	if (m_curve != nullptr && m_curve->curvePoints.size() > 0)
	{
		m_model_mat = glm::mat4(1.0f);

		glm::quat q = m_curve->quats[pointIndex];

		glm::mat4 rotmat = glm::toMat4(q);

		m_model_mat = glm::translate(m_model_mat, m_curve->curvePoints[pointIndex]);
		
		m_model_mat *= rotmat;
	}
}

void Aircraft_Animation::reset()
{
	move(0);
	setAnimationVariables();
}

float Aircraft_Animation::ease(float t)
{
	return (sin((t * PI) - (PI / 2)) + 1) / 2;
	/*
	if (t == 0.0f)
	{
		return v0;
	}

	if (t > 0.0f && t < t1)
	{
		return v = v0 - (t / t1);
	}

	if (t > t1 && t < adjt2)
	{
		return v = v0;
	}

	if (t > adjt2 && t < 1.0f)
	{
		return v = v0 * ((t - t2) / (1.0f - t2));
	}

	return .01f;
	*/
}

void Aircraft_Animation::setAnimationVariables()
{
	Curve* c = m_curve;

	t1 = c->t1;
	t2 = c->t2;
	ttotal = c->ttotal;

	tmid = ttotal - t1 - t2;

	timeSegmentSize = c->table.size() / ttotal;

	tcurrent = 0.0f;

	tinc = timeSegmentSize;

	v0 = 1.0f / timeSegmentSize;

	adjt2 = ttotal - t2;

	pointIndex = 0;
	distanceCovered = 0.0f;

	initV0();

	maxVelocity = getMaxVelocity();
}

void Aircraft_Animation::checkAcceleration()
{
	if (tcurrent > 0.0f && tcurrent < t1)
	{
		acc = 1.0f;
		return;
	}

	if (tcurrent > t1&& tcurrent < adjt2)
	{
		acc = 0.0f;
		return;
	}

	if (tcurrent > adjt2&& tcurrent < 1.0f)
	{
		-1.0f;
		return;
	}

	acc = .01f;
}

float Aircraft_Animation::getDistance()
{
	if (tcurrent > 0.0f && tcurrent <= t1)
	{
		return v0 * ((tcurrent * tcurrent) / (2 * t1));
	}

	if (tcurrent > t1&& tcurrent <= t2)
	{
		return v0 * t1 / 2 + (v0 * (tcurrent - t2));
	}

	if (tcurrent > t2&& tcurrent <= 1.0)
	{
		return v0 * t1 / 2 + v0 * (t2 - t1) + v0 * (1 - ((tcurrent - t2) / (2 * (1 - t2))) * (tcurrent - t2));
	}

	return .0f;
}

float Aircraft_Animation::getMaxVelocity()
{
	return (2.0f / (adjt2 - t1 + 1.0));
}

void Aircraft_Animation::initV0()
{
	v0 = 1.0f / ((t1 / 2.0f) + (adjt2 - t1) + ((1.0f - ((1.0f - adjt2) / (2.0f * (1.0f - adjt2)))) * (1.0f - adjt2)));
}

void Aircraft_Animation::calcV0()
{
	v0 = distanceCovered / ((t1 / 2.0f) + (adjt2 - t1) + ((1.0f - ((tcurrent - adjt2) / (2.0f * (1.0f - adjt2)))) * (tcurrent - adjt2)));
}