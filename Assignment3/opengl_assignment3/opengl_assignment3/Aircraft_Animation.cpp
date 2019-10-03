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
	m_animation_curve = animation_curve;
	reset();
}

void Aircraft_Animation::update(float delta_time)
{
	if (!m_animation_curve->aircraftIsMoving)
	{
		return;
	}

	checkAcceleration();

	//std::cout << ease(tcurrent * delta_time) << std::endl;

	//tcurrent += ease(tcurrent * delta_time + delta_time);

	//tcurrent += v0;// *delta_time;
	pointIndex = tcurrent;
	//pointIndex++;

	//std::cout << tcurrent << " " << timeSegmentSize << std::endl;

	if (pointIndex >= m_animation_curve->curve_points_pos.size())
	{
		pointIndex = 0;
		tcurrent = 0;
	}
	std::cout << getDistance() << std::endl;
	pointIndex = m_animation_curve->getNextPoint(getDistance());
	std::cout << pointIndex << std::endl;
	move(pointIndex);
}

void Aircraft_Animation::move(int pointIndex)
{
	m_model_mat = glm::mat4(1.0f);

	if (m_animation_curve != nullptr && m_animation_curve->curve_points_pos.size() > 0)
	{
		m_model_mat = glm::translate(m_model_mat, m_animation_curve->curve_points_pos[pointIndex]);
	}
}

void Aircraft_Animation::reset()
{
	move(0);
	setAnimationVariables();
}

float Aircraft_Animation::ease(float t)
{
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
}

void Aircraft_Animation::setAnimationVariables()
{
	Curve* c = m_animation_curve;

	t1 = c->t1;
	t2 = c->t2;
	ttotal = c->ttotal;

	tmid = ttotal - t1 - t2;

	timeSegmentSize = c->table.size() / ttotal;

	tcurrent = 0.0f;
	tinc = 1.0f;
	v0 = 1.0f / timeSegmentSize;

	std::cout << "v0: " << v0 << std::endl;

	adjt2 = ttotal - t2;
}

void Aircraft_Animation::checkAcceleration()
{
	if (tcurrent > 0.0f && tcurrent < t1)
	{
		acc = 1.0f;
		return;
	}

	if (tcurrent > t1 && tcurrent < adjt2)
	{
		acc = 0.0f;
		return;
	}

	if (tcurrent > adjt2 && tcurrent < 1.0f)
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

	if (tcurrent > t1 && tcurrent <= t2)
	{
		return v0 * t1 / 2 + (v0 * (tcurrent - t2));
	}

	if (tcurrent > t2 && tcurrent <= 1.0)
	{
		return v0 * t1 / 2 + v0 * (t2 - t1) + v0 * (1 - ((tcurrent - t2) / (2 * (1 - t2))) * (tcurrent - t2));
	}

	return .01f;
}