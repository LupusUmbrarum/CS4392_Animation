#pragma once

#include <vector>
#include <iostream>

 #define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Curve.h"

class Aircraft_Animation
{
public:
	// ttotal is the total amount of time it should take to traverse the loop,
	// tcurrent is the current time that we're at, and tinc is the value that 
	// tcurrent will be incremented by
	float t1 = 0.5f, t2 = 0.5f, ttotal = 5.0f, tmid = 4.0f, tcurrent = 0.0f, tinc = 0.0f, v0 = 0.0f, v = 0.0f;
	// adjt2 is the adjusted time for t2, it equals ttotal - t2, giving the actual time position of t2
	float adjt2;
	float distanceCovered = 0.0f;
	float maxVelocity = 0.0f, currentVelocity = 0.0f;

	float acc = 0.0f;

	float timeSegmentSize = 0;

	int pointIndex = 0;
	float total_moving_time = 10;

	Aircraft_Animation();
	~Aircraft_Animation();

	void init();
	void init(Curve* animation_curve);

	void update(float delta_time);

	void reset();
	void move(int pointIndex);
	float ease(float t);
	void setAnimationVariables();
	void checkAcceleration();
	float getDistance();

	float getMaxVelocity();
	void initV0();
	void calcV0();
	glm::mat4 get_model_mat() { return m_model_mat; };

private:
	glm::mat4 m_model_mat;
	Curve* m_curve = nullptr;

};


