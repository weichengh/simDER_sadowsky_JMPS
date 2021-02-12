#ifndef EXTERNALGRAVITYFORCE_H
#define EXTERNALGRAVITYFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalGravityForce
{
public:
	externalGravityForce(elasticRod &m_rod, timeStepper &m_stepper, Vector3d m_gVector);
	~externalGravityForce();
	void computeFg();
	void computeJg();

	void setGravity();
	Vector3d gVector;
	
private:
	elasticRod *rod;
	timeStepper *stepper;
    VectorXd massGravity;
};

#endif
