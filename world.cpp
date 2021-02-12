#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	speed = m_inputData.GetScalarOpt("speed");
	width = m_inputData.GetScalarOpt("width");
	epsilon = m_inputData.GetScalarOpt("epsilon");
	inputShear = m_inputData.GetScalarOpt("inputShear");

	totalShear = 0.0;
	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) + 0.5);
    eta_par = 2.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) - 0.5);

    distance = 0.0;

    inputShear = inputShear * RodLength;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER";
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
		return;

	if ( timeStep % 10 != 0)
	{
	//	return;
	}

	if (currentTime > 2.0)
	{
		Vector3d xMid = rod->getVertex(rod->nv / 2);
		computeReactionForce();

		outfile << totalShear/RodLength << " " << xMid(2)/RodLength << endl;
	}
}

void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, theta, width);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = pow(rodRadius, 3) * width / 12.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);

//	m_bendingForce = new elasticBendingForce(*rod, *stepper);
//	m_twistingForce = new elasticTwistingForce(*rod, *stepper);

	m_elasticRibbonForce = new elasticRibbonForce(*rod, *stepper, epsilon);
	m_inertialForce = new inertialForce(*rod, *stepper);	
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	reactionForce = VectorXd::Zero(rod->ndof);
}

// Setup geometry
void world::rodGeometry()
{
    double deltaLength = RodLength / (numVertices - 1);

    vertices = MatrixXd(numVertices, 3);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = rand()%10 * 1e-8;
    	vertices(i, 1) = deltaLength * i;
    	vertices(i, 2) = rand()%10 * 1e-8;
    } 

    // initial theta should be zeros
    theta = VectorXd::Zero(numVertices - 1);
}

void world::rodBoundaryCondition()
{
	// Apply boundary condition
	rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	rod->setThetaBoundaryCondition(rod->getTheta(0), 0);
	rod->setVertexBoundaryCondition(rod->getVertex(1), 1);

	rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 2), rod->nv - 2);
	rod->setThetaBoundaryCondition(rod->getTheta(rod->nv - 2), rod->nv - 2);
	rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 1), rod->nv - 1);
}
	

void world::updateTimeStep()
{

	// compress
	if (distance <= RodLength/2)
	{
		Vector3d x_0 = rod->getVertex(0);
		Vector3d x_1 = rod->getVertex(1);

		x_0(1) = x_0(1) + deltaTime * 0.03;
		x_1(1) = x_1(1) + deltaTime * 0.03;

		rod->setVertexBoundaryCondition(x_0, 0);
		rod->setVertexBoundaryCondition(x_1, 1);

		Vector3d x_2 = rod->getVertex(rod->nv - 2);
		Vector3d x_3 = rod->getVertex(rod->nv - 1);

		x_2(1) = x_2(1) - deltaTime * 0.03;
		x_3(1) = x_3(1) - deltaTime * 0.03;	

		rod->setVertexBoundaryCondition(x_2, rod->nv - 2);
		rod->setVertexBoundaryCondition(x_3, rod->nv - 1);	

		distance = distance + 2 * deltaTime * 0.03;
	}

	// shear
	if (currentTime > 2.0 && totalShear < inputShear)
	{
		Vector3d x_0 = rod->getVertex(0);
		Vector3d x_1 = rod->getVertex(1);

		x_0(0) = x_0(0) + deltaTime * speed;
		x_1(0) = x_1(0) + deltaTime * speed;

		rod->setVertexBoundaryCondition(x_0, 0);
		rod->setVertexBoundaryCondition(x_1, 1);

		Vector3d x_2 = rod->getVertex(rod->nv - 2);
		Vector3d x_3 = rod->getVertex(rod->nv - 1);

		x_2(0) = x_2(0) - deltaTime * speed;
		x_3(0) = x_3(0) - deltaTime * speed;	

		rod->setVertexBoundaryCondition(x_2, rod->nv - 2);
		rod->setVertexBoundaryCondition(x_3, rod->nv - 1);	

		totalShear = totalShear + 2 * deltaTime * speed;
	}

	// compute
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();

		m_elasticRibbonForce->computeFribbon();
		m_elasticRibbonForce->computeJribbon();

		if (currentTime > 1.0)
		{
			m_gravityForce->setGravity();

			m_gravityForce->gVector(0) = 10.0;
			m_gravityForce->gVector(1) = 0.0;
			m_gravityForce->gVector(2) = 0.0;
		}

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();
		
		m_dampingForce->computeFd();
		m_dampingForce->computeJd();

		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);
		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rod->updateTimeStep();

	if (render) cout << "time: " << currentTime << " iter=" << iter << endl;

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}


}

int world::simulationRunning()
{
	if (timeStep<Nstep) 
		return 1;
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] / RodLength;
}

double world::getBoundaryCoordination_left(int i)
{
	int j = i / 4;
	int k = i - 4 * j;

	return	(rod->x[i] + 0.5 * rod->m1(j, k) * width) / RodLength;
}

double world::getBoundaryCoordination_right(int i)
{
	int j = i / 4;
	int k = i - 4 * j;

	return	(rod->x[i] - 0.5 * rod->m1(j,k) * width) / RodLength;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

void world::computeReactionForce()
{
	reactionForce = VectorXd::Zero(rod->ndof);

	reactionForce = m_inertialForce->forceVec - m_dampingForce->forceVec - m_elasticRibbonForce->forceVec - m_stretchForce->forceVec;
}