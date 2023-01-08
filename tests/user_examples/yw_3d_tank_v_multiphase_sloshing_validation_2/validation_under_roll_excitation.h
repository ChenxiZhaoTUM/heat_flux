/*
* @file 	validation_under_roll_excitation.h
*/

#ifndef	VALIDATION_UNDER_ROLL_EXCITATION_H
#define VALIDATION_UNDER_ROLL_EXCITATION_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/*@brief Basic geometry parameters and numerical setup.
*/

Real resolution_ref = 0.006;   /* Initial particle spacing*/


/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vecd(-0.6, -0.2,-0.2), Vecd(0.6, 0.4,0.2));


/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;         /*Fluid density*/
Real rho0_a = 1.0;            /*Air density*/
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_max = 2.0 * sqrt(gravity_g*0.0612);	    /* Characteristic velocity*/
Real c_f = 10.0 *U_max;					/*Reference sound speed*/
Real mu_f = 1.01e-3;
Real mu_a = 17.9e-6;

Real length_scale = 1.0;
Vecd translation(0, 0.12, 0);
/*
@Brief material properties of the solid,
*/

std::string fuel_tank_outer = "./input/3D_grotle_tank_outer_03.STL";
std::string fuel_tank_inner = "./input/3D_grotle_tank_inner.STL";
std::string water_05 = "./input/3D_grotle_water_0255.STL";
std::string air_05 = "./input/3D_grotle_air_0255.STL";

/*
The Tank.
*/
class Tank : public ComplexShape
{
public:
	explicit Tank(const std::string &shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */

		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale, "OuterWall");
		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale, "InnerWall");

	}
};
/*
The Water
*/
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};
/*
The Air
*/
class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
	}
};
/*
External Excitation
*/
class SloshMaking : public solid_dynamics::BaseMotionConstraint
{

	Vecd getDisplacement(const Real &time, const Vecd &pos_, const Vecd &pos0_)
	{
		Vecd displacement(0);
		Vecd rotation(0);
		Real angular = 0;
		angular = getAngular(time);

		rotation[0] = pos0_[0] * cos(angular) - pos0_[1] * sin(angular);
		rotation[1] = pos0_[0] * sin(angular) + pos0_[1] * cos(angular);
		rotation[2] = pos0_[2];
		return rotation;
	}

	Vecd getVelocity(const Real &time, const Vecd &pos_, const Vecd &pos0_)
	{
		Vecd velocity(0);
		Real x = pos0_[0];
		Real y = pos0_[1];
		Real local_radius = sqrt(pow(y, 2.0) + pow(x, 2.0));
		Real angular = atan2(y, x);

		Real angular_velocity = -2.0*PI*0.916*0.6*PI*cos(time * 2.0*PI*0.916*0.6) / 60.0;

		velocity[0] = angular_velocity * local_radius*cos(angular + PI / 2);
		velocity[1] = angular_velocity * local_radius*sin(angular + PI / 2);
		velocity[2] = 0.0;
		return velocity;
	}

	Real getAngular(const Real &time )
	{
		Real angular = 0;
		angular = -PI * sin(time * 2.0*PI*0.916*0.6) / 60.0;
		return angular;
	}

public:
	SloshMaking(BodyPartByParticle &constrained_region)
		:solid_dynamics::BaseMotionConstraint(constrained_region)
	{}
	void update(size_t index_i, Real dt=0.0)
	{
		Real time = GlobalStaticVariables::physical_time_;
		Real angular = getAngular(time);
		n_[index_i][0] = n0_[index_i][0] * cos(angular) - n0_[index_i][1] * sin(angular);
		n_[index_i][1] = n0_[index_i][0] * sin(angular) + n0_[index_i][1] * cos(angular);
		n_[index_i][2] = n0_[index_i][2];
		acc_[index_i] = Vecd(0.0,0.0,0.0);
		vel_[index_i] = getVelocity(time, pos_[index_i], pos0_[index_i]);
		pos_[index_i] = getDisplacement(time, pos_[index_i], pos0_[index_i]);
	};
};

class TankObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit TankObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		positions_.push_back(Vecd(-0.198, 0.0,0.0));
	}
};

#endif //VALIDATION_UNDER_ROLL_EXCITATION_H