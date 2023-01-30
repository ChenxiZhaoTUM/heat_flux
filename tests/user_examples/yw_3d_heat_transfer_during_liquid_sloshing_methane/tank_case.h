/*
* @file 	waterentry_elastic_case.h
*/

#ifndef	TANK_CASE_H
#define TANK_CASE_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.2, -0.05,-0.2), Vec3d(0.2, 1.0,0.2));

Real resolution_ref = 0.0075;   /* Initial particle spacing*/
/*
Material properties of the fluid.
*/
Real rho0_l = 422.6;         /*Liquid methane density*/
Real rho0_n = 2.601;            /*Nitrogen density*/
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_l = 2.0* sqrt(gravity_g * 0.5);	/**< Characteristic velocity. */
Real U_n = 2.0* sqrt(gravity_g * 0.5);  
Real c_f = 10.0 * SMAX(U_l, U_n);	/**< Reference sound speed. */
Real f = 1.0;
Real a = 0.1;
Real c_p_methane = 3.1397e3;
Real c_p_nitrogen = 1.040e3;
Real k_methane = 0.1846;
Real k_nitrogen = 0.01231;
Real diffusion_coff_methane = k_methane / (c_p_methane * rho0_l);
Real diffusion_coff_nitrogen = k_nitrogen / (c_p_nitrogen * rho0_n);
Real mu_methane = 118.3e-6;
Real mu_nitrogen = 8.39e-6;
Real length_scale = 1.0;
Vec3d translation(0, 0.0, 0);

std::string fuel_tank_outer = "./input/tank_outer.STL";
std::string fuel_tank_inner = "./input/tank_inner.STL";
std::string water_05 = "./input/water_05.STL";
std::string air_05 = "./input/gas_05.STL";
std::string probe_s1_shape = "./input/ProbeS1.STL";
std::string probe_s2_shape = "./input/ProbeS2.STL";
std::string probe_s3_shape = "./input/ProbeS3.STL";

/*
The Tank.
*/
class Tank : public ComplexShape
{
public:
	explicit Tank(const std::string &shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */

		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale,"OuterWall");
		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale,"InnerWall");
	}
};

/*
The Water.
*/
class LiquidBlock : public ComplexShape
{
public:
	explicit LiquidBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};

/*
The Air.
*/
class GasBlock : public ComplexShape
{
public:
	explicit GasBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
	}
};
/*
External Excitation.
*/
class VariableGravity : public Gravity
{
	Real time_ = 0;
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g, 0.0)) {};
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		time_= GlobalStaticVariables::physical_time_;

		global_acceleration_[0] = 4.0*PI*PI*f*f*a*sin(2 * PI*f*time_);
		
		return global_acceleration_;
	}
};

/*
Sensors: S1, S2, and S3
*/
class ProbeS1 : public ComplexShape
{
public:
	explicit ProbeS1(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
	}
};

class ProbeS2 : public ComplexShape
{
public:
	explicit ProbeS2(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_2(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
	}
};

class ProbeS3 : public ComplexShape
{
public:
	explicit ProbeS3(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_3(0.0,0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
	}
};

//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion water body
//----------------------------------------------------------------------
class ThermoLiquidBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	ThermoLiquidBodyMaterial()
		: DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_l, c_f, mu_methane)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff_methane);
	};
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion air body
//----------------------------------------------------------------------
class ThermoGasBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	ThermoGasBodyMaterial() : DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_n, c_f,mu_nitrogen)
	{
		// only default property is given, as no heat transfer within solid considered here.
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff_nitrogen);
	};
};

//----------------------------------------------------------------------
//	Application dependent water body initial condition
//----------------------------------------------------------------------
class ThermoLiquidBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoLiquidBodyInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		species_n_[phi_][index_i] = 111.15;
		thermal_conductivity_[index_i] = k_methane;
	};
};
//----------------------------------------------------------------------
//	Application dependent air body initial condition
//----------------------------------------------------------------------
class ThermoGasBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoGasBodyInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		species_n_[phi_][index_i] = 131.15;
		thermal_conductivity_[index_i] = k_nitrogen;
	};
};

//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
class ThermalRelaxationComplex
	: public TwoPhaseRelaxationOfAllDiffusionSpeciesRK2<
	TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<
	FluidParticles, WeaklyCompressibleFluid, FluidParticles, WeaklyCompressibleFluid>>
{
public:
	explicit ThermalRelaxationComplex(ComplexRelation& body_complex_relation)
		: TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~ThermalRelaxationComplex() {};
};

class LiquidTemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit LiquidTemperatureObserverParticleGenerator(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
			positions_.push_back(Vecd(0.0, 0.35, 0.0));
	}
};


class GasTemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit GasTemperatureObserverParticleGenerator(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		positions_.push_back(Vecd(0.0, 0.65, 0.0));
	}
};

#endif //TANK_CASE_H