/**
 * @file 	dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------


Real dp = 0.0125;	/**< Initial reference particle spacing. */

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;					 /**< Reference density of fluid. */
Real rho0_v = 1.226;						 /**< Reference density of air. */
Real c_p_liquid = 4.179;
Real c_p_vapor = 1.012;
Real k_liquid = 0.620;
Real k_vapor = 0.0254;
Real diffusion_coff_liquid = k_liquid / (c_p_liquid * rho0_f);
Real diffusion_coff_vapor = k_vapor / (c_p_vapor * rho0_v);

Real c_f = 1.0;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric elements used in shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(40*dp, 0.0));
	water_block_shape.push_back(Vecd(40*dp, 40*dp));
	water_block_shape.push_back(Vecd(80*dp, 40*dp));
	water_block_shape.push_back(Vecd(80*dp, 0.0));
	water_block_shape.push_back(Vecd(40*dp, 0.0));
	return water_block_shape;
}

std::vector<Vecd> createAirBlockShape()
{
	std::vector<Vecd> air_block_shape;
	air_block_shape.push_back(Vecd(0.0, 0.0));
	air_block_shape.push_back(Vecd(0.0, 40*dp));
	air_block_shape.push_back(Vecd(40*dp, 40*dp));
	air_block_shape.push_back(Vecd(40*dp, 0.0));
	air_block_shape.push_back(Vecd(0.0, 0.0));

	return air_block_shape;
}
//----------------------------------------------------------------------
//	cases-dependent geometric shape for water block.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for air block.
//----------------------------------------------------------------------
class AirBlock : public MultiPolygonShape
{
public:
	explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createAirBlockShape(), ShapeBooleanOps::add);
	}
};


//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion fluid body
//----------------------------------------------------------------------
class ThermoWaterBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	ThermoWaterBodyMaterial()
		: DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_f,c_f)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff_liquid);
	};
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion solid body
//----------------------------------------------------------------------
class ThermoAirBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	ThermoAirBodyMaterial() : DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_v,c_f)
	{
		// only default property is given, as no heat transfer within solid considered here.
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi",diffusion_coff_vapor);
	};
};
//----------------------------------------------------------------------
//	Application dependent solid body initial condition
//----------------------------------------------------------------------
class ThermoAirBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoAirBodyInitialCondition(SPHBody &sph_body)
		: DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		species_n_[phi_][index_i] = 0.0;
		thermal_conductivity_[index_i] = k_vapor;
	};
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermoWaterBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoWaterBodyInitialCondition(SPHBody &sph_body)
		: DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
			species_n_[phi_][index_i] = 1.0;
			thermal_conductivity_[index_i] = k_liquid;
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
	explicit ThermalRelaxationComplex(ComplexRelation &body_complex_relation)
		: TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~ThermalRelaxationComplex() {};
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class TemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit TemperatureObserverParticleGenerator(SPHBody &sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 80;
		
		for(int i = 0;i< number_of_observation_points;i++)
		{
			positions_.push_back(Vecd(i*dp,20*dp));
		}
	}
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(80*dp, 40*dp));
	SPHSystem sph_system(system_domain_bounds, dp);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineAdaptation<SPHAdaptation>(1.0, 1.0);
	water_block.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoWaterBodyMaterial>();
	water_block.generateParticles<ParticleGeneratorLattice>();

	FluidBody vapor_block(sph_system, makeShared<AirBlock>("AirBody"));
	vapor_block.defineAdaptation<SPHAdaptation>(1.0, 1.0);
	vapor_block.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoAirBodyMaterial>();
	vapor_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
	temperature_observer.generateParticles<TemperatureObserverParticleGenerator>();

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_vapor_complex(water_block, { &vapor_block });
	ComplexRelation vapor_water_complex(vapor_block, { &water_block });
	ContactRelation temperature_observer_contact(temperature_observer, RealBodyVector{ &water_block,&vapor_block });

	//BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */

	SimpleDynamics<ThermoAirBodyInitialCondition> thermo_vapor_initial_condition(vapor_block);
	SimpleDynamics<ThermoWaterBodyInitialCondition> thermo_water_initial_condition(water_block);

	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(water_block);
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_vapor(vapor_block);

	//PeriodicConditionUsingCellLinkedList periodic_condition_water(water_block, water_block.getBodyShapeBounds(), yAxis);
	//PeriodicConditionUsingCellLinkedList periodic_condition_vapor(vapor_block, vapor_block.getBodyShapeBounds(), yAxis);
	ThermalRelaxationComplex thermal_relaxation_complex_water(water_vapor_complex);
	ThermalRelaxationComplex thermal_relaxation_complex_vapor(vapor_water_complex);

	/*ReducedQuantityRecording<ReduceDynamics<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		write_external_heat_water(io_environment, water_block, 2000);
	ReducedQuantityRecording<ReduceDynamics<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		write_external_heat_vapor(io_environment, vapor_block, 1006.43);*/
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	/*RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
	/*RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);*/
	ObservedQuantityRecording<Real> write_temperature("Phi", io_environment, temperature_observer_contact);

	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	//periodic_condition_water.update_cell_linked_list_.parallel_exec();
	//periodic_condition_vapor.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
		water_block.updateCellLinkedList();
		vapor_block.updateCellLinkedList();
		water_vapor_complex.updateConfiguration();
		vapor_water_complex.updateConfiguration();


		//fluid_observer_contact.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	int ite=0.0;
	Real end_time = 1.0;
	Real Output_Time = 0.1 * end_time;
	Real Observe_time =  0.05*Output_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_fluid_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------

	thermo_water_initial_condition.parallel_exec();
	thermo_vapor_initial_condition.parallel_exec();
	body_states_recording.writeToFile(0);
	write_temperature.writeToFile(0);

	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------

	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < Output_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Observe_time)
			{
				Real dt_thermal_water = get_thermal_time_step_water.parallel_exec();
				Real dt_thermal_air = get_thermal_time_step_vapor.parallel_exec();
				dt = SMIN(dt_thermal_water, dt_thermal_air);

				thermal_relaxation_complex_vapor.parallel_exec(dt);
				thermal_relaxation_complex_water.parallel_exec(dt);

				if (ite % 100== 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				ite++;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;
			/** Update cell linked list and configuration. */
			//periodic_condition_water.bounding_.parallel_exec();
			//periodic_condition_vapor.bounding_.parallel_exec();
			water_block.updateCellLinkedListWithParticleSort(100);
			vapor_block.updateCellLinkedListWithParticleSort(100);
			//periodic_condition_water.update_cell_linked_list_.parallel_exec();

			//periodic_condition_vapor.update_cell_linked_list_.parallel_exec();
			water_vapor_complex.updateConfiguration();
			vapor_water_complex.updateConfiguration();
			temperature_observer_contact.updateConfiguration();
			write_temperature.writeToFile();
			//write_temperature.writeToFile();
			time_instance = tick_count::now();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		
		
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}

	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;

	return 0;
};
