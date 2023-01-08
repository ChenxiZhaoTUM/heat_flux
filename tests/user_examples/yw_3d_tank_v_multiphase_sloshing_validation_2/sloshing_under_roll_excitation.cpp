/**
 * @file 	sloshing_under_roll_excitation.cpp
 * @brief 	Sloshing in Marine LNG Fuel Tank under Roll Excitation
 */
#include "sphinxsys.h"
#include "validation_under_roll_excitation.h"
 /*@brief Namespace cite here.
 */
using namespace SPH;



/*
Main program starts here.
*/
int main(int ac, char* av[])
{
	/* Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	system.setReloadParticles(true);
	/* Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);
	//handle command line arguments
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	/* Output environment. */
	IOEnvironment in_output(system);

	/*
	@Brief creating body, materials and particles for the cylinder.
	*/


	SolidBody tank(system, makeShared<Tank>("Tank"));
	tank.defineParticlesAndMaterial<SolidParticles, Solid>(); 
	//tank.defineBodyLevelSetShape()->writeLevelSet(in_output);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? tank.generateParticles<ParticleGeneratorReload>(in_output, tank.getName())
		: tank.generateParticles<ParticleGeneratorLattice>();
	tank.addBodyStateForRecording<Vecd>("NormalDirection");

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<Vecd>("Position");

	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f);
	air_block.generateParticles<ParticleGeneratorLattice>();
	air_block.addBodyStateForRecording<Real>("Pressure");

	ObserverBody tank_observer(system, "Tankobserver");
	tank_observer.generateParticles<TankObserverParticleGenerator>();

	InnerRelation tank_inner(tank);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_tank_to_vtp(in_output, {&tank});
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_particle_reload_files(in_output,  tank , "Tank" );
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner tank_relaxation_step_complex(tank_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_tank_particles.parallel_exec(0.25);
		tank_relaxation_step_complex.SurfaceBounding().parallel_exec();
		write_tank_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			tank_relaxation_step_complex.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
				write_tank_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_tank_particle_reload_files.writeToFile(0);
		return 0;
	}

	ContactRelation water_block_contact(water_block, { &tank });
	ContactRelation air_block_contact(air_block, { &tank });
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });
	ContactRelation tank_observer_contact(tank_observer, { &tank });
	ObservedQuantityRecording<Vecd> write_tank_move("Position", in_output, tank_observer_contact);
	ObservedQuantityRecording<Vecd> write_tank_nom("NormalDirection", in_output, tank_observer_contact);
	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);

	InteractionDynamics<solid_dynamics::CorrectConfiguration>		tank_corrected_configuration(tank_inner);
	SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g,0.0)));
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, makeShared<Gravity>(Vecd(0.0, -gravity_g,0.0)));
	/* Fluid dynamics */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(water_block_contact,water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_block_contact,air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_block_contact,air_water_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize>  fluid_advection_time_step(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize>  air_advection_time_step(air_block, U_max);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(air_block);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall> fluid_pressure_relaxation(water_block_contact,water_air_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall> fluid_density_relaxation(water_block_contact,water_air_complex);
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
		air_pressure_relaxation(air_block_contact, air_water_complex,2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		air_density_relaxation(air_block_contact, air_water_complex);

	BodyRegionByParticle wave_maker(tank, makeShared<Tank>("SloshingMaking"));
	SimpleDynamics<SloshMaking, BodyRegionByParticle> slosh_making(wave_maker);
	InteractionDynamics<InterpolatingAQuantity<Vecd>>
		interpolation_observer_position(tank_observer_contact, "Position", "Position");


	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	tank_corrected_configuration.parallel_exec();
	inner_normal_direction.parallel_exec();
	/** computing linear reproducing configuration for the tank. */

	write_real_body_states.writeToFile(0);
	write_tank_move.writeToFile(0);
	write_tank_nom.writeToFile(0);
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 17;			/**< End time. */
	Real D_Time = 0.025;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */
	Real dt_a = 0.0;				/**< Default acoustic time step sizes for air. */
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			initialize_a_water_step.parallel_exec();
			initialize_a_air_step.parallel_exec();

			Real Dt_f = fluid_advection_time_step.parallel_exec();
			Real Dt_a = air_advection_time_step.parallel_exec();
			Dt = SMIN(Dt_f, Dt_a);
			fluid_density_by_summation.parallel_exec();
			air_density_by_summation.parallel_exec();
			air_transport_correction.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = fluid_acoustic_time_step.parallel_exec();
				dt_a = air_acoustic_time_step.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_a),Dt);
				/* Fluid pressure relaxation */
				fluid_pressure_relaxation.parallel_exec(dt);
				air_pressure_relaxation.parallel_exec(dt);
				fluid_density_relaxation.parallel_exec(dt);
				air_density_relaxation.parallel_exec(dt);
				slosh_making.parallel_exec(dt);
				interpolation_observer_position.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}

			/** screen output, write body reduced values and restart files  */
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */

			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_contact.updateConfiguration();
			water_air_complex.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_block_contact.updateConfiguration();
			air_water_complex.updateConfiguration();
			tank.updateCellLinkedList();
			tank_observer_contact.updateConfiguration();

		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		write_real_body_states.writeToFile();
		write_tank_move.writeToFile();
		write_tank_nom.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
