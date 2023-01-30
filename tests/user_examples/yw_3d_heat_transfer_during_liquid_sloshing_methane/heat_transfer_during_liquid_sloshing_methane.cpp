/**
 * @file 	heat_transfer_during_liquid_sloshing.cpp
 * @brief 	Heat Transfer During Liquid Sloshing
 */
#include "sphinxsys.h"
#include "tank_case.h"
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
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);

	/*
	@Brief creating body, materials and particles for the tank, water, air, and sensors.
	*/
	SolidBody tank(system, makeShared<Tank>("Tank"));
	tank.defineParticlesAndMaterial<SolidParticles, Solid>();
	//tank.defineBodyLevelSetShape()->writeLevelSet(in_output);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? tank.generateParticles<ParticleGeneratorReload>(in_output, tank.getName())
		: tank.generateParticles<ParticleGeneratorLattice>();
	
	FluidBody methane_block(system, makeShared<LiquidBlock>("LiquidMethane"));
	methane_block.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoLiquidBodyMaterial>();
	methane_block.generateParticles<ParticleGeneratorLattice>();

	FluidBody nitrogen_block(system, makeShared<GasBlock>("Nitrogen"));
	nitrogen_block.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoGasBodyMaterial>();
	nitrogen_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody liquid_temperature_observer(system, "LiquidTemperatureObserver");
	liquid_temperature_observer.generateParticles<LiquidTemperatureObserverParticleGenerator>();

	ObserverBody gas_temperature_observer(system, "GasTemperatureObserver");
	gas_temperature_observer.generateParticles<GasTemperatureObserverParticleGenerator>();

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
		BodyStatesRecordingToVtp write_tank_to_vtp(in_output, { &tank });
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_particle_reload_files(in_output, tank, "Tank");
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_tank_particles.parallel_exec(0.25);
		tank_relaxation_step_inner.SurfaceBounding().parallel_exec();
		write_tank_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			tank_relaxation_step_inner.parallel_exec();
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

	ContactRelation methane_block_contact(methane_block,{ &tank });
	ContactRelation nitrogen_block_contact(nitrogen_block, { &tank });
	ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, { &methane_block });
	ContactRelation gas_temperature_observer_contact(gas_temperature_observer, { &nitrogen_block });
	ComplexRelation methane_nitrogen_complex(methane_block, { &nitrogen_block });
	ComplexRelation nitrogen_methane_complex(nitrogen_block, { &methane_block });

	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);
	ObservedQuantityRecording<Real> write_temperature_liquid("Phi", in_output, liquid_temperature_observer_contact);
	ObservedQuantityRecording<Real> write_temperature_gas("Phi", in_output, gas_temperature_observer_contact);
	
	SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
	SimpleDynamics<TimeStepInitialization> initialize_a_liquid_step(methane_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_gas_step(nitrogen_block, makeShared<VariableGravity>());
	/* Fluid dynamics */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(methane_block_contact, methane_nitrogen_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(nitrogen_block_contact, nitrogen_methane_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(nitrogen_block_contact, nitrogen_methane_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_nitrogen(nitrogen_block_contact, nitrogen_methane_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_methane(methane_block_contact, methane_nitrogen_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(methane_block, U_l);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> air_advection_time_step(nitrogen_block, U_n);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(methane_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(nitrogen_block);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> fluid_pressure_relaxation(methane_block_contact, methane_nitrogen_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> fluid_density_relaxation(methane_block_contact, methane_nitrogen_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
		gas_pressure_relaxation(nitrogen_block_contact, nitrogen_methane_complex, 2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		gas_density_relaxation(nitrogen_block_contact, nitrogen_methane_complex);

	SimpleDynamics<ThermoGasBodyInitialCondition> thermo_gas_initial_condition(nitrogen_block);
	SimpleDynamics<ThermoLiquidBodyInitialCondition> thermo_liquid_initial_condition(methane_block);

	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity_methane(methane_nitrogen_complex.getInnerRelation());
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity_nitrogen(nitrogen_methane_complex.getInnerRelation());

	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(methane_block);
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_air(nitrogen_block);

	ThermalRelaxationComplex thermal_relaxation_complex_water(methane_nitrogen_complex);
	ThermalRelaxationComplex thermal_relaxation_complex_air(nitrogen_methane_complex);
	
	BodyRegionByCell probe_s1(methane_block, makeShared<ProbeS1>("PorbeS1"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_1(in_output, probe_s1);
	BodyRegionByCell probe_s2(methane_block, makeShared<ProbeS2>("PorbeS2"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_2(in_output, probe_s2);
	BodyRegionByCell probe_s3(methane_block, makeShared<ProbeS3>("PorbeS3"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_3(in_output, probe_s3);
	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		water_average_temperature(in_output, methane_block,"Phi");
	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		air_average_temperature(in_output, nitrogen_block, "Phi");
	ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>>
		nitrogen_rate_of_heat_flow(in_output,nitrogen_block, "HeatFlux");
	ReducedQuantityRecording<ReduceAverage<QuantitySummation<Real>>>
		methane_rate_of_heat_flow(in_output, methane_block, "HeatFlux");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_nitrogen_total_mass(in_output,nitrogen_block, "MassiveMeasure");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_methane_total_mass(in_output, methane_block, "MassiveMeasure");
	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	inner_normal_direction.parallel_exec();
	/** computing linear reproducing configuration for the tank. */
	write_real_body_states.writeToFile(0);
	thermo_liquid_initial_condition.parallel_exec(0);
	thermo_gas_initial_condition.parallel_exec(0);
	probe_1.writeToFile(0);
	probe_2.writeToFile(0);
	probe_3.writeToFile(0);
	water_average_temperature.writeToFile(0);
	air_average_temperature.writeToFile(0);
	nitrogen_rate_of_heat_flow.writeToFile(0);
	methane_rate_of_heat_flow.writeToFile(0);
	compute_methane_total_mass.writeToFile(0);
	compute_nitrogen_total_mass.writeToFile(0);
	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		methane_block.updateCellLinkedList();
		nitrogen_block.updateCellLinkedList();
		methane_nitrogen_complex.updateConfiguration();
		methane_block_contact.updateConfiguration();
		nitrogen_methane_complex.updateConfiguration();
		nitrogen_block_contact.updateConfiguration();
	}

	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0;			/**< End time. */
	Real D_Time = 0.05;	/**< time stamps for output. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */

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
			initialize_a_liquid_step.parallel_exec();
			initialize_a_gas_step.parallel_exec();
			Real Dt_f = fluid_advection_time_step.parallel_exec();
			Real Dt_a = air_advection_time_step.parallel_exec();
			Real Dt = SMIN(Dt_f, Dt_a);
			fluid_density_by_summation.parallel_exec();
			air_density_by_summation.parallel_exec();
			air_transport_correction.parallel_exec();
			viscous_acceleration_nitrogen.parallel_exec();
			viscous_acceleration_methane.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = fluid_acoustic_time_step.parallel_exec();
				Real dt_a = air_acoustic_time_step.parallel_exec();
				Real dt_thermal_water = get_thermal_time_step_water.parallel_exec();
				Real dt_thermal_air = get_thermal_time_step_air.parallel_exec();
				dt = SMIN(SMIN(dt_f,dt_thermal_water),SMIN(dt_thermal_air, dt_a), Dt);
				/* Fluid pressure relaxation */
				fluid_pressure_relaxation.parallel_exec(dt);
				gas_pressure_relaxation.parallel_exec(dt);
				/* Fluid density relaxation */
				fluid_density_relaxation.parallel_exec(dt);
				gas_density_relaxation.parallel_exec(dt);
				thermal_relaxation_complex_air.parallel_exec(dt);
				thermal_relaxation_complex_water.parallel_exec(dt);
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

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			methane_block.updateCellLinkedListWithParticleSort(100);
			methane_block_contact.updateConfiguration();
			methane_nitrogen_complex.updateConfiguration();

			nitrogen_block.updateCellLinkedListWithParticleSort(100);
			nitrogen_block_contact.updateConfiguration();
			nitrogen_methane_complex.updateConfiguration();

			liquid_temperature_observer_contact.updateConfiguration();
			gas_temperature_observer_contact.updateConfiguration();

			int num = methane_nitrogen_complex.getInnerRelation().base_particles_.pos_.size();
			int count_num = 0;
			for (size_t k = 0; k < methane_nitrogen_complex.contact_configuration_.size(); ++k)
			{
				for (int i = 0; i != num; ++i)
				{
					Neighborhood& neighborhood = methane_nitrogen_complex.contact_configuration_[k][i];
					if (neighborhood.current_size_ != 0)
					{
						count_num += 1;
					}
				}
			}
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::string output_folder_ = "./output";
				std::string filefullpath = output_folder_ + "/" + "number_of_contact_particles" + ".dat";
				std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
				out_file_ << "\n";
				out_file_ << "Time " << "Num" << std::endl;
				out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num << "  " << std::endl;
			}
		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity_methane.parallel_exec();
		compute_vorticity_nitrogen.parallel_exec();
		write_real_body_states.writeToFile();
		probe_1.writeToFile();
		probe_2.writeToFile();
		probe_3.writeToFile();
		write_temperature_liquid.writeToFile();
		write_temperature_gas.writeToFile();
		water_average_temperature.writeToFile();
		air_average_temperature.writeToFile();
		nitrogen_rate_of_heat_flow.writeToFile();
		methane_rate_of_heat_flow.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
/*std::string output_folder_ = "./output";
+		std::string filefullpath = output_folder_ + "/" + "check_particle_infomation_" + ".dat";
+		std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
+		out_file_ << "\n";
+		out_file_ << "index_i, " << "position_x, " << "position_y, " << "volume, " << "smoothing_length_ratio, " << std::endl;
+		for (size_t i = 0; i != base_particles_->total_real_particles_; i++) {
+			out_file_ << i << "  " << base_particles_->pos_n_[i][0] << "  " << base_particles_->pos_n_[i][1] << "  "
+				<< base_particles_->Vol_[i] << " " << h_ratio_[i] << std::endl;
+		}*/