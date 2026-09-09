program package_interface

	use IFPORT
	
	use kind_parameters
	use global_data
	use computational_domain_class
	use chemical_properties_class
	use thermophysical_properties_class
	use solver_options_class
	use computational_mesh_class
	use mpi_communications_class
	use data_manager_class
	use boundary_conditions_class
	use field_scalar_class
	use field_vector_class
	use data_save_class
	use data_io_class
	use post_processor_manager_class
    use supplementary_routines

	implicit none
	
	type(computational_domain)					:: problem_domain
	type(data_manager)							:: problem_data_manager
	type(mpi_communications)					:: problem_mpi_support

	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics

	type(solver_options)						:: problem_solver_options
	
	type(computational_mesh)		,target		:: problem_mesh
	type(boundary_conditions)		,target		:: problem_boundaries
	type(field_scalar_cons)			,target		:: p, T, rho
	type(field_vector_cons)			,target		:: v, Y
	
	type(post_processor_manager)				:: problem_post_proc_manager
	
	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	character(len=5000)			:: initial_work_dir
	character(len=100)			:: work_dir
	character(len=20)			:: solver_name, coordinate_system, setup, setup_boundary
	
	real(dp)					:: u_inflow, inflow_radius
	real(dp)					:: ambient_temperature, ambient_pressure, ignition_temperature
    real(dp)                    :: delta_x
    real(dp)					:: domain_length, domain_height
    
	real(dp)	,dimension(3)	:: cell_size
	integer		,dimension(3,2)	:: utter_loop, observation_slice
	integer		,dimension(3,2)	:: upstream_transducer, downstream_transducer, summation_region
	integer	:: transducer_offset, front_lp, HPC, LPC, Ign_zone
	integer :: log_unit, ierr	
	
	integer	:: task1, task2, task3, task4, task5, task6, task7
	logical	:: stop_flag

	integer	:: dimensions, species_number
	integer	:: i,j,k, specie_number
	
	
	ierr = getcwd(initial_work_dir)
	
	do task1 = 1, 1
	do task2 = 1, 1

		work_dir = 'DGun_fill_acetylene_2D'
        
        u_inflow	= 2.0391_dp * task1
        work_dir	= trim(work_dir) // trim(fold_sep) //  trim(str_e(real(task1,dp))) //'_mps'

        inflow_radius = 5e-02_dp
        
		ierr = system('mkdir '// work_dir)
		
		select case(task2)
			case(1)
				work_dir = trim(work_dir) // trim(fold_sep) //  'dx_1.0e-03'
				delta_x	 = 5.0e-03_dp		
		end select
	
		ierr = system('echo d | xcopy .\' // trim(task_setup_folder) // ' .\' // trim(work_dir) // trim(fold_sep) // trim(task_setup_folder) // ' /E/K' )
		ierr = chdir(work_dir)
		
		open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
	
        domain_length	= 1.024_dp
        domain_height	= 1.024_dp * 0.25_dp
       
		problem_domain			= computational_domain_c(	dimensions 			=	2				,	&
															cells_number 		=	(/nint((domain_length)/delta_x),nint((domain_height)/delta_x),1/),			&
															coordinate_system	=	'cartesian',   &													
															lengths				=	reshape((/	0.0_dp,0.0_dp,0.0_dp,					&
																								domain_length, domain_height,0.005_dp/),(/3,2/)),	&
															axis_names			=	(/'x','y','z'/))
	
		problem_chemistry		= chemical_properties_c(chemical_mechanism_file_name	= 'ACETYLENE_Varatharajan.txt'	,	&	
														default_enhanced_efficiencies	= 1.0_dp			,	&
														E_act_units						= 'J.mol')
					
		problem_thermophysics	= thermophysical_properties_c(	chemistry					= problem_chemistry			,	&
																thermo_data_file_name		= 'FFCM-1_thermo.txt'		,	&	
																transport_data_file_name	= 'FFCM-1_transport.txt'	,	&	
																molar_masses_data_file_name	= 'molar_masses.dat')
	
		problem_mpi_support		= mpi_communications_c(problem_domain)
        
		problem_solver_options = solver_options_c(solver_name					= 'CABARET', &
													hydrodynamics_flag			= .true.		, &
													heat_transfer_flag			= .true.		, &
													molecular_diffusion_flag	= .true.		, &
													viscosity_flag				= .true.		, &
													chemical_reaction_flag		= .false.		, &
                                                    grav_acc                    = (/0.0_dp, 9.8_dp, 0.0_dp/)    , &
													CFL_flag					= .true.		, &
													CFL_coefficient				= 0.25_dp	, &
													initial_time_step			= 1e-08_dp	)	

		problem_data_manager	= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics,problem_solver_options)
	
			call problem_data_manager%create_boundary_conditions(	problem_boundaries				, &	
																	number_of_boundary_types	= 3	, &
																	default_boundary			= 1)
																
			call problem_data_manager%create_computational_mesh(problem_mesh)
			call problem_data_manager%create_scalar_field(p		,'pressure'						, 'p')
			call problem_data_manager%create_scalar_field(T		,'temperature'					, 'T')
			call problem_data_manager%create_scalar_field(rho	,'density'						, 'rho')
	
			call problem_data_manager%create_vector_field(v		,'velocity'						,'v'		,'spatial')
			call problem_data_manager%create_vector_field(Y		,'specie_molar_concentration'	,'Y'		,'chemical')	
	
			cell_size				= problem_mesh%get_cell_edges_length()
			utter_loop				= problem_domain%get_global_utter_cells_bounds()

			problem_post_proc_manager = post_processor_manager_c(problem_data_manager,number_post_processors = 0)
	
			problem_data_save	= data_save_c(	problem_data_manager	,	&
												visible_fields_names	= (/'pressure'						,&
																			'temperature'					,&
																			'density'						,&
																			'velocity'						,&
																			'specie_molar_concentration'	,&
																			'diffusivity'					,&
																			'thermal_conductivity'			,&
																			'viscosity'						,&
																			'velocity_production_viscosity'	,&
																			'energy_production_chemistry'	,&
																			'pressure_dynamic'/) ,	&
												save_time				= 1.0_dp		,	&
												save_time_units			= 'milliseconds'		,	&
												save_format				= 'tecplot'			,	&
												data_save_folder		= 'data_save'		,	&
												debug_flag				= .false.)
	
			problem_data_io		= data_io_c(problem_data_manager				    ,	&
											check_time			= 5000.0_dp	        ,	&
											check_time_units	= 'microseconds'	,	&
											data_output_folder	= 'data_output')

		!****************************** Setting initial conditions *************************************

		!# Ambient pressure and temperature		
            
            p%cells(:,:,:)	        = 1.0_dp * 101325.0_dp
            T%cells(:,:,:)	        = 300.0_dp
            
       		Y%pr(3)%cells(:,:,:)	= 1.0_dp			 
                    
		!***********************************************************************************************
			
			call problem_thermophysics%change_field_units_mole_to_dimless(Y)

		!****************************** Setting boundary conditions ************************************
			
			call problem_boundaries%create_boundary_type (	type_name				= 'outlet'	,	&
															farfield_pressure		= 101325.0_dp               ,	&
															farfield_temperature	= 300.0_dp	                ,	&
															farfield_velocity		= 0.0_dp   			        ,	&	
                											farfield_species_names	= (/'N2'/)					,	&	
															farfield_concentrations	= (/1.0_dp/)	,	&
															priority				= 1)					

			call problem_boundaries%create_boundary_type (	type_name				= 'inlet'			        ,	&
															farfield_pressure		= 101325.0_dp               ,	&
															farfield_temperature	= 300.0_dp	                ,	&
															farfield_velocity		= u_inflow			        ,	&	
															farfield_species_names	= (/'C2H2','O2'/)       ,	&	
															farfield_concentrations	= (/1.0_dp,1.0_dp/)	    ,	&	
															priority				= 2)	
            
			call problem_boundaries%create_boundary_type (	type_name				= 'wall'			        ,	&
															slip					= .false.	,	&		
															conductive				= .false.	,	&		
															wall_temperature		= 0.0_dp	,	&		
															wall_conductivity_ratio	= 0.0_dp	,	&		
															priority				= 3)           

        
            problem_boundaries%bc_markers(:,utter_loop(2,1),:)	= 3	
			problem_boundaries%bc_markers(:,utter_loop(2,2),:)	= 3
			problem_boundaries%bc_markers(utter_loop(1,1),:,:)	= 3	
			problem_boundaries%bc_markers(utter_loop(1,2),:,:)	= 1	

            do j = utter_loop(2,1), utter_loop(2,2)
                if (abs (j - 0.5*nint(domain_height/delta_x)) < nint(inflow_radius/delta_x)) then
  			        problem_boundaries%bc_markers(utter_loop(1,1),j,:)	= 2 
                end if  
            end do

		!***********************************************************************************************
	

		!****************************** Writing problem short description ******************************
	
		write(log_unit,'(A)') 'General description: Numerical modeling of ultra-lean hydrogen-air rising flame ball in terrestrial gravity conditions.'
		write(log_unit,'(A)') 'Main aim: To model flame dynamics in ultra-lean 4-6% H2-Air mixture using low-mach FDS solver.' 
		write(log_unit,'(A)') '			 Provide clear interpretation of the main stages of combustion process and analyse combustion wave structure and stability.'
		write(log_unit,'(A)') '			 Analyze geometrical factor effect on the flame evolution (axisymmetric vs cartesian).'
		write(log_unit,'(A)') '			 Obtain combustibility limits for different radii of ignition zone and different geometries.'
		write(log_unit,'(A)') '			 Perform comparison between numerical results and experimental ones, obtained by Golub et.al.'
		write(log_unit,'(A)') 'Problem setup: 4-6% H2-Air combustible mixture is ignited near horisontal solid wall. Garivation induced boyancy in vertical direction determines vertical rise.'
		write(log_unit,'(A)') 'Solver setup: FDS low mach solver, Keromnes oxidation scheme, gravity enabled, cartesian/cylindrical geometry.'
		write(log_unit,'(A)') 'Validation and comparison: Golub experiments on rising flames in large volumes.'
		write(log_unit,'(A)') '--------------------------------------------------------------------------'
	
		!***********************************************************************************************
														
		call problem_data_io%output_all_data(0.0_dp,stop_flag,make_output	= .true.)												
		call problem_data_save%save_all_data(0.0_dp,stop_flag,make_save		= .true.)
	
		ierr = chdir(initial_work_dir)	
		
		continue

	end do
	end do

end program