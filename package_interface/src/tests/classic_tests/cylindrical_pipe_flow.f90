!================================================================================
! CYLINDRICAL PIPE FLOW SIMULATION INTERFACE
!================================================================================
!
! PROGRAM: package_interface
!
! DESCRIPTION:
!   This program sets up a 2D axisymmetric simulation of gas flow 
!   in a cylindrical pipe. The computational domain uses cylindrical
!   coordinates (r, z) where:
!   - r (radial): from axis (r=0) to pipe wall (r=R)
!   - z (axial): from inlet (z=0) to outlet (z=L)
!
! KEY FEATURES:
!   - 2D axisymmetric formulation in cylindrical coordinates
!   - Symmetry boundary condition at pipe axis (r=0)
!   - No-slip wall boundary at pipe wall (r=R)
!   - Inlet velocity profile at entrance (z=0)
!   - Pressure outlet at exit (z=L)
!   - Support for laminar and turbulent flow regimes
!
! PHYSICAL SYSTEM:
!   Simulates gas flow (air or mixture) through a cylindrical pipe
!   Typical applications: pipe flow validation, pressure drop studies,
!   boundary layer development
!
!================================================================================

program package_interface

    !==========================================
    ! MODULE IMPORTS
    !==========================================
    use iso_c_binding, only: c_int, c_size_t, c_char, c_ptr, c_null_char, c_associated
    use kind_parameters            ! Defines precision kinds (dp, sp, etc.)
    use global_data                ! Global constants and parameters
    use nrg_build_info, only: write_nrg_source_revision
    use computational_domain_class ! Domain definition and management
    use chemical_properties_class  ! Chemical kinetics and species data
    use thermophysical_properties_class ! Thermodynamic and transport properties
    use solver_options_class       ! Numerical solver configuration
    use problem_control_class      ! Problem-specific active controls
    use computational_mesh_class   ! Mesh generation and management
    use mpi_communications_class   ! Parallel communication routines
    use data_manager_class         ! Central data management
    use boundary_conditions_class  ! Boundary condition specification
    use field_scalar_class         ! Scalar field operations (p, T, rho, etc.)
    use field_vector_class         ! Vector field operations (v, Y, etc.)
    use data_save_class            ! Data output and saving
    use data_io_class              ! Input/output operations
    use post_processor_manager_class ! Post-processing and monitoring
    use supplementary_routines     ! Utility functions and helper routines

    implicit none

    ! C runtime bindings used only for changing/querying the process working
    ! directory.  File-tree creation/copy is delegated to portable `cmake -E`
    ! commands below.
    interface
#ifdef WIN
        function c_getcwd(buffer, maxlen) bind(C, name="_getcwd") result(ptr)
            import :: c_ptr, c_char, c_int
            character(kind=c_char) :: buffer(*)
            integer(c_int), value :: maxlen
            type(c_ptr) :: ptr
        end function c_getcwd

        function c_chdir(path) bind(C, name="_chdir") result(status)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
#else
        function c_getcwd(buffer, maxlen) bind(C, name="getcwd") result(ptr)
            import :: c_ptr, c_char, c_size_t
            character(kind=c_char) :: buffer(*)
            integer(c_size_t), value :: maxlen
            type(c_ptr) :: ptr
        end function c_getcwd

        function c_chdir(path) bind(C, name="chdir") result(status)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
#endif
    end interface
    
    !==========================================
    ! PRIMARY SIMULATION OBJECTS
    !==========================================
    type(computational_domain)              :: problem_domain        ! Computational domain
    type(data_manager)                      :: problem_data_manager  ! Central data coordinator
    type(mpi_communications)                :: problem_mpi_support   ! MPI communication handler
    
    ! Physical properties objects (with TARGET attribute for pointer associations)
    type(chemical_properties)        ,target :: problem_chemistry      ! Chemical reaction data
    type(thermophysical_properties)  ,target :: problem_thermophysics  ! Thermodynamic properties
    
    ! Solver configuration
    type(solver_options)                     :: problem_solver_options ! Numerical method settings
    type(flame_stabilization_control)        :: problem_flame_stabilization
    type(problem_controls)                   :: problem_controls_setup
    
    ! Solution fields (with TARGET attribute)
    type(computational_mesh)         ,target :: problem_mesh          ! Computational grid
    type(boundary_conditions)        ,target :: problem_boundaries    ! Boundary conditions
    type(field_scalar_cons)          ,target :: p, T, rho, E_f, v_s   ! Scalar fields:
                                                                      ! p - pressure
                                                                      ! T - temperature
                                                                      ! rho - density
                                                                      ! E_f - flame energy?
                                                                      ! v_s - speed of sound
    type(field_vector_cons)          ,target :: v, Y                  ! Vector fields:
                                                                      ! v - velocity (3 components)
                                                                      ! Y - species mass fractions
    
    ! Post-processing and I/O
    type(post_processor_manager)             :: problem_post_proc_manager  ! Monitoring and analysis
    type(data_io)                            :: problem_data_io        ! Data input/output
    type(data_save)                          :: problem_data_save      ! Solution file writing
    
    !==========================================
    ! GEOMETRIC AND DOMAIN PARAMETERS
    !==========================================
    real(dp)    ,dimension(3)   :: cell_size           ! Cell dimensions (dx, dy, dz)
    integer     ,dimension(3,2) :: utter_loop          ! Global computational bounds
    integer     ,dimension(3,2) :: observation_slice   ! Region for monitoring
    integer     ,dimension(3,2) :: summation_region    ! Region for integral calculations
    integer                     :: transducer_offset   ! Offset for transducer placement
    
    !==========================================
    ! CHEMICAL AND PHYSICAL PARAMETERS
    !==========================================
    integer                     :: species_number       ! Number of chemical species
    real(dp)                    :: pipe_radius          ! Pipe radius [m]
    real(dp)                    :: pipe_length          ! Pipe length [m]
    real(dp)                    :: CFL_coeff            ! CFL stability coefficient
    real(dp)                    :: delta_r              ! Radial resolution [m]
    real(dp)                    :: delta_z              ! Axial resolution [m]
    real(dp)                    :: inlet_velocity       ! Inlet velocity [m/s]
    real(dp)                    :: inlet_temperature    ! Inlet temperature [K]
    real(dp)                    :: outlet_pressure      ! Outlet pressure [Pa]
    real(dp)                    :: wall_temperature     ! Wall temperature [K]
    
    !==========================================
    ! FILE AND DIRECTORY PATHS
    !==========================================
    character(len=10)           :: string              ! Temporary string buffer
    character(len=4096)         :: initial_work_dir    ! Initial working directory
    character(len=1024)         :: work_dir            ! Current working directory
    character(len=20)           :: solver_name         ! Solver type identifier
    character(len=20)           :: coordinate_system   ! Coordinate system type
    character(len=30)           :: mech_name           ! Chemical mechanism name
    character(len=30)           :: mech_file           ! Chemical mechanism file
    character(len=30)           :: thermo_file         ! Thermodynamic data file
    character(len=30)           :: transdata_file      ! Transport data file
    
    !==========================================
    ! CONTROL AND STATUS VARIABLES
    !==========================================
    logical :: stop_flag          ! Simulation termination flag
    integer :: i, j, k, spec      ! Loop indices
    integer :: ierr               ! Error status indicator
    integer :: io_unit            ! File I/O unit
    integer :: log_unit           ! Log file unit
    
    !==========================================
    ! PROGRAM EXECUTION BEGINS
    !==========================================
    
    ! Get initial working directory for later return
    call get_current_directory(initial_work_dir)
    
    !================================================================
    ! SINGLE CASE SETUP (No parametric loops for this test)
    !================================================================
    
    ! Initialize working directory structure
    work_dir = 'cylindrical_pipe_flow'  ! Main results directory
    
    ! Create directory tree
    call ensure_directory(work_dir)
    
    !------------------------------------------------
    ! GEOMETRY PARAMETERS
    !------------------------------------------------
    pipe_radius     = 0.01_dp       ! 1 cm pipe radius
    pipe_length     = 0.1_dp        ! 10 cm pipe length
    delta_r         = 5.0e-04_dp    ! 0.5 mm radial resolution
    delta_z         = 2.0e-03_dp    ! 2.0 mm axial resolution
    
    !------------------------------------------------
    ! FLOW PARAMETERS
    !------------------------------------------------
    inlet_velocity   = 10.0_dp      ! 10 m/s inlet velocity
    inlet_temperature = 300.0_dp    ! 300 K inlet temperature
    outlet_pressure  = 101325.0_dp  ! Atmospheric pressure at outlet
    wall_temperature = 300.0_dp     ! Isothermal wall at 300 K
    
    !------------------------------------------------
    ! COORDINATE SYSTEM AND SOLVER
    !------------------------------------------------
    coordinate_system = 'cylindrical'
    solver_name = 'fds_low_mach'
    
    ! Add subdirectories for organization
    work_dir = trim(work_dir) // trim(fold_sep) // trim(coordinate_system)
    call ensure_directory(work_dir)
    
    work_dir = trim(work_dir) // trim(fold_sep) // trim(solver_name)
    call ensure_directory(work_dir)
    
    ! Copy setup files to working directory and change to it.
    call copy_directory_tree(trim(task_setup_folder), &
                             trim(work_dir) // trim(fold_sep) // trim(task_setup_folder))
    call change_directory(work_dir)
    
    ! Open log file for this configuration
    open(newunit = log_unit, file = problem_setup_log_file, status = 'replace', form = 'formatted')
    call write_nrg_source_revision(log_unit, 'problem setup generation')
    
    !================================================================
    ! DOMAIN DEFINITION (2D CYLINDRICAL)
    !================================================================
    ! In cylindrical coordinates:
    ! - Dimension 1 (x) -> radial direction r
    ! - Dimension 2 (y) -> axial direction z
    ! - Dimension 3 (z) -> azimuthal (not used in 2D axisymmetric)
    
    problem_domain = computational_domain_c(  &
        dimensions         = 2,                                     &
        cells_number       = (/int(pipe_radius/delta_r), int(pipe_length/delta_z), 1/), &
        coordinate_system  = coordinate_system,                     &
        lengths            = reshape((/0.0_dp, 0.0_dp, 0.0_dp,      &
                                      pipe_radius, pipe_length, 0.0_dp/),(/3,2/)), &
        axis_names         = (/'r','z','theta'/) )
    
    !================================================================
    ! CHEMICAL AND THERMOPHYSICAL PROPERTIES INITIALIZATION
    !================================================================
    ! Using air as working fluid (N2-O2 mixture)
    mech_name      = 'AIR_SIMPLIFIED'
    mech_file      = 'AIR.txt'
    thermo_file    = 'AIR_THERMO.txt'
    transdata_file = 'AIR_TRANSDATA.txt'
    
    problem_chemistry = chemical_properties_c( &
        chemical_mechanism_file_name     = mech_file,        &
        default_enhanced_efficiencies    = 1.0_dp,           &
        E_act_units                      = 'cal.mol')
    
    problem_thermophysics = thermophysical_properties_c( &
        chemistry                    = problem_chemistry,   &
        thermo_data_file_name        = thermo_file,        &
        transport_data_file_name     = transdata_file,     &
        molar_masses_data_file_name  = 'molar_masses.dat')
    
    !================================================================
    ! SOLVER OPTIONS CONFIGURATION
    !================================================================
    problem_solver_options = solver_options_c( &
        solver_name                 = solver_name,                              &
        hydrodynamics_flag          = .true., &      ! Solve momentum equations
        heat_transfer_flag          = .true., &      ! Solve energy equation
        molecular_diffusion_flag    = .false., &     ! No species diffusion for single fluid
        viscosity_flag              = .true., &      ! Include viscous effects
        thermal_radiation_flag      = .false., &     ! No radiation
        chemical_reaction_flag      = .false., &     ! No chemistry for pure air
        grav_acc                    = (/0.0_dp, 0.0_dp, 0.0_dp/), &  ! No gravity
        additional_particles_phases = 0, &           ! No particle phases
        CFL_flag                    = .true.,  &     ! Use CFL condition
        CFL_coefficient             = 0.5_dp,  &     ! CFL safety factor
        initial_time_step           = 1e-07_dp)      ! Initial Δt [s]
    
    ! No flame stabilization for pipe flow
    problem_flame_stabilization = flame_stabilization_control_c( &
        mode = 'none')
    
    problem_controls_setup = problem_controls_c( &
        flame_stabilization = problem_flame_stabilization)
    
    !================================================================
    ! MPI AND DATA MANAGEMENT SETUP
    !================================================================
    problem_mpi_support   = mpi_communications_c(problem_domain)
    problem_data_manager  = data_manager_c(problem_domain, problem_mpi_support, &
                                           problem_chemistry, problem_thermophysics, &
                                           problem_solver_options, problem_controls_setup)
    call problem_controls_setup%write_log(log_unit)
    
    ! Create boundary conditions (4 types: symmetry, wall, inlet, outlet)
    call problem_data_manager%create_boundary_conditions( &
        problem_boundaries, number_of_boundary_types = 4, default_boundary = 1)
    
    ! Create computational mesh
    call problem_data_manager%create_computational_mesh(problem_mesh)
    
    ! Create scalar solution fields
    call problem_data_manager%create_scalar_field(p,   'pressure',    'p')
    call problem_data_manager%create_scalar_field(T,   'temperature', 'T')
    call problem_data_manager%create_scalar_field(rho, 'density',     'rho')
    
    ! Create vector solution fields
    call problem_data_manager%create_vector_field(v, 'velocity', 'v', 'spatial')
    call problem_data_manager%create_vector_field(Y, 'specie_mass_fraction', 'Y', 'chemical')
    
    ! Get geometric information
    cell_size = problem_mesh%get_cell_edges_length()
    utter_loop = problem_domain%get_global_utter_cells_bounds()
    
    ! Define monitoring regions (centerline and near-wall)
    observation_slice(:,1) = (/1, 1, 1/)
    observation_slice(:,2) = utter_loop(:,2)
    
    summation_region(:,1) = (/1, 1, 1/)
    summation_region(:,2) = utter_loop(:,2)
    
    !================================================================
    ! POST-PROCESSING SETUP
    !================================================================
    problem_post_proc_manager = post_processor_manager_c(problem_data_manager, number_post_processors = 2)
  
    ! Post-processor 1: Flow diagnostics along centerline
    call problem_post_proc_manager%create_post_processor( &
        problem_data_manager,                            &
        post_processor_name  = 'centerline_monitor',     &
        operations_number    = 3,                        &
        save_time            = 1.0_dp,                   &
        save_time_units      = 'milliseconds',           &
        post_processor_title = 'Centerline flow diagnostics')
    
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager, 1, 'velocity', 'max', &
        operation_area = observation_slice, grad_projection = 2)
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager, 1, 'pressure', 'average', &
        operation_area = observation_slice)
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager, 1, 'temperature', 'average', &
        operation_area = observation_slice)
    
    ! Post-processor 2: Wall shear stress monitoring
    call problem_post_proc_manager%create_post_processor( &
        problem_data_manager,                            &
        post_processor_name  = 'wall_monitor',           &
        operations_number    = 2,                        &
        save_time            = 1.0_dp,                   &
        save_time_units      = 'milliseconds',           &
        post_processor_title = 'Wall boundary diagnostics')
    
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager, 2, 'velocity', 'min', &
        operation_area = observation_slice, grad_projection = 1)
    call problem_post_proc_manager%create_post_processor_operation( &
        problem_data_manager, 2, 'pressure', 'transducer', &
        operation_area_distance = (/0.0_dp, pipe_length*0.5_dp, 0.0_dp/))
    
    !================================================================
    ! DATA SAVING CONFIGURATION
    !================================================================
    problem_data_save = data_save_c( &
        problem_data_manager, &
        visible_fields_names = [ character(len=40) :: &
            'pressure',                       &
            'pressure_dynamic',               &
            'temperature',                    &
            'density',                        &
            'velocity',                       &
            'velocity_magnitude',             &
            'vorticity',                      &
            'viscosity',                      &
            'mach_number'                     &
        ], &
        save_time         = 5.0_dp,          &   ! Save interval
        save_time_units   = 'milliseconds',  &   ! Time units for saving
        save_format       = 'tecplot',       &   ! Output format
        data_save_folder  = 'data_save',     &   ! Output directory
        dataset_name      = 'Cylindrical_Pipe_Flow', &
        debug_flag        = .false.)              ! Debug mode off
    
    !================================================================
    ! DATA I/O CONFIGURATION
    !================================================================
    problem_data_io = data_io_c( &
        problem_data_manager,    &
        check_time         = 10.0_dp,        &  ! Checkpoint interval
        check_time_units   = 'milliseconds', &  ! Time units for checkpoints
        data_output_folder = 'data_output')     ! Checkpoint directory
    
    !================================================================
    ! INITIAL CONDITIONS SETUP
    !================================================================
    
    ! Set uniform ambient conditions throughout the domain
    T%cells(:,:,:)   = inlet_temperature
    p%cells(:,:,:)   = outlet_pressure
    rho%cells(:,:,:) = 1.225_dp  ! Air density at STP
    
    ! Initialize velocity field to zero (fluid at rest initially)
    v%pr(1)%cells(:,:,:) = 0.0_dp  ! Radial velocity
    v%pr(2)%cells(:,:,:) = 0.0_dp  ! Axial velocity
    v%pr(3)%cells(:,:,:) = 0.0_dp  ! Azimuthal velocity
    
    ! Initialize species (pure air: 23.3% O2, 76.7% N2 by mass)
    if (associated(Y)) then
        do spec = 1, problem_chemistry%species_number
            Y%pr(spec)%cells(:,:,:) = 0.0_dp
        end do
        ! Set O2 and N2 mass fractions (assuming species order from mechanism)
        ! Adjust indices based on actual mechanism file
        do i = utter_loop(1,1), utter_loop(1,2)
            do j = utter_loop(2,1), utter_loop(2,2)
                do k = utter_loop(3,1), utter_loop(3,2)
                    ! Find species indices dynamically if needed
                    ! For now, set dummy values
                end do
            end do
        end do
    end if
    
    !================================================================
    ! BOUNDARY CONDITIONS SETUP
    !================================================================
    ! Boundary markers in 2D cylindrical (r-z):
    ! - Left (r=0): symmetry axis (marker 1)
    ! - Right (r=R): pipe wall (marker 2)
    ! - Bottom (z=0): inlet (marker 3)
    ! - Top (z=L): outlet (marker 4)
    
    ! Boundary type 1: Symmetry plane at pipe axis (r=0)
    call problem_boundaries%create_boundary_type( &
        type_name               = 'symmetry_plane',    &
        priority                = 1)
    
    ! Boundary type 2: No-slip isothermal wall at pipe wall (r=R)
    call problem_boundaries%create_boundary_type( &
        type_name               = 'wall',              &
        slip                    = .false.,             &  ! No-slip condition
        conductive              = .true.,              &  ! Heat transfer enabled
        wall_temperature        = wall_temperature,    &  ! Fixed wall temperature
        wall_conductivity_ratio = 1.0_dp,              &  ! Same conductivity as fluid
        priority                = 2)
    
    ! Boundary type 3: Velocity inlet at pipe entrance (z=0)
    call problem_boundaries%create_boundary_type( &
        type_name               = 'inlet',             &
        farfield_pressure       = outlet_pressure,     &  ! Reference pressure
        farfield_temperature    = inlet_temperature,   &  ! Inlet temperature
        farfield_velocity       = (/0.0_dp, inlet_velocity, 0.0_dp/), &  ! Axial velocity only
        farfield_species_names  = [character(len=5) :: 'O2','N2'], &
        farfield_concentrations = (/0.233_dp, 0.767_dp/), &  ! Air composition by mass
        priority                = 3)
    
    ! Boundary type 4: Pressure outlet at pipe exit (z=L)
    call problem_boundaries%create_boundary_type( &
        type_name               = 'outlet',            &
        farfield_pressure       = outlet_pressure,     &  ! Fixed outlet pressure
        farfield_temperature    = inlet_temperature,   &  ! Reference temperature
        farfield_velocity       = (/0.0_dp, 0.0_dp, 0.0_dp/), &  ! Extrapolated velocity
        farfield_species_names  = [character(len=5) :: 'O2','N2'], &
        farfield_concentrations = (/0.233_dp, 0.767_dp/), &
        priority                = 4)
    
    ! Apply boundary markers to domain boundaries
    ! Left boundary (r=0): symmetry axis
    problem_boundaries%bc_markers(utter_loop(1,1),:,:) = 1
    
    ! Right boundary (r=R): pipe wall
    problem_boundaries%bc_markers(utter_loop(1,2),:,:) = 2
    
    ! Bottom boundary (z=0): inlet
    problem_boundaries%bc_markers(:,utter_loop(2,1),:) = 3
    
    ! Top boundary (z=L): outlet
    problem_boundaries%bc_markers(:,utter_loop(2,2),:) = 4
    
    !================================================================
    ! LOG FILE ENTRY
    !================================================================
    write(log_unit,'(A)') 'General description: '
    write(log_unit,'(A)') '  2D axisymmetric simulation of gas flow in a cylindrical pipe'
    write(log_unit,'(A)') 'Main aim: '
    write(log_unit,'(A)') '  Study laminar/turbulent pipe flow development, pressure drop, and boundary layers'
    write(log_unit,'(A)') 'Problem setup:'
    write(log_unit,'(A,I0,A,F8.4,A)') '  Pipe radius: ', nint(pipe_radius*1000), ' mm (', pipe_radius, ' m)'
    write(log_unit,'(A,I0,A,F8.4,A)') '  Pipe length: ', nint(pipe_length*1000), ' mm (', pipe_length, ' m)'
    write(log_unit,'(A,F6.3,A)') '  Radial resolution: ', delta_r*1000, ' mm'
    write(log_unit,'(A,F6.3,A)') '  Axial resolution: ', delta_z*1000, ' mm'
    write(log_unit,'(A,I0,A,I0,A)') '  Grid size: ', utter_loop(1,2)-utter_loop(1,1)+1, ' x ', &
                                    utter_loop(2,2)-utter_loop(2,1)+1, ' cells'
    write(log_unit,'(A)') 'Boundary conditions:'
    write(log_unit,'(A)') '  r=0: Symmetry axis'
    write(log_unit,'(A)') '  r=R: No-slip isothermal wall'
    write(log_unit,'(A,F6.1,A)') '  z=0: Velocity inlet (', inlet_velocity, ' m/s)'
    write(log_unit,'(A,F8.1,A)') '  z=L: Pressure outlet (', outlet_pressure/101325.0_dp, ' atm)'
    write(log_unit,'(A)') 'Solver setup: '
    write(log_unit,'(A)') '  FDS low-Mach number solver with cylindrical coordinates'
    write(log_unit,'(A)') 'Validation and comparison:'
    write(log_unit,'(A)') '  Compare with analytical Poiseuille flow solution for laminar regime'
    write(log_unit,'(A)') '----------------------------------------------------------'
    
    close(log_unit)
    
    !================================================================
    ! INITIAL OUTPUT AND CLEANUP
    !================================================================
    ! Output initial solution and save data
    call problem_data_io%output_all_data(0.0_dp, stop_flag, make_output = .true.)
    call problem_data_save%save_all_data(0.0_dp, stop_flag, make_save   = .true.)
    
    ! Return to initial directory for next configuration
    call change_directory(initial_work_dir)
    
    print *, '>>> Case generation complete.'
    print *, '>>> Case directory: ', trim(work_dir)
    print *, '>>> To run the simulation:'
    print *, '   cd ', trim(work_dir)
    print *, '   ../computing_module.exe --num_threads=4'
    
end program package_interface

!================================================================
! CONTAINED PROCEDURES
!================================================================

contains

    !--------------------------------------------------------------------------
    ! Return the process working directory without relying on compiler-specific
    ! modules such as IFPORT.
    !--------------------------------------------------------------------------
    subroutine get_current_directory(directory)
        character(len=*), intent(out) :: directory
        character(kind=c_char) :: buffer(4096)
        type(c_ptr) :: ptr
        integer :: idx, max_copy

        buffer = c_null_char
#ifdef WIN
        ptr = c_getcwd(buffer, int(size(buffer), c_int))
#else
        ptr = c_getcwd(buffer, int(size(buffer), c_size_t))
#endif
        if (.not. c_associated(ptr)) then
            error stop 'Unable to obtain current working directory.'
        end if

        directory = ''
        max_copy = min(len(directory), size(buffer))
        do idx = 1, max_copy
            if (buffer(idx) == c_null_char) exit
            directory(idx:idx) = transfer(buffer(idx), directory(idx:idx))
        end do

        if (idx > len(directory) .and. buffer(min(idx,size(buffer))) /= c_null_char) then
            error stop 'Current working directory path exceeds internal buffer.'
        end if
    end subroutine get_current_directory

    !--------------------------------------------------------------------------
    ! Change the process working directory through the platform C runtime.
    !--------------------------------------------------------------------------
    subroutine change_directory(directory)
        character(len=*), intent(in) :: directory
        character(kind=c_char), allocatable :: c_path(:)
        integer(c_int) :: status
        integer :: idx, path_length

        path_length = len_trim(directory)
        allocate(c_path(path_length + 1))
        do idx = 1, path_length
            c_path(idx) = transfer(directory(idx:idx), c_path(idx))
        end do
        c_path(path_length + 1) = c_null_char

        status = c_chdir(c_path)
        if (status /= 0_c_int) then
            write(*,'(A)') 'Unable to change working directory to: ' // trim(directory)
            error stop 'Directory change failed.'
        end if
    end subroutine change_directory

    !--------------------------------------------------------------------------
    ! Run a CMake -E filesystem command and fail immediately on errors.
    !--------------------------------------------------------------------------
    subroutine run_filesystem_command(command)
        character(len=*), intent(in) :: command
        integer :: cmdstat, exitstat

        cmdstat = 0
        exitstat = 0
        call execute_command_line(command, wait=.true., exitstat=exitstat, cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat /= 0) then
            write(*,'(A)') 'Filesystem command failed: ' // trim(command)
            write(*,'(A,I0,A,I0)') 'cmdstat=', cmdstat, ', exitstat=', exitstat
            error stop 'Case-directory setup failed.'
        end if
    end subroutine run_filesystem_command

    subroutine ensure_directory(directory)
        character(len=*), intent(in) :: directory
        character(len=:), allocatable :: command

        command = 'cmake -E make_directory "' // trim(directory) // '"'
        call run_filesystem_command(command)
    end subroutine ensure_directory

    subroutine copy_directory_tree(source_directory, destination_directory)
        character(len=*), intent(in) :: source_directory
        character(len=*), intent(in) :: destination_directory
        character(len=:), allocatable :: command

        command = 'cmake -E copy_directory "' // trim(source_directory) // '" "' // &
                  trim(destination_directory) // '"'
        call run_filesystem_command(command)
    end subroutine copy_directory_tree

end program package_interface
!================================================================================
! END OF PROGRAM
!================================================================================
