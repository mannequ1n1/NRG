module turbulence_forcing_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: int64

    use kind_parameters, only: dp
    use global_data, only: pi
    use data_manager_class
    use field_pointers
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use mpi_communications_class
    use problem_control_class, only: turbulence_forcing_control

    implicit none

    private
    public :: turbulence_forcing_solver, turbulence_forcing_solver_c

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 48271_int64

    type :: turbulence_forcing_solver
        private
        logical :: enabled = .false.
        logical :: forcing_state_initialized = .false.
        type(turbulence_forcing_control) :: control
        type(field_vector_cons_pointer) :: acceleration_source
        type(field_vector_cons), pointer :: acceleration_source_store => null()
        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary
        type(mpi_communications) :: mpi_support
        integer(int64) :: rng_state = 1_int64
        real(dp) :: forcing_angle = 0.0_dp
        real(dp) :: forcing_phase = 0.0_dp
        real(dp) :: next_update_time = 0.0_dp
    contains
        procedure :: solve
        procedure :: is_enabled
        procedure, private :: update_forcing_state
        procedure, private :: random_uniform
        procedure, private :: cell_is_in_region
    end type turbulence_forcing_solver

    interface turbulence_forcing_solver_c
        module procedure constructor
    end interface turbulence_forcing_solver_c

contains

    type(turbulence_forcing_solver) function constructor(manager)
        type(data_manager), intent(inout) :: manager

        constructor%control = manager%problem_controls_config%get_turbulence_forcing()
        constructor%enabled = constructor%control%is_enabled()
        constructor%domain = manager%domain
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%mpi_support = manager%mpi_communications

        allocate(constructor%acceleration_source_store)
        call manager%create_vector_field( &
            constructor%acceleration_source_store, &
            'velocity_production_turbulence', 'v_prod_turbulence', 'spatial')
        constructor%acceleration_source%v_ptr => &
            constructor%acceleration_source_store
        call zero_acceleration(constructor%acceleration_source%v_ptr)

        if (.not. constructor%enabled) return
        call constructor%control%validate_problem_compatibility( &
            constructor%domain%get_domain_dimensions())
        constructor%rng_state = int(constructor%control%get_random_seed(), int64)
        constructor%rng_state = modulo(constructor%rng_state, rng_modulus - 1_int64) + 1_int64
        constructor%next_update_time = constructor%control%get_start_time()
    end function constructor


    subroutine solve(this, step_start_time, time_step)
        class(turbulence_forcing_solver), intent(inout) :: this
        real(dp), intent(in) :: step_start_time, time_step

        integer :: i, j, k, dimensions
        integer, dimension(3,2) :: cell_loop
        real(dp) :: control_start, control_end, step_end, overlap
        real(dp) :: sample_time, ramp_factor, amplitude, wave_number
        real(dp) :: x, y, phase_argument, radius

        call zero_acceleration(this%acceleration_source%v_ptr)
        if (.not. this%enabled) return

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Turbulence forcing: time step must be finite and positive'
        end if

        control_start = this%control%get_start_time()
        control_end = control_start + this%control%get_duration()
        step_end = step_start_time + time_step
        overlap = max(0.0_dp, min(step_end, control_end) - &
            max(step_start_time, control_start))
        if (overlap <= 0.0_dp) return

        sample_time = max(control_start, min(0.5_dp*(step_start_time + step_end), control_end))
        call this%update_forcing_state(sample_time)

        if (this%control%get_ramp_time() > 0.0_dp) then
            ramp_factor = min(max((sample_time-control_start) / &
                this%control%get_ramp_time(), 0.0_dp), 1.0_dp)
        else
            ramp_factor = 1.0_dp
        end if

        amplitude = this%control%get_acceleration_amplitude() * &
            ramp_factor * overlap / time_step
        wave_number = 2.0_dp*pi / this%control%get_integral_scale()
        radius = this%control%get_radius()
        dimensions = this%domain%get_domain_dimensions()
        cell_loop = this%domain%get_local_inner_cells_bounds()

        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (this%boundary%bc_ptr%bc_markers(i,j,k) /= 0) cycle
                    if (radius > 0.0_dp .and. .not. this%cell_is_in_region(i,j,k)) cycle

                    x = this%mesh%mesh_ptr%mesh(1,i,j,k)
                    y = this%mesh%mesh_ptr%mesh(2,i,j,k)
                    phase_argument = wave_number * ( &
                        sin(this%forcing_angle)*x + &
                        cos(this%forcing_angle)*y) + this%forcing_phase

                    ! A vector perpendicular to the planar wave vector gives a
                    ! divergence-free single-mode forcing field, preserving the
                    ! useful structure of the former FDS-local perturbation.
                    this%acceleration_source%v_ptr%pr(1)%cells(i,j,k) = &
                        amplitude*cos(this%forcing_angle)*cos(phase_argument)
                    this%acceleration_source%v_ptr%pr(2)%cells(i,j,k) = &
                        -amplitude*sin(this%forcing_angle)*cos(phase_argument)
                    if (dimensions >= 3) &
                        this%acceleration_source%v_ptr%pr(3)%cells(i,j,k) = 0.0_dp
                end do
            end do
        end do

        call this%mpi_support%exchange_conservative_vector_field( &
            this%acceleration_source%v_ptr)
    end subroutine solve


    logical function is_enabled(this)
        class(turbulence_forcing_solver), intent(in) :: this
        is_enabled = this%enabled
    end function is_enabled


    subroutine update_forcing_state(this, sample_time)
        class(turbulence_forcing_solver), intent(inout) :: this
        real(dp), intent(in) :: sample_time

        real(dp) :: interval

        interval = this%control%get_update_interval()
        if (.not. this%forcing_state_initialized) then
            this%forcing_angle = 2.0_dp*pi*this%random_uniform()
            this%forcing_phase = 2.0_dp*pi*this%random_uniform()
            this%forcing_state_initialized = .true.
            this%next_update_time = this%control%get_start_time() + interval
        end if

        do while (sample_time >= this%next_update_time)
            this%forcing_angle = 2.0_dp*pi*this%random_uniform()
            this%forcing_phase = 2.0_dp*pi*this%random_uniform()
            this%next_update_time = this%next_update_time + interval
        end do
    end subroutine update_forcing_state


    real(dp) function random_uniform(this)
        class(turbulence_forcing_solver), intent(inout) :: this

        this%rng_state = modulo(rng_multiplier*this%rng_state, rng_modulus)
        if (this%rng_state <= 0_int64) this%rng_state = 1_int64
        random_uniform = real(this%rng_state,dp) / real(rng_modulus,dp)
    end function random_uniform


    logical function cell_is_in_region(this, i, j, k)
        class(turbulence_forcing_solver), intent(in) :: this
        integer, intent(in) :: i, j, k

        integer :: dimensions, dim
        real(dp), dimension(3) :: center
        real(dp) :: distance_squared, radius

        dimensions = this%domain%get_domain_dimensions()
        center = this%control%get_center()
        radius = this%control%get_radius()
        distance_squared = 0.0_dp
        do dim = 1, dimensions
            distance_squared = distance_squared + &
                (this%mesh%mesh_ptr%mesh(dim,i,j,k) - center(dim))**2
        end do
        cell_is_in_region = distance_squared <= radius*radius
    end function cell_is_in_region


    subroutine zero_acceleration(field)
        type(field_vector_cons), intent(inout) :: field
        integer :: dim

        do dim = 1, size(field%pr)
            field%pr(dim)%cells = 0.0_dp
        end do
    end subroutine zero_acceleration

end module turbulence_forcing_solver_class
