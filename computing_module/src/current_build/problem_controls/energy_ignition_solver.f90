module energy_ignition_solver_class

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    use kind_parameters, only: dp
    use global_data, only: pi
    use data_manager_class
    use field_pointers
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use mpi_communications_class
    use problem_control_class, only: energy_ignition_control

    implicit none

    private
    public :: energy_ignition_solver, energy_ignition_solver_c

    type :: energy_ignition_solver
        private
        logical :: enabled = .false.
        type(energy_ignition_control) :: control
        type(field_scalar_cons_pointer) :: energy_source
        type(field_scalar_cons), pointer :: energy_source_store => null()
        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary
        type(mpi_communications) :: mpi_support
        real(dp) :: source_region_volume = 0.0_dp
    contains
        procedure :: solve
        procedure :: is_enabled
        procedure, private :: cell_is_in_region
        procedure, private :: physical_cell_volume
    end type energy_ignition_solver

    interface energy_ignition_solver_c
        module procedure constructor
    end interface energy_ignition_solver_c

contains

    type(energy_ignition_solver) function constructor(manager)
        type(data_manager), intent(inout) :: manager

        integer :: i, j, k
        integer, dimension(3,2) :: cell_loop
        real(dp) :: local_volume

        constructor%control = manager%problem_controls_config%get_energy_ignition()
        constructor%enabled = constructor%control%is_enabled()
        constructor%domain = manager%domain
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%mpi_support = manager%mpi_communications

        allocate(constructor%energy_source_store)
        call manager%create_scalar_field( &
            constructor%energy_source_store, &
            'energy_production_ignition', 'E_f_prod_ignition')
        constructor%energy_source%s_ptr => constructor%energy_source_store
        constructor%energy_source%s_ptr%cells = 0.0_dp

        if (.not. constructor%enabled) return

        cell_loop = constructor%domain%get_local_inner_cells_bounds()
        local_volume = 0.0_dp
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (constructor%boundary%bc_ptr%bc_markers(i,j,k) /= 0) cycle
                    if (.not. constructor%cell_is_in_region(i,j,k)) cycle
                    local_volume = local_volume + &
                        constructor%physical_cell_volume(i,j,k)
                end do
            end do
        end do

        call constructor%mpi_support%global_sum_real( &
            local_volume, constructor%source_region_volume)
        if (.not. ieee_is_finite(constructor%source_region_volume) .or. &
            constructor%source_region_volume <= 0.0_dp) then
            error stop 'Energy ignition: configured source region contains no active cells'
        end if
    end function constructor


    subroutine solve(this, step_start_time, time_step)
        class(energy_ignition_solver), intent(inout) :: this
        real(dp), intent(in) :: step_start_time, time_step

        integer :: i, j, k
        integer, dimension(3,2) :: cell_loop
        real(dp) :: control_start, control_end, step_end
        real(dp) :: overlap, average_power, source_value

        this%energy_source%s_ptr%cells = 0.0_dp
        if (.not. this%enabled) return

        if (.not. ieee_is_finite(time_step) .or. time_step <= 0.0_dp) then
            error stop 'Energy ignition: time step must be finite and positive'
        end if

        control_start = this%control%get_start_time()
        control_end = control_start + this%control%get_duration()
        step_end = step_start_time + time_step
        overlap = max(0.0_dp, min(step_end, control_end) - &
            max(step_start_time, control_start))
        if (overlap <= 0.0_dp) return

        ! Average the source over the current numerical time step.  This makes
        ! the integrated deposited energy independent of whether a step crosses
        ! the start or end of the prescribed ignition interval.
        average_power = this%control%get_total_energy() / &
            this%control%get_duration() * overlap / time_step
        source_value = average_power / this%source_region_volume

        cell_loop = this%domain%get_local_inner_cells_bounds()
        do k = cell_loop(3,1), cell_loop(3,2)
            do j = cell_loop(2,1), cell_loop(2,2)
                do i = cell_loop(1,1), cell_loop(1,2)
                    if (this%boundary%bc_ptr%bc_markers(i,j,k) /= 0) cycle
                    if (this%cell_is_in_region(i,j,k)) then
                        this%energy_source%s_ptr%cells(i,j,k) = source_value
                    end if
                end do
            end do
        end do
    end subroutine solve


    logical function is_enabled(this)
        class(energy_ignition_solver), intent(in) :: this
        is_enabled = this%enabled
    end function is_enabled


    logical function cell_is_in_region(this, i, j, k)
        class(energy_ignition_solver), intent(in) :: this
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


    real(dp) function physical_cell_volume(this, i, j, k)
        class(energy_ignition_solver), intent(in) :: this
        integer, intent(in) :: i, j, k

        integer :: dimensions
        real(dp) :: base_volume, radius
        character(len=20) :: coordinate_system

        dimensions = this%domain%get_domain_dimensions()
        coordinate_system = this%domain%get_coordinate_system_name()
        base_volume = this%mesh%mesh_ptr%get_cell_volume()
        radius = this%mesh%mesh_ptr%mesh(1,i,j,k)

        select case (trim(coordinate_system))
        case ('cartesian')
            physical_cell_volume = base_volume
        case ('cylindrical')
            ! Axisymmetric NRG coordinates omit the azimuthal coordinate.
            physical_cell_volume = 2.0_dp*pi*max(radius,0.0_dp)*base_volume
        case ('spherical')
            if (dimensions /= 1) then
                error stop 'Energy ignition: spherical total-energy normalization currently requires 1D radial geometry'
            end if
            physical_cell_volume = 4.0_dp*pi*max(radius,0.0_dp)**2*base_volume
        case default
            error stop 'Energy ignition: unsupported coordinate system'
        end select
    end function physical_cell_volume

end module energy_ignition_solver_class
