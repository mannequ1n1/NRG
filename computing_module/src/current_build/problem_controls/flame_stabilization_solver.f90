module flame_stabilization_solver_class

    use kind_parameters, only: dp
    use data_manager_class
    use field_pointers
    use computational_domain_class
    use computational_mesh_class
    use boundary_conditions_class
    use chemical_properties_class
    use problem_control_class, only: flame_stabilization_control
    use supplementary_routines

    implicit none

    private
    public :: flame_stabilization_solver, flame_stabilization_solver_c


    type :: flame_stabilization_runtime_state
        real(dp), allocatable :: time_hist(:), front_coord_hist(:)
        real(dp), allocatable :: diag_time_hist(:), diag_front_coord_hist(:)
        character(len=200) :: data_table_filename = ''
        character(len=200) :: chem_table_filename = ''
        integer :: track_counter = 0
        integer :: correction_counter = 0
        integer :: stabilization_counter = 0
        integer :: hist_count = 0
        integer :: diag_hist_count = 0
        integer :: same_sign_error_counter = 0
        integer :: post_flamelet_hold_counter = 0
        integer :: measurement_attempt = 0
        integer :: flame_loc_unit = -1
        real(dp) :: previous_correction_time = -huge(1.0_dp)
        real(dp) :: filtered_velocity_save = 0.0_dp
        real(dp) :: diag_filtered_velocity_save = 0.0_dp
        real(dp) :: adaptive_gain = 5.0e-02_dp
        real(dp) :: inlet_velocity_target = 0.0_dp
        real(dp) :: inlet_velocity_applied = 0.0_dp
        real(dp) :: ramp_start_time = 0.0_dp
        real(dp) :: ramp_start_velocity = 0.0_dp
        real(dp) :: active_inlet_ramp_time = 2.0e-03_dp
        real(dp) :: previous_control_velocity = 0.0_dp
        real(dp) :: previous_control_inlet_velocity = 0.0_dp
        real(dp) :: front_reference_coord = 0.0_dp
        real(dp) :: heat_release_peak_save = 0.0_dp
        real(dp) :: bracket_u_a = 0.0_dp
        real(dp) :: bracket_v_a = 0.0_dp
        real(dp) :: bracket_u_b = 0.0_dp
        real(dp) :: bracket_v_b = 0.0_dp
        logical :: initialized = .false.
        logical :: inlet_velocity_initialized = .false.
        logical :: output_initialized = .false.
        logical :: flamelet_output_written = .false.
        logical :: sl_output_written = .false.
        logical :: have_previous_control_point = .false.
        logical :: has_bracket = .false.
        logical :: front_reference_initialized = .false.
        integer :: control_stage = 0
        real(dp) :: stabilized_inlet_velocity = 0.0_dp
        real(dp) :: measurement_inlet_velocity = 0.0_dp
        real(dp) :: measurement_start_time = 0.0_dp
        real(dp) :: measurement_start_coord = 0.0_dp
        real(dp) :: sl_displacement_save = 0.0_dp
        real(dp) :: measurement_velocity_save = 0.0_dp
        real(dp) :: measurement_r2_save = 0.0_dp
        real(dp) :: measurement_rms_save = 0.0_dp
        real(dp) :: measurement_split_slope_diff_save = 0.0_dp
        real(dp) :: current_measurement_delta = 0.0_dp
    end type flame_stabilization_runtime_state

    type :: flame_stabilization_solver
        private
        logical :: enabled = .false.
        type(flame_stabilization_control) :: flame_stabilization

        type(computational_domain) :: domain
        type(computational_mesh_pointer) :: mesh
        type(boundary_conditions_pointer) :: boundary
        type(chemical_properties_pointer) :: chem
        type(field_scalar_cons_pointer) :: T
        type(field_scalar_cons_pointer) :: E_f_prod_chem
        type(field_vector_cons_pointer) :: Y

        integer :: load_counter = 0
        real(dp) :: inlet_velocity = 0.0_dp

        logical :: flamelet_output_pending = .false.
        character(len=200) :: requested_data_table_filename = ''
        character(len=200) :: requested_chem_table_filename = ''
        type(flame_stabilization_runtime_state) :: state
    contains
        procedure :: is_enabled => flame_stabilization_solver_is_enabled
        procedure :: solve
        procedure :: get_inlet_velocity
        procedure :: consume_flamelet_output_request
        procedure, private :: set_inlet_velocity
    end type flame_stabilization_solver

    interface flame_stabilization_solver_c
        module procedure constructor
    end interface flame_stabilization_solver_c

contains

    type(flame_stabilization_solver) function constructor(manager, load_counter, initial_inlet_velocity)
        type(data_manager), intent(inout) :: manager
        integer, intent(in) :: load_counter
        real(dp), intent(in) :: initial_inlet_velocity

        type(field_scalar_cons_pointer) :: scal_ptr
        type(field_vector_cons_pointer) :: vect_ptr
        type(field_tensor_cons_pointer) :: tens_ptr
        integer :: dimensions, number_of_boundary_types
        integer :: inlet_count, outlet_count, bound_number
        character(len=20) :: boundary_type_name

        constructor%flame_stabilization = &
            manager%problem_controls_config%get_flame_stabilization()
        constructor%enabled = constructor%flame_stabilization%is_enabled()
        constructor%load_counter = load_counter
        constructor%inlet_velocity = initial_inlet_velocity

        if (.not. constructor%enabled) return

        dimensions = manager%domain%get_domain_dimensions()
        inlet_count = 0
        outlet_count = 0
        number_of_boundary_types = &
            manager%boundary_conditions_pointer%bc_ptr%get_boundary_types()
        do bound_number = 1, number_of_boundary_types
            boundary_type_name = manager%boundary_conditions_pointer%bc_ptr% &
                boundary_types(bound_number)%get_type_name()
            select case (boundary_type_name)
                case ('inlet')
                    inlet_count = inlet_count + 1
                case ('outlet')
                    outlet_count = outlet_count + 1
            end select
        end do
        call constructor%flame_stabilization%validate_problem_compatibility( &
            manager%solver_options%get_chemical_reaction_flag(), &
            dimensions, inlet_count, outlet_count)

        constructor%domain = manager%domain
        constructor%mesh%mesh_ptr => manager%computational_mesh_pointer%mesh_ptr
        constructor%boundary%bc_ptr => manager%boundary_conditions_pointer%bc_ptr
        constructor%chem%chem_ptr => manager%chemistry%chem_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'temperature')
        constructor%T%s_ptr => scal_ptr%s_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'specie_mass_fraction')
        constructor%Y%v_ptr => vect_ptr%v_ptr

        call manager%get_cons_field_pointer_by_name( &
            scal_ptr, vect_ptr, tens_ptr, 'energy_production_chemistry')
        constructor%E_f_prod_chem%s_ptr => scal_ptr%s_ptr
    end function constructor

    logical function flame_stabilization_solver_is_enabled(this)
        class(flame_stabilization_solver), intent(in) :: this
        flame_stabilization_solver_is_enabled = this%enabled
    end function flame_stabilization_solver_is_enabled

    real(dp) function get_inlet_velocity(this)
        class(flame_stabilization_solver), intent(in) :: this
        get_inlet_velocity = this%inlet_velocity
    end function get_inlet_velocity

    subroutine consume_flamelet_output_request(this, requested, data_table_filename, chem_table_filename)
        class(flame_stabilization_solver), intent(inout) :: this
        logical, intent(out) :: requested
        character(len=*), intent(out) :: data_table_filename
        character(len=*), intent(out) :: chem_table_filename

        requested = this%flamelet_output_pending
        data_table_filename = ''
        chem_table_filename = ''

        if (requested) then
            data_table_filename = trim(this%requested_data_table_filename)
            chem_table_filename = trim(this%requested_chem_table_filename)
            this%flamelet_output_pending = .false.
        end if
    end subroutine consume_flamelet_output_request

    subroutine set_inlet_velocity(this, inlet_velocity)
        class(flame_stabilization_solver), intent(inout) :: this
        real(dp), intent(in) :: inlet_velocity

        integer :: bound_number, boundary_types

        this%inlet_velocity = inlet_velocity
        boundary_types = this%boundary%bc_ptr%get_boundary_types()
        do bound_number = 1, boundary_types
            if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
                call this%boundary%bc_ptr%boundary_types(bound_number)%set_farfield_velocity(inlet_velocity)
            end if
        end do
    end subroutine set_inlet_velocity

    !--------------------------------------------------------------------------
    ! One-dimensional flame anchoring / laminar-burning-velocity controller.
    ! The numerical algorithm is intentionally kept equivalent to the former
    ! FDS-owned flame-stabilization implementation.  This subsolver owns the
    ! controller state and boundary intervention, while the gas solver remains
    ! responsible for solver-specific flamelet/chemistry table output.
    !--------------------------------------------------------------------------
    subroutine solve(this, time, stabilized)
        class(flame_stabilization_solver), intent(inout) :: this
        real(dp), intent(in) :: time
        logical, intent(out) :: stabilized

        integer, parameter :: max_hist_size = 120
        integer, parameter :: max_diag_hist_size = 200
        integer, parameter :: min_hist_for_control = 12
        integer, parameter :: min_hist_for_capture = 6
        integer, parameter :: stable_required_count = 50
        integer, parameter :: flamelet_required_count = 50
        integer, parameter :: post_flamelet_hold_required_count = 100

        real(dp), parameter :: time_control_capture = 2.0e-04_dp
        real(dp), parameter :: response_settle_capture = 2.0e-04_dp
        real(dp), parameter :: inlet_ramp_time_fine = 2.0e-03_dp
        real(dp), parameter :: inlet_ramp_time_capture = 1.0e-03_dp
        real(dp), parameter :: controller_gain_initial = 2.5e-01_dp
        real(dp), parameter :: controller_gain_min = 2.0e-02_dp
        real(dp), parameter :: controller_gain_max = 5.0e-01_dp
        real(dp), parameter :: controller_gain_capture = 5.0e-01_dp
        real(dp), parameter :: controller_max_fraction = 5.0e-02_dp
        real(dp), parameter :: controller_max_fraction_capture = 1.0e-01_dp
        real(dp), parameter :: controller_max_fraction_emergency = 2.5e-01_dp
        real(dp), parameter :: min_abs_velocity_step = 2.0e-05_dp
        real(dp), parameter :: min_abs_velocity_step_capture = 1.0e-03_dp
        real(dp), parameter :: emergency_min_fraction = 1.0e-01_dp
        real(dp), parameter :: velocity_tolerance_on = 2.0e-05_dp
        real(dp), parameter :: velocity_tolerance_off = 5.0e-06_dp
        real(dp), parameter :: filter_alpha = 1.0e-01_dp
        real(dp), parameter :: secant_relaxation = 2.0e-01_dp
        real(dp), parameter :: feedback_sign_default = -1.0_dp
        real(dp), parameter :: heat_release_cut_fraction = 1.0e-08_dp
        real(dp), parameter :: heat_release_valid_relative = 1.0e-07_dp
        real(dp), parameter :: heat_release_valid_absolute = 1.0e-20_dp
        real(dp), parameter :: position_relaxation_time = 5.0e-02_dp
        real(dp), parameter :: position_tolerance_cells = 4.0_dp
        real(dp), parameter :: capture_position_tolerance_cells = 12.0_dp
        real(dp), parameter :: capture_velocity_threshold = 2.0e-02_dp
        real(dp), parameter :: outlet_guard_cells = 20.0_dp
        real(dp), parameter :: outlet_guard_fraction = 1.5e-01_dp
        real(dp), parameter :: min_secant_du = 5.0e-05_dp
        real(dp), parameter :: max_response_slope_abs = 1.0e+03_dp
        real(dp), parameter :: persistent_error_factor = 5.0_dp
        real(dp), parameter :: ramp_settle_fraction = 1.0e-03_dp
        real(dp), parameter :: tiny_weight = tiny(1.0_dp)
        logical, parameter :: use_secant_control = .false.

        integer, parameter :: workflow_flamelet_sl = 1
        integer, parameter :: workflow_anchor_observation = 2
        integer :: active_workflow
        logical :: enable_flamelet_output
        logical :: enable_drift_measurement
        logical :: pause_after_sl_measurement
        integer, parameter :: stage_anchor_control = 0
        integer, parameter :: stage_flamelet_ready = 1
        integer, parameter :: stage_measurement_ramp = 2
        integer, parameter :: stage_drift_measurement = 3
        integer, parameter :: stage_measurement_done = 4
        integer, parameter :: stage_measurement_failed = 5
        integer, parameter :: min_hist_for_sl_measurement = 40
        real(dp), parameter :: measurement_delta_fraction = 3.0e-02_dp
        real(dp), parameter :: measurement_delta_sign = 1.0_dp
        real(dp), parameter :: measurement_delta_min = 5.0e-03_dp
        real(dp), parameter :: measurement_delta_max_fraction = 3.0e-01_dp
        real(dp), parameter :: measurement_delta_growth = 2.0_dp
        real(dp) :: measurement_max_duration
        real(dp), parameter :: measurement_no_motion_displacement_cells = 0.5_dp
        integer, parameter :: measurement_max_attempts = 5
        real(dp), parameter :: measurement_ramp_time = 1.0e-03_dp
        real(dp), parameter :: measurement_settle_time = 1.0e-03_dp
        real(dp) :: measurement_min_displacement_cells
        real(dp), parameter :: measurement_min_duration = 5.0e-03_dp
        real(dp), parameter :: measurement_r2_min = 9.95e-01_dp
        real(dp), parameter :: measurement_split_slope_rel_tol = 2.5e-01_dp
        real(dp), parameter :: measurement_split_slope_abs_tol = 2.0e-04_dp
        real(dp), parameter :: measurement_residual_cells = 1.0_dp

        integer :: dimensions, species_number, boundary_types
        integer :: front_axis
        integer :: i, j, k, dim, bound_number, specie_number, specie_index
        integer :: H2_index, H_index
        integer :: cons_inner_loop(3,2)
        integer :: active_track_number

        real(dp) :: cell_size(3), cell_volume
        real(dp) :: current_flame_location(3)
        real(dp) :: heat_release_centroid(3)
        real(dp) :: H_centroid(3)
        real(dp) :: Tgrad_centroid(3)
        real(dp) :: current_front_coord
        real(dp) :: flame_velocity_lsq, flame_velocity_filtered
        real(dp) :: diag_flame_velocity_lsq, diag_flame_velocity_filtered
        real(dp) :: control_velocity, position_error, position_velocity, position_control_error
        real(dp) :: front_spread, heat_release_integral
        real(dp) :: heat_release_max, heat_release_valid_limit, H_max, Tgrad_max
        real(dp) :: measured_inlet_velocity, proposed_inlet_velocity
        real(dp) :: du_raw, du_limited, max_velocity_step, X_H2
        real(dp) :: time_delay, time_track, time_control, response_settle_time, inlet_ramp_time
        real(dp) :: time_control_effective, response_settle_effective
        real(dp) :: coord(3), weight, qdot_cut, qdot, hval, tgrad
        real(dp) :: sum_weight, sum_s, sum_s2
        real(dp) :: sum_H_weight, sum_Tgrad_weight
        real(dp) :: s_coord
        real(dp) :: ramp_elapsed, ramp_fraction, ramp_residual, ramp_tolerance
        real(dp) :: target_step_for_log
        real(dp) :: sl_displacement, linear_r2, linear_rms, split_slope_diff
        real(dp) :: measurement_elapsed, measurement_displacement, measurement_delta
        logical :: measurement_linear_ok
        real(dp) :: front_position_tolerance, capture_position_tolerance
        real(dp) :: domain_front_min, domain_front_max, domain_front_length
        real(dp) :: outlet_guard_distance, outlet_distance
        real(dp) :: front_safe_min, front_safe_max
        real(dp), allocatable :: farfield_concentrations(:), concs(:)
        character(len=10), allocatable :: farfield_species_names(:)
        character(len=5) :: axis_names(3)
        character(len=20) :: flame_data_file
        character(len=500) :: av_header
        character(len=100) :: chemical_mechanism
        logical :: found_inlet_farfield, trace_success, flame_detected, control_performed
        logical :: measurement_enabled, ramp_settled, reset_history_after_log
        logical :: capture_mode, emergency_mode


        stabilized = .false.
        if (.not. this%enabled) return

        if (this%flame_stabilization%is_anchor()) then
            active_workflow = workflow_anchor_observation
        else if (this%flame_stabilization%is_laminar_burning_velocity()) then
            active_workflow = workflow_flamelet_sl
        else
            error stop 'Flame stabilization called with disabled/unknown problem-control mode'
        end if

        enable_flamelet_output = (active_workflow == workflow_flamelet_sl)
        enable_drift_measurement = (active_workflow == workflow_flamelet_sl)
        pause_after_sl_measurement = (active_workflow == workflow_flamelet_sl)

        dimensions = this%domain%get_domain_dimensions()
        axis_names = this%domain%get_axis_names()
        species_number = this%chem%chem_ptr%species_number
        boundary_types = this%boundary%bc_ptr%get_boundary_types()
        cons_inner_loop = this%domain%get_local_inner_cells_bounds()
        cell_size = this%mesh%mesh_ptr%get_cell_edges_length()

        front_axis = 1
        if (front_axis > dimensions) front_axis = 1

        cell_volume = 1.0_dp
        do dim = 1, dimensions
            cell_volume = cell_volume * cell_size(dim)
        end do
        front_position_tolerance = max(position_tolerance_cells * cell_size(front_axis), 1.0e-8_dp)
        capture_position_tolerance = max(capture_position_tolerance_cells * cell_size(front_axis), &
            front_position_tolerance)
        domain_front_min = (real(cons_inner_loop(front_axis,1),dp) - 0.5_dp) * cell_size(front_axis)
        domain_front_max = (real(cons_inner_loop(front_axis,2),dp) - 0.5_dp) * cell_size(front_axis)
        domain_front_length = max(domain_front_max - domain_front_min, cell_size(front_axis))
        outlet_guard_distance = max(outlet_guard_cells * cell_size(front_axis), &
            outlet_guard_fraction * domain_front_length)

        H2_index = this%chem%chem_ptr%get_chemical_specie_index('H2')
        H_index = this%chem%chem_ptr%get_chemical_specie_index('H')

        if (.not. allocated(this%state%time_hist)) then
            allocate(this%state%time_hist(max_hist_size), this%state%front_coord_hist(max_hist_size))
            this%state%time_hist = 0.0_dp
            this%state%front_coord_hist = 0.0_dp
        end if

        if (.not. allocated(this%state%diag_time_hist)) then
            allocate(this%state%diag_time_hist(max_diag_hist_size), this%state%diag_front_coord_hist(max_diag_hist_size))
            this%state%diag_time_hist = 0.0_dp
            this%state%diag_front_coord_hist = 0.0_dp
        end if

        if (.not. this%state%initialized) then
            allocate(concs(species_number))
            concs = 0.0_dp

            found_inlet_farfield = .false.
            do bound_number = 1, boundary_types
                if (this%boundary%bc_ptr%boundary_types(bound_number)%get_type_name() == 'inlet') then
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_concentrations(farfield_concentrations)
                    call this%boundary%bc_ptr%boundary_types(bound_number)%get_farfield_species_names(farfield_species_names)
                    found_inlet_farfield = .true.
                    exit
                end if
            end do

            if (found_inlet_farfield .and. allocated(farfield_species_names)) then
                do specie_number = 1, size(farfield_species_names)
                    specie_index = this%chem%chem_ptr%get_chemical_specie_index(farfield_species_names(specie_number))
                    if (specie_index >= 1 .and. specie_index <= species_number) then
                        concs(specie_index) = farfield_concentrations(specie_number)
                    end if
                end do
            end if

            if (sum(concs) > 0.0_dp .and. H2_index >= 1 .and. H2_index <= species_number) then
                X_H2 = concs(H2_index) / sum(concs) * 100.0_dp
            else
                X_H2 = 0.0_dp
            end if

            chemical_mechanism = trim(this%chem%chem_ptr%get_chemical_mechanism())
            this%state%data_table_filename = 'H2-Air_flamelet_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'
            this%state%chem_table_filename = 'H2-Air_chem_table_' // trim(chemical_mechanism) // '_' // &
                trim(str_r(X_H2)) // '_pcnt_' // trim(str_e(cell_size(1))) // '_dx.dat'

            this%state%previous_correction_time = time
            this%state%adaptive_gain = controller_gain_initial
            this%state%initialized = .true.
        end if

        if (.not. this%state%inlet_velocity_initialized) then
            this%state%inlet_velocity_target = this%inlet_velocity
            this%state%inlet_velocity_applied = this%inlet_velocity
            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
            this%state%ramp_start_time = time
            this%state%inlet_velocity_initialized = .true.
        end if

        time_delay = this%flame_stabilization%get_time_delay()
        time_track = this%flame_stabilization%get_time_track()
        time_control = this%flame_stabilization%get_time_control()
        response_settle_time = this%flame_stabilization%get_response_settle_time()
        measurement_max_duration = this%flame_stabilization%get_measurement_max_duration()
        measurement_min_displacement_cells = &
            this%flame_stabilization%get_measurement_min_displacement_cells()
        inlet_ramp_time = this%state%active_inlet_ramp_time

        if (inlet_ramp_time > 0.0_dp) then
            ramp_elapsed = max(time - this%state%ramp_start_time, 0.0_dp)
            ramp_fraction = min(ramp_elapsed / inlet_ramp_time, 1.0_dp)
        else
            ramp_fraction = 1.0_dp
        end if
        this%state%inlet_velocity_applied = this%state%ramp_start_velocity + &
            ramp_fraction * (this%state%inlet_velocity_target - this%state%ramp_start_velocity)
        call this%set_inlet_velocity(this%state%inlet_velocity_applied)

        ramp_residual = abs(this%state%inlet_velocity_target - this%state%inlet_velocity_applied)
        ramp_tolerance = max(ramp_settle_fraction * max(abs(this%state%inlet_velocity_target), min_abs_velocity_step), &
            1.0e-12_dp)
        ramp_settled = (ramp_residual <= ramp_tolerance .or. ramp_fraction >= 1.0_dp)

        if ((time - time_delay) / time_track <= real(this%state%track_counter + 1, dp)) return

        associate (T => this%T%s_ptr, &
                   Y => this%Y%v_ptr, &
                   E_f_prod_chem => this%E_f_prod_chem%s_ptr, &
                   bc => this%boundary%bc_ptr)

            if (.not. this%state%output_initialized) call initialize_output_file()

            heat_release_max = 0.0_dp
            H_max = 0.0_dp
            Tgrad_max = 0.0_dp
            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
                heat_release_max = max(heat_release_max, qdot)

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = abs(Y%pr(H_index)%cells(i,j,k))
                    H_max = max(H_max, hval)
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                Tgrad_max = max(Tgrad_max, tgrad)
            end do
            end do
            end do

            this%state%heat_release_peak_save = max(this%state%heat_release_peak_save, heat_release_max)
            heat_release_valid_limit = max(heat_release_valid_absolute, &
                heat_release_valid_relative * this%state%heat_release_peak_save)

            heat_release_centroid = 0.0_dp
            H_centroid = 0.0_dp
            Tgrad_centroid = 0.0_dp
            sum_weight = 0.0_dp
            sum_H_weight = 0.0_dp
            sum_Tgrad_weight = 0.0_dp
            sum_s = 0.0_dp
            sum_s2 = 0.0_dp
            heat_release_integral = 0.0_dp

            qdot_cut = heat_release_cut_fraction * heat_release_max

            do k = cons_inner_loop(3,1), cons_inner_loop(3,2)
            do j = cons_inner_loop(2,1), cons_inner_loop(2,2)
            do i = cons_inner_loop(1,1), cons_inner_loop(1,2)
                if (bc%bc_markers(i,j,k) /= 0) cycle

                coord = cell_center_coordinates(i,j,k)
                s_coord = coord(front_axis)

                qdot = max(E_f_prod_chem%cells(i,j,k), 0.0_dp)
                if (qdot > qdot_cut) then
                    weight = qdot * cell_volume
                    heat_release_centroid = heat_release_centroid + weight * coord
                    sum_weight = sum_weight + weight
                    sum_s = sum_s + weight * s_coord
                    sum_s2 = sum_s2 + weight * s_coord * s_coord
                    heat_release_integral = heat_release_integral + weight
                end if

                if (H_index >= 1 .and. H_index <= species_number) then
                    hval = max(Y%pr(H_index)%cells(i,j,k), 0.0_dp)
                    if (hval > 0.0_dp) then
                        weight = hval * cell_volume
                        H_centroid = H_centroid + weight * coord
                        sum_H_weight = sum_H_weight + weight
                    end if
                end if

                tgrad = temperature_gradient_norm(i,j,k)
                if (tgrad > 0.0_dp) then
                    weight = tgrad * cell_volume
                    Tgrad_centroid = Tgrad_centroid + weight * coord
                    sum_Tgrad_weight = sum_Tgrad_weight + weight
                end if
            end do
            end do
            end do

            trace_success = .false.
            flame_detected = .false.
            current_flame_location = 0.0_dp
            front_spread = 0.0_dp

            if (sum_weight > tiny_weight) then
                current_flame_location = heat_release_centroid / sum_weight
                front_spread = max(sum_s2 / sum_weight - (sum_s / sum_weight)**2, 0.0_dp)
                front_spread = sqrt(front_spread)
                trace_success = .true.
                flame_detected = (heat_release_max >= heat_release_valid_limit)
            else if (sum_H_weight > tiny_weight) then
                current_flame_location = H_centroid / sum_H_weight
                trace_success = .true.
            else if (sum_Tgrad_weight > tiny_weight) then
                current_flame_location = Tgrad_centroid / sum_Tgrad_weight
                trace_success = .true.
            end if

            if (.not. trace_success) then
                this%state%track_counter = this%state%track_counter + 1
                return
            end if

            current_front_coord = current_flame_location(front_axis)
            if (flame_detected .and. .not. this%state%front_reference_initialized) then
                this%state%front_reference_coord = current_front_coord
                this%state%front_reference_initialized = .true.
            end if

            front_safe_min = domain_front_min + outlet_guard_distance
            front_safe_max = domain_front_max - outlet_guard_distance
            if (front_safe_min >= front_safe_max) then
                front_safe_min = domain_front_min
                front_safe_max = domain_front_max
            end if

            diag_flame_velocity_lsq = 0.0_dp
            diag_flame_velocity_filtered = 0.0_dp
            if (flame_detected) then
                call append_diagnostic_history(time, current_front_coord)
                diag_flame_velocity_lsq = diagnostic_least_squares_velocity()
                if (this%state%diag_hist_count <= 2) then
                    diag_flame_velocity_filtered = diag_flame_velocity_lsq
                    this%state%diag_filtered_velocity_save = diag_flame_velocity_filtered
                else
                    diag_flame_velocity_filtered = (1.0_dp - filter_alpha) * this%state%diag_filtered_velocity_save + &
                        filter_alpha * diag_flame_velocity_lsq
                    this%state%diag_filtered_velocity_save = diag_flame_velocity_filtered
                end if
            else
                this%state%diag_hist_count = 0
                this%state%diag_filtered_velocity_save = 0.0_dp
            end if

            control_performed = .false.
            reset_history_after_log = .false.
            du_raw = 0.0_dp
            du_limited = 0.0_dp
            target_step_for_log = 0.0_dp
            flame_velocity_lsq = diag_flame_velocity_lsq
            flame_velocity_filtered = diag_flame_velocity_filtered
            position_error = 0.0_dp
            position_control_error = 0.0_dp
            position_velocity = 0.0_dp
            control_velocity = 0.0_dp
            sl_displacement = this%state%sl_displacement_save
            linear_r2 = this%state%measurement_r2_save
            linear_rms = this%state%measurement_rms_save
            split_slope_diff = this%state%measurement_split_slope_diff_save
            measurement_elapsed = 0.0_dp
            measurement_displacement = 0.0_dp
            measurement_linear_ok = .false.
            capture_mode = .false.
            emergency_mode = .false.

            if (this%state%front_reference_initialized) then
                position_error = current_front_coord - this%state%front_reference_coord
                position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = (position_control_error /= 0.0_dp) .or. &
                    (outlet_distance < outlet_guard_distance)
                emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                    (abs(position_control_error) > capture_position_tolerance)
            end if

            if (capture_mode) then
                response_settle_effective = response_settle_capture
                time_control_effective = time_control_capture
            else
                response_settle_effective = response_settle_time
                time_control_effective = time_control
            end if

            measurement_enabled = flame_detected .and. this%state%front_reference_initialized .and. ramp_settled .and. &
                ((this%state%control_stage == stage_anchor_control) .or. &
                 (this%state%control_stage == stage_flamelet_ready)) .and. &
                ((time - this%state%previous_correction_time) >= response_settle_effective)

            if (measurement_enabled) then
                call append_front_history(time, current_front_coord)
                flame_velocity_lsq = least_squares_velocity()

                if (this%state%hist_count <= 2) then
                    flame_velocity_filtered = flame_velocity_lsq
                    this%state%filtered_velocity_save = flame_velocity_filtered
                else
                    flame_velocity_filtered = (1.0_dp - filter_alpha) * this%state%filtered_velocity_save + &
                        filter_alpha * flame_velocity_lsq
                    this%state%filtered_velocity_save = flame_velocity_filtered
                end if

                position_error = current_front_coord - this%state%front_reference_coord
                position_control_error = safe_window_error(current_front_coord, front_safe_min, front_safe_max)
                position_velocity = position_control_error / max(position_relaxation_time, time_track)
                control_velocity = flame_velocity_filtered + position_velocity
                outlet_distance = domain_front_max - current_front_coord
                capture_mode = (position_control_error /= 0.0_dp) .or. &
                    (abs(flame_velocity_filtered) > capture_velocity_threshold) .or. &
                    (outlet_distance < outlet_guard_distance)
                emergency_mode = (outlet_distance < 0.5_dp * outlet_guard_distance) .or. &
                    (abs(position_control_error) > capture_position_tolerance)
                if (capture_mode) then
                    time_control_effective = time_control_capture
                else
                    time_control_effective = time_control
                end if
            else if (this%state%control_stage == stage_anchor_control) then
                call clear_control_history()
                if (.not. flame_detected) this%state%has_bracket = .false.
            else if (.not. flame_detected) then
                call clear_control_history()
                this%state%has_bracket = .false.
            end if

            if (this%state%control_stage == stage_measurement_ramp) then
                if (flame_detected .and. ramp_settled .and. &
                    (time - this%state%ramp_start_time) >= measurement_settle_time) then
                    call clear_control_history()
                    this%state%measurement_start_time = time
                    this%state%measurement_start_coord = current_front_coord
                    this%state%control_stage = stage_drift_measurement
                end if
            else if (this%state%control_stage == stage_drift_measurement) then
                if (flame_detected .and. ramp_settled) then
                    call append_front_history(time, current_front_coord)
                    flame_velocity_lsq = least_squares_velocity()
                    flame_velocity_filtered = flame_velocity_lsq
                    this%state%measurement_velocity_save = flame_velocity_lsq
                    control_velocity = flame_velocity_lsq
                    measurement_elapsed = time - this%state%measurement_start_time
                    measurement_displacement = current_front_coord - this%state%measurement_start_coord
                    call drift_linearity_diagnostics(linear_r2, linear_rms, split_slope_diff)
                    this%state%measurement_r2_save = linear_r2
                    this%state%measurement_rms_save = linear_rms
                    this%state%measurement_split_slope_diff_save = split_slope_diff
                    sl_displacement = this%state%inlet_velocity_target - flame_velocity_lsq
                    this%state%sl_displacement_save = sl_displacement
                    measurement_linear_ok = (this%state%hist_count >= min_hist_for_sl_measurement) .and. &
                        (measurement_elapsed >= measurement_min_duration) .and. &
                        (abs(measurement_displacement) >= measurement_min_displacement_cells * cell_size(front_axis)) .and. &
                        (linear_r2 >= measurement_r2_min) .and. &
                        (linear_rms <= measurement_residual_cells * cell_size(front_axis)) .and. &
                        (split_slope_diff <= max(measurement_split_slope_abs_tol, &
                            measurement_split_slope_rel_tol * max(abs(flame_velocity_lsq), velocity_tolerance_on)))
                    if (measurement_linear_ok) then
                        this%state%measurement_velocity_save = flame_velocity_lsq
                        this%state%control_stage = stage_measurement_done
                        this%state%stabilization_counter = stable_required_count
                        call write_laminar_velocity_once()
                    else if (measurement_elapsed >= measurement_max_duration .and. &
                        abs(measurement_displacement) < measurement_no_motion_displacement_cells * cell_size(front_axis)) then
                        if (this%state%measurement_attempt < measurement_max_attempts) then
                            this%state%measurement_attempt = this%state%measurement_attempt + 1
                            this%state%current_measurement_delta = min(measurement_delta_growth * max(this%state%current_measurement_delta, &
                                measurement_delta_min), measurement_delta_max_fraction * &
                                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step))
                            this%state%measurement_inlet_velocity = max(this%state%stabilized_inlet_velocity + &
                                measurement_delta_sign * this%state%current_measurement_delta, 0.0_dp)
                            this%state%inlet_velocity_target = this%state%measurement_inlet_velocity
                            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
                            this%state%ramp_start_time = time
                            this%state%active_inlet_ramp_time = measurement_ramp_time
                            this%state%control_stage = stage_measurement_ramp
                            this%state%sl_displacement_save = 0.0_dp
                            this%state%measurement_velocity_save = 0.0_dp
                            this%state%measurement_r2_save = 0.0_dp
                            this%state%measurement_rms_save = 0.0_dp
                            this%state%measurement_split_slope_diff_save = 0.0_dp
                            call clear_control_history()
                        else
                            this%state%control_stage = stage_measurement_failed
                        end if
                    end if
                else
                    call clear_control_history()
                end if
            end if

            if (measurement_enabled .and. &
                this%state%hist_count >= merge(min_hist_for_capture, min_hist_for_control, capture_mode)) then
                if ((time - this%state%previous_correction_time) >= time_control_effective) then
                    if (control_action_needed()) then
                        measured_inlet_velocity = this%state%inlet_velocity_target

                        call update_adaptive_gain(control_velocity)
                        call update_bracket(measured_inlet_velocity, control_velocity)
                        call choose_new_inlet_target(measured_inlet_velocity, control_velocity, &
                            proposed_inlet_velocity, du_raw, du_limited)

                        if (abs(proposed_inlet_velocity - measured_inlet_velocity) > 0.0_dp) then
                            this%state%inlet_velocity_target = proposed_inlet_velocity
                            target_step_for_log = this%state%inlet_velocity_target - measured_inlet_velocity
                            if (capture_mode) then
                                this%state%active_inlet_ramp_time = inlet_ramp_time_capture
                            else
                                this%state%active_inlet_ramp_time = inlet_ramp_time_fine
                            end if
                            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
                            this%state%ramp_start_time = time
                            this%state%previous_correction_time = time
                            this%state%correction_counter = this%state%correction_counter + 1
                            control_performed = .true.
                            reset_history_after_log = .true.
                        end if

                        this%state%previous_control_inlet_velocity = measured_inlet_velocity
                        this%state%previous_control_velocity = control_velocity
                        this%state%have_previous_control_point = .true.
                    end if
                end if
            end if

            if (this%state%control_stage == stage_anchor_control) then
                if (measurement_enabled .and. ramp_settled .and. this%state%hist_count >= min_hist_for_control .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    this%state%stabilization_counter = this%state%stabilization_counter + 1
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    this%state%stabilization_counter = 0
                end if
            else if (this%state%control_stage == stage_flamelet_ready) then
                if (measurement_enabled .and. ramp_settled .and. this%state%hist_count >= min_hist_for_control .and. &
                    abs(diag_flame_velocity_filtered) < velocity_tolerance_off .and. &
                    position_control_error == 0.0_dp) then
                    this%state%post_flamelet_hold_counter = this%state%post_flamelet_hold_counter + 1
                    this%state%stabilization_counter = this%state%post_flamelet_hold_counter
                else if (control_performed .or. .not. flame_detected .or. &
                    abs(diag_flame_velocity_filtered) > velocity_tolerance_on .or. &
                    position_control_error /= 0.0_dp) then
                    this%state%post_flamelet_hold_counter = 0
                    this%state%stabilization_counter = 0
                end if
            else if (this%state%control_stage == stage_measurement_done) then
                this%state%stabilization_counter = stable_required_count
            end if

            if (this%state%control_stage == stage_anchor_control .or. this%state%control_stage == stage_flamelet_ready) then
                flame_velocity_lsq = diag_flame_velocity_lsq
                flame_velocity_filtered = diag_flame_velocity_filtered
            else if (this%state%control_stage == stage_measurement_done) then
                flame_velocity_lsq = this%state%measurement_velocity_save
                flame_velocity_filtered = this%state%measurement_velocity_save
                control_velocity = this%state%measurement_velocity_save
                sl_displacement = this%state%sl_displacement_save
                linear_r2 = this%state%measurement_r2_save
                linear_rms = this%state%measurement_rms_save
                split_slope_diff = this%state%measurement_split_slope_diff_save
            end if

            active_track_number = this%state%track_counter
            call write_tracking_line()

            if (reset_history_after_log) call clear_control_history()

            this%state%track_counter = this%state%track_counter + 1

            if (this%state%control_stage == stage_anchor_control .and. &
                this%state%stabilization_counter >= flamelet_required_count) then
                if (active_workflow == workflow_flamelet_sl) then
                    if (enable_flamelet_output) call request_flamelet_tables_once()
                    this%state%control_stage = stage_flamelet_ready
                    this%state%stabilization_counter = 0
                    this%state%post_flamelet_hold_counter = 0
                    call clear_control_history()
                else
                    stabilized = .false.
                end if

            else if (this%state%control_stage == stage_flamelet_ready .and. &
                this%state%post_flamelet_hold_counter >= post_flamelet_hold_required_count) then
                if (enable_drift_measurement .and. dimensions == 1) then
                    call start_drift_measurement_ramp(time)
                else if (pause_after_sl_measurement) then
                    stabilized = .true.
                end if

            else if (this%state%control_stage == stage_measurement_done) then
                if (pause_after_sl_measurement) stabilized = .true.
            end if

        end associate

    contains

        subroutine initialize_output_file()
            integer :: local_dim

            write(flame_data_file,'(A,I0,A)') 'av_flame_data_', this%load_counter, '.dat'
            open(newunit = this%state%flame_loc_unit, file = flame_data_file, status = 'replace', form = 'formatted')

            av_header = 'VARIABLES="time" '
            do local_dim = 1, dimensions
                av_header = trim(av_header) // '"xf_' // trim(axis_names(local_dim)) // '" '
            end do
            av_header = trim(av_header) // &
                '"Vfl_lsq" "Vfl_filtered" "Vfl_diag_lsq" "Vfl_diag_filtered" "Vfl_measurement_lsq" "Vcontrol" "pos_error" "x_ref" '
            av_header = trim(av_header) // '"U_in_applied" "U_in_target" '
            av_header = trim(av_header) // '"dU_target" "adaptive_gain" "hist_count" ' // &
                '"measurement_on" "bracket_on" "capture_on" "emergency_on" "flame_detected" '
            av_header = trim(av_header) // '"front_spread" "Qint" "Qmax" "Qvalid" "Hmax" "Tgradmax" '
            av_header = trim(av_header) // '"stage" "SL_disp" "lin_R2" "lin_RMS" "split_dV" '
            av_header = trim(av_header) // '"corr_count" "stab_count" "track_count"'

            write(this%state%flame_loc_unit,'(A)') trim(av_header)
            this%state%output_initialized = .true.
        end subroutine initialize_output_file

        function cell_center_coordinates(ii,jj,kk) result(xc)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: xc(3)

            xc = 0.0_dp
            xc(1) = (real(ii,dp) - 0.5_dp) * cell_size(1)
            if (dimensions >= 2) xc(2) = (real(jj,dp) - 0.5_dp) * cell_size(2)
            if (dimensions >= 3) xc(3) = (real(kk,dp) - 0.5_dp) * cell_size(3)
        end function cell_center_coordinates

        function temperature_gradient_norm(ii,jj,kk) result(grad_norm)
            integer, intent(in) :: ii, jj, kk
            real(dp) :: grad_norm
            real(dp) :: g2, gd

            g2 = 0.0_dp
            if (dimensions >= 1) then
                if (ii > cons_inner_loop(1,1) .and. ii < cons_inner_loop(1,2)) then
                    gd = (this%T%s_ptr%cells(ii+1,jj,kk) - this%T%s_ptr%cells(ii-1,jj,kk)) / (2.0_dp * cell_size(1))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 2) then
                if (jj > cons_inner_loop(2,1) .and. jj < cons_inner_loop(2,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj+1,kk) - this%T%s_ptr%cells(ii,jj-1,kk)) / (2.0_dp * cell_size(2))
                    g2 = g2 + gd * gd
                end if
            end if
            if (dimensions >= 3) then
                if (kk > cons_inner_loop(3,1) .and. kk < cons_inner_loop(3,2)) then
                    gd = (this%T%s_ptr%cells(ii,jj,kk+1) - this%T%s_ptr%cells(ii,jj,kk-1)) / (2.0_dp * cell_size(3))
                    g2 = g2 + gd * gd
                end if
            end if
            grad_norm = sqrt(g2)
        end function temperature_gradient_norm

        subroutine append_front_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (this%state%hist_count < max_hist_size) then
                this%state%hist_count = this%state%hist_count + 1
                this%state%time_hist(this%state%hist_count) = t_new
                this%state%front_coord_hist(this%state%hist_count) = s_new
            else
                this%state%time_hist(1:max_hist_size-1) = this%state%time_hist(2:max_hist_size)
                this%state%front_coord_hist(1:max_hist_size-1) = this%state%front_coord_hist(2:max_hist_size)
                this%state%time_hist(max_hist_size) = t_new
                this%state%front_coord_hist(max_hist_size) = s_new
            end if
        end subroutine append_front_history

        subroutine append_diagnostic_history(t_new, s_new)
            real(dp), intent(in) :: t_new, s_new

            if (this%state%diag_hist_count < max_diag_hist_size) then
                this%state%diag_hist_count = this%state%diag_hist_count + 1
                this%state%diag_time_hist(this%state%diag_hist_count) = t_new
                this%state%diag_front_coord_hist(this%state%diag_hist_count) = s_new
            else
                this%state%diag_time_hist(1:max_diag_hist_size-1) = this%state%diag_time_hist(2:max_diag_hist_size)
                this%state%diag_front_coord_hist(1:max_diag_hist_size-1) = this%state%diag_front_coord_hist(2:max_diag_hist_size)
                this%state%diag_time_hist(max_diag_hist_size) = t_new
                this%state%diag_front_coord_hist(max_diag_hist_size) = s_new
            end if
        end subroutine append_diagnostic_history

        subroutine clear_control_history()
            this%state%hist_count = 0
            this%state%filtered_velocity_save = 0.0_dp
            this%state%time_hist = 0.0_dp
            this%state%front_coord_hist = 0.0_dp
        end subroutine clear_control_history

        function least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (this%state%hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(this%state%time_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            s_av = sum(this%state%front_coord_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = 1, this%state%hist_count
                numerator = numerator + (this%state%time_hist(n) - t_av) * (this%state%front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity

        function diagnostic_least_squares_velocity() result(vfit)
            real(dp) :: vfit
            integer :: n
            real(dp) :: t_av, s_av, numerator, denominator

            if (this%state%diag_hist_count < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = sum(this%state%diag_time_hist(1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            s_av = sum(this%state%diag_front_coord_hist(1:this%state%diag_hist_count)) / real(this%state%diag_hist_count, dp)
            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = 1, this%state%diag_hist_count
                numerator = numerator + (this%state%diag_time_hist(n) - t_av) * &
                    (this%state%diag_front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%diag_time_hist(n) - t_av)**2
            end do

            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function diagnostic_least_squares_velocity

        subroutine drift_linearity_diagnostics(r2_out, rms_out, split_slope_diff_out)
            real(dp), intent(out) :: r2_out, rms_out, split_slope_diff_out
            integer :: n, mid
            real(dp) :: v_all, t_av, s_av, intercept, ss_tot, ss_res, residual
            real(dp) :: v_first, v_second

            r2_out = 0.0_dp
            rms_out = huge(1.0_dp)
            split_slope_diff_out = huge(1.0_dp)
            if (this%state%hist_count < 4) return

            v_all = least_squares_velocity()
            t_av = sum(this%state%time_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            s_av = sum(this%state%front_coord_hist(1:this%state%hist_count)) / real(this%state%hist_count, dp)
            intercept = s_av - v_all * t_av
            ss_tot = 0.0_dp
            ss_res = 0.0_dp
            do n = 1, this%state%hist_count
                residual = this%state%front_coord_hist(n) - (intercept + v_all * this%state%time_hist(n))
                ss_res = ss_res + residual * residual
                ss_tot = ss_tot + (this%state%front_coord_hist(n) - s_av)**2
            end do
            rms_out = sqrt(ss_res / real(this%state%hist_count, dp))
            if (ss_tot > tiny(ss_tot)) then
                r2_out = max(0.0_dp, 1.0_dp - ss_res / ss_tot)
            else
                r2_out = 0.0_dp
            end if

            mid = this%state%hist_count / 2
            v_first = least_squares_velocity_range(1, mid)
            v_second = least_squares_velocity_range(mid + 1, this%state%hist_count)
            split_slope_diff_out = abs(v_second - v_first)
        end subroutine drift_linearity_diagnostics

        function least_squares_velocity_range(n_first, n_last) result(vfit)
            integer, intent(in) :: n_first, n_last
            real(dp) :: vfit
            integer :: n, n_local
            real(dp) :: t_av, s_av, numerator, denominator

            n_local = n_last - n_first + 1
            if (n_local < 2) then
                vfit = 0.0_dp
                return
            end if

            t_av = 0.0_dp
            s_av = 0.0_dp
            do n = n_first, n_last
                t_av = t_av + this%state%time_hist(n)
                s_av = s_av + this%state%front_coord_hist(n)
            end do
            t_av = t_av / real(n_local, dp)
            s_av = s_av / real(n_local, dp)

            numerator = 0.0_dp
            denominator = 0.0_dp
            do n = n_first, n_last
                numerator = numerator + (this%state%time_hist(n) - t_av) * (this%state%front_coord_hist(n) - s_av)
                denominator = denominator + (this%state%time_hist(n) - t_av)**2
            end do
            if (denominator > tiny(denominator)) then
                vfit = numerator / denominator
            else
                vfit = 0.0_dp
            end if
        end function least_squares_velocity_range

        pure function safe_window_error(value, lower_bound, upper_bound) result(error_out)
            real(dp), intent(in) :: value, lower_bound, upper_bound
            real(dp) :: error_out

            if (value < lower_bound) then
                error_out = value - lower_bound
            else if (value > upper_bound) then
                error_out = value - upper_bound
            else
                error_out = 0.0_dp
            end if
        end function safe_window_error

        logical function control_action_needed()
            control_action_needed = (abs(flame_velocity_filtered) > velocity_tolerance_on) .or. &
                (abs(position_control_error) > 0.0_dp .and. abs(control_velocity) > velocity_tolerance_on)
        end function control_action_needed

        subroutine update_adaptive_gain(v_current)
            real(dp), intent(in) :: v_current

            if (this%state%have_previous_control_point) then
                if (v_current * this%state%previous_control_velocity < 0.0_dp) then
                    this%state%same_sign_error_counter = 0
                    this%state%adaptive_gain = max(0.7_dp * this%state%adaptive_gain, controller_gain_min)
                else
                    this%state%same_sign_error_counter = this%state%same_sign_error_counter + 1
                    if (this%state%same_sign_error_counter >= 2) then
                        this%state%adaptive_gain = min(1.10_dp * this%state%adaptive_gain, controller_gain_max)
                    else
                        this%state%adaptive_gain = min(1.03_dp * this%state%adaptive_gain, controller_gain_max)
                    end if
                end if
            else
                this%state%same_sign_error_counter = 0
            end if
        end subroutine update_adaptive_gain

        subroutine update_bracket(u_current, v_current)
            real(dp), intent(in) :: u_current, v_current

            if (abs(v_current) <= velocity_tolerance_off) return

            if (this%state%have_previous_control_point) then
                if (abs(u_current - this%state%previous_control_inlet_velocity) >= min_secant_du .and. &
                    v_current * this%state%previous_control_velocity < 0.0_dp) then
                    this%state%bracket_u_a = this%state%previous_control_inlet_velocity
                    this%state%bracket_v_a = this%state%previous_control_velocity
                    this%state%bracket_u_b = u_current
                    this%state%bracket_v_b = v_current
                    this%state%has_bracket = .true.
                end if
            end if

            if (this%state%has_bracket) then
                if (v_current * this%state%bracket_v_a > 0.0_dp) then
                    this%state%bracket_u_a = u_current
                    this%state%bracket_v_a = v_current
                else if (v_current * this%state%bracket_v_b > 0.0_dp) then
                    this%state%bracket_u_b = u_current
                    this%state%bracket_v_b = v_current
                end if
                if (this%state%bracket_v_a * this%state%bracket_v_b > 0.0_dp) this%state%has_bracket = .false.
                if (abs(this%state%bracket_u_a - this%state%bracket_u_b) < min_secant_du) this%state%has_bracket = .false.
            end if
        end subroutine update_bracket

        subroutine choose_new_inlet_target(u_current, v_current, u_new, du_unlimited, du_final)
            real(dp), intent(in) :: u_current, v_current
            real(dp), intent(out) :: u_new, du_unlimited, du_final
            real(dp) :: dU, dV, response_slope, secant_step, bracket_target
            real(dp) :: gain_effective, max_fraction_effective, min_step_effective

            if (emergency_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_emergency
                min_step_effective = max(min_abs_velocity_step_capture, &
                    emergency_min_fraction * max(abs(u_current), min_abs_velocity_step))
            else if (capture_mode) then
                gain_effective = controller_gain_capture
                max_fraction_effective = controller_max_fraction_capture
                min_step_effective = min_abs_velocity_step_capture
            else
                gain_effective = this%state%adaptive_gain
                max_fraction_effective = controller_max_fraction
                min_step_effective = min_abs_velocity_step
            end if

            max_velocity_step = max(max_fraction_effective * &
                max(abs(u_current), min_abs_velocity_step), min_step_effective)

            du_unlimited = feedback_sign_default * gain_effective * v_current
            if (emergency_mode) then
                if (abs(du_unlimited) < min_step_effective) then
                    du_unlimited = sign(min_step_effective, &
                        feedback_sign_default * max(abs(v_current), velocity_tolerance_on))
                end if
            else if (this%state%has_bracket .and. .not. capture_mode) then
                bracket_target = 0.5_dp * (this%state%bracket_u_a + this%state%bracket_u_b)
                du_unlimited = bracket_target - u_current
            else if (use_secant_control .and. this%state%have_previous_control_point .and. .not. capture_mode) then
                dU = u_current - this%state%previous_control_inlet_velocity
                dV = v_current - this%state%previous_control_velocity
                if (abs(dU) >= min_secant_du .and. abs(dV) > velocity_tolerance_off) then
                    response_slope = dV / dU
                    if (response_slope * feedback_sign_default < 0.0_dp .and. &
                        abs(response_slope) > 1.0e-12_dp .and. &
                        abs(response_slope) < max_response_slope_abs) then
                        secant_step = -secant_relaxation * v_current / response_slope
                        if (abs(secant_step) <= 5.0_dp * max_velocity_step) then
                            du_unlimited = secant_step
                        end if
                    end if
                end if
            end if

            du_final = min(max(du_unlimited, -max_velocity_step), max_velocity_step)

            if (abs(control_velocity) > persistent_error_factor * velocity_tolerance_on .or. &
                abs(position_control_error) > 0.0_dp) then
                if (abs(du_final) < min_step_effective) then
                    if (du_unlimited /= 0.0_dp) then
                        du_final = sign(min_step_effective, du_unlimited)
                    else
                        du_final = feedback_sign_default * sign(min_step_effective, v_current)
                    end if
                end if
            end if

            u_new = max(u_current + du_final, 0.0_dp)
        end subroutine choose_new_inlet_target

        subroutine request_flamelet_tables_once()
            if (.not. this%state%flamelet_output_written) then
                this%requested_chem_table_filename = trim(this%state%chem_table_filename)
                this%requested_data_table_filename = trim(this%state%data_table_filename)
                this%flamelet_output_pending = .true.
                this%state%flamelet_output_written = .true.
            end if
        end subroutine request_flamelet_tables_once

        subroutine write_laminar_velocity_once()
            integer :: sl_unit
            character(len=200) :: sl_file_name

            if (.not. this%state%sl_output_written) then
                write(sl_file_name,'(A,I0,A)') 'laminar_flame_velocity_', this%load_counter, '.dat'
                open(newunit = sl_unit, file = sl_file_name, status = 'replace', form = 'formatted')
                write(sl_unit,'(A)') &
                    'VARIABLES="time" "U_in" "Vfl_measurement" "SL_disp" "lin_R2" "lin_RMS" "split_dV" "attempt"'
                write(sl_unit,'(100E20.12)') time, this%state%measurement_inlet_velocity, this%state%measurement_velocity_save, &
                    this%state%sl_displacement_save, this%state%measurement_r2_save, this%state%measurement_rms_save, &
                    this%state%measurement_split_slope_diff_save, real(this%state%measurement_attempt,dp)
                close(sl_unit)
                this%state%sl_output_written = .true.
            end if
        end subroutine write_laminar_velocity_once

        subroutine start_drift_measurement_ramp(t_now)
            real(dp), intent(in) :: t_now

            this%state%stabilized_inlet_velocity = this%state%inlet_velocity_target
            this%state%measurement_attempt = 1
            this%state%current_measurement_delta = max(measurement_delta_fraction * &
                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step), measurement_delta_min)
            this%state%current_measurement_delta = min(this%state%current_measurement_delta, measurement_delta_max_fraction * &
                max(abs(this%state%stabilized_inlet_velocity), min_abs_velocity_step))
            this%state%measurement_inlet_velocity = max(this%state%stabilized_inlet_velocity + &
                measurement_delta_sign * this%state%current_measurement_delta, 0.0_dp)
            this%state%inlet_velocity_target = this%state%measurement_inlet_velocity
            this%state%ramp_start_velocity = this%state%inlet_velocity_applied
            this%state%ramp_start_time = t_now
            this%state%active_inlet_ramp_time = measurement_ramp_time
            this%state%sl_displacement_save = 0.0_dp
            this%state%measurement_velocity_save = 0.0_dp
            this%state%measurement_r2_save = 0.0_dp
            this%state%measurement_rms_save = 0.0_dp
            this%state%measurement_split_slope_diff_save = 0.0_dp
            this%state%control_stage = stage_measurement_ramp
            this%state%stabilization_counter = 0
            this%state%post_flamelet_hold_counter = 0
            call clear_control_history()
        end subroutine start_drift_measurement_ramp

        subroutine write_tracking_line()
            real(dp) :: measurement_flag, bracket_flag, capture_flag, emergency_flag, flame_detected_flag

            if (measurement_enabled) then
                measurement_flag = 1.0_dp
            else
                measurement_flag = 0.0_dp
            end if

            if (this%state%has_bracket) then
                bracket_flag = 1.0_dp
            else
                bracket_flag = 0.0_dp
            end if

            if (capture_mode) then
                capture_flag = 1.0_dp
            else
                capture_flag = 0.0_dp
            end if

            if (emergency_mode) then
                emergency_flag = 1.0_dp
            else
                emergency_flag = 0.0_dp
            end if

            if (flame_detected) then
                flame_detected_flag = 1.0_dp
            else
                flame_detected_flag = 0.0_dp
            end if

            write(this%state%flame_loc_unit,'(100E20.12)') &
                time, current_flame_location(1:dimensions), &
                flame_velocity_lsq, flame_velocity_filtered, diag_flame_velocity_lsq, &
                diag_flame_velocity_filtered, this%state%measurement_velocity_save, &
                control_velocity, position_error, this%state%front_reference_coord, &
                this%state%inlet_velocity_applied, this%state%inlet_velocity_target, &
                target_step_for_log, this%state%adaptive_gain, real(this%state%hist_count,dp), measurement_flag, bracket_flag, &
                capture_flag, emergency_flag, flame_detected_flag, &
                front_spread, heat_release_integral, heat_release_max, heat_release_valid_limit, H_max, Tgrad_max, &
                real(this%state%control_stage,dp), sl_displacement, linear_r2, linear_rms, split_slope_diff, &
                real(this%state%correction_counter,dp), real(this%state%stabilization_counter,dp), real(active_track_number,dp)
        end subroutine write_tracking_line

    end subroutine solve

end module flame_stabilization_solver_class
