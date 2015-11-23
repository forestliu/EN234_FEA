subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
!  use Mesh
  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
  use Printparameters, only : user_print_units                    ! Unit numbers
  use Printparameters, only : user_print_parameters               ! List of user supplied parameters
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  lmn
  integer ::  status
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
  real (prec)  :: vol_averaged_stress(6)
  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element
  real (prec) :: J_integral_value


!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!  element state variables (if they exist) in an element.   The element is specified by the user.
!  The first six state variables (which are usually the stresses) are printed along with the strains.
!
!
   write(*,*) 'user print is called!'
   allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)

   if (status/=0) then
      write(IOW,*) ' Error in subroutine user_print'
      write(IOW,*) ' Unable to allocate memory for state variables '
      stop
   endif

   lmn = int(user_print_parameters(1))     ! The element number

   call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
                                                       n_state_vars_per_intpt,vol_averaged_stress)
    if (TIME<1.d-12) then
      if (n_state_vars_per_intpt<6) then
        write(user_print_units(1),'(A)') 'VARIABLES = e11,s11'
      else
         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
      endif
    endif

   if (n_state_vars_per_intpt<6) then
      write(*,*) 'usr print'
      write(user_print_units(1),'(2(1x,D12.5))') vol_averaged_strain(1),vol_averaged_stress(1)
   else
      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
      write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
   endif
!
!
    !write(*,*) 'debug'
!    call compute_J_integral(J_integral_value)
!    write(user_print_units(1),'(A)') 'J_integral = '
!    write(user_print_units(1),'(D12.5)') J_integral_value
end subroutine user_print

subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_vars,length_output_array, &
                                                                       n_state_vars_per_intpt,vol_averaged_stress)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNbdx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_stress(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i,j,a
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  new_strain(6),dev_strain(6),vol_strain,stress(6)
    real (prec)  ::  E,xnu,Y,dedt0,m,q1,q2,q3,fn,en,sn,fc,ff              ! Material properties
    real (prec)  ::  eta,deta,JJ,F(3,3),dF(3,3),F_mid(3,3),dL(3,3),dLb(3,3)
    real (prec)  ::  dNdy(length_node_array,3),dNbdy(length_node_array,3),dR(3,3),dW(3,3),deps(3,3),I3(3,3)
    real (prec)  ::  dof_new(length_dof_array),inv_F(3,3),inv_F_mid(3,3),temp,temp3(3,3)
    real (prec)  ::  s0(6),eps_b_mat0,Vf0,s_new(3,3),kstress(6)
    real (prec)  ::  eps_b_mat,Vf


    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif

    E = element_properties(1)
    xnu = element_properties(2)
    Y = element_properties(3)
    dedt0 = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    sn = element_properties(10)
    en = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)
    eta = 0.d0
    deta = 0.d0
    dNbdy = 0.d0
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        !dof_new = 0.d0
        !dof_new = dof_total + dof_increment
        !F = 0.d0
        !F(1,1) = 1.d0
        !F(2,2) = 1.d0
        !F(3,3) = 1.d0
        dF = 0.d0
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i = 1,3 ! Compute deformation gradient tensor F
            do j = 1,3
                do a= 1,n_nodes
                    !F(i,j) = F(i,j) +  dof_new(3*(a-1)+i)*dNdx(a,j)
                    dF(i,j) = dF(i,j) +  dof_increment(3*(a-1)+i)*dNdx(a,j)
                    F_mid(i,j) = F_mid(i,j) + (dof_total(3*(a-1)+i) + &
                            0.5d0*dof_increment(3*(a-1)+i))*dNdx(a,j)
                end do
            end do
        end do

        dNdy = 0.d0
        !call invert_small(F,inv_F,JJ)
        call invert_small(F_mid,inv_F_mid,JJ)
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F_mid(1:3,1:3))
        dL(1:3,1:3) = matmul(dF(1:3,1:3),inv_F_mid(1:3,1:3))

        ! compute various volume integration
        eta = eta + w(kint) * determinant * JJ
        deta = deta + w(kint) * determinant * JJ * (dL(1,1) + dL(2,2) + dL(3,3))
        dNbdy(1:n_nodes,1:3) = dNbdy(1:n_nodes,1:3) + dNdy(1:n_nodes,1:3) * w(kint) * JJ * determinant
    end do
    deta = deta / eta
    dNbdy(1:n_nodes,1:3) = dNbdy(1:n_nodes,1:3) / eta

    I3 = 0.D0
    I3(1,1) = 1.d0
    I3(2,2) = 1.d0
    I3(3,3) = 1.d0
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        dF = 0.d0
        F_mid = 0.d0
        F_mid(1,1) = 1.d0
        F_mid(2,2) = 1.d0
        F_mid(3,3) = 1.d0
        do i = 1,3 ! Compute deformation gradient tensor F
            do j = 1,3
                do a= 1,n_nodes
                    !F(i,j) = F(i,j) +  dof_new(3*(a-1)+i)*dNdx(a,j)
                    dF(i,j) = dF(i,j) +  dof_increment(3*(a-1)+i)*dNdx(a,j)
                    F_mid(i,j) = F_mid(i,j) + (dof_total(3*(a-1)+i) + &
                            0.5d0*dof_increment(3*(a-1)+i))*dNdx(a,j)
                end do
            end do
        end do

        call invert_small(F_mid,inv_F_mid,JJ)
        dNdy = 0.d0
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F_mid(1:3,1:3))
        dLb(1:3,1:3) = matmul(dF(1:3,1:3),inv_F_mid(1:3,1:3))
        temp = dLb(1,1) + dLb(2,2) + dLb(3,3)
        dLb(1,1) = dLb(1,1) + (deta-temp)/3.d0
        dLb(2,2) = dLb(2,2) + (deta-temp)/3.d0
        dLb(3,3) = dLb(3,3) + (deta-temp)/3.d0

        deps = (dLb(1:3,1:3) + transpose(dLb(1:3,1:3)))/2.d0
        dW = (dLb(1:3,1:3) - transpose(dLb(1:3,1:3)))/2.d0
        call invert_small(I3(1:3,1:3)-dW(1:3,1:3)/2.d0,temp3(1:3,1:3),temp)
        dR = matmul(temp3(1:3,1:3), I3(1:3,1:3) + dW(1:3,1:3)/2.d0)

        s0(1:6) = initial_state_variables(8*(kint-1)+1:8*(kint-1)+6)
        eps_b_mat0 = initial_state_variables(8*(kint-1)+7)
        Vf0 = initial_state_variables(8*kint)
!        if (.not. isnan(s0(1))) then
!        write(*,*) 's0 =',s0
!        end if
        call get_stress_gurson(kstress(1:6),eps_b_mat,Vf,deps(1:3,1:3),s0(1:6),eps_b_mat0,Vf0, &
                                dR(1:3,1:3),E,xnu,Y, dedt0,m,q1,q2,q3,fn,en,sn,fc,ff)
        updated_state_variables(8*(kint-1)+1:8*(kint-1)+6) = kstress(1:6)
        updated_state_variables(8*(kint-1)+7) = eps_b_mat
        updated_state_variables(8*kint) = Vf
        stress(1:6) = kstress(1:6) / JJ

        vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + (stress(1:6))*w(kint)*determinant
        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif

        el_vol = el_vol + w(kint)*determinant

    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_stress = vol_averaged_stress/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)

    return




end subroutine compute_element_volume_average_3D

subroutine compute_J_integral(J_integral_value)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    real (prec), intent( out )  ::  J_integral_value


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node

    integer      :: status
    integer      :: iof
    integer      :: lmn               ! Element number
    integer      :: lmn_start,lmn_end ! First and last crack tip element
    integer      :: i                 ! Loop counter

!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:),dof_new(:)                            ! accumulated DOF, using usual element storage convention

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)
    real (prec) :: x_int(2), r0, r, dqdx1, dqdx2, du1dx2, du2dx2, strain_energy

    integer      :: n_points,kint

    real (prec)  ::  strain(3), dstrain(3), new_strain(3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    !real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    !real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties


    !
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(2,length_coord_array/2), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(3,length_dof_array), stat=status)

  !  Write your code to calculate the J integral here

  !  You will need to loop over the crack tip elements, and sum the contribution to the J integral from each element.
  !
  !  You can access the first and last crack tip element using
  !
    !write(*,*) 'believe'
      lmn_start = zone_list(2)%start_element
  !
     lmn_end = zone_list(2)%end_element

    J_integral_value = 0.d0
    do lmn = lmn_start , lmn_end

      !  The two subroutines below extract data for elements and nodes (see module Mesh.f90 for the source code for these subroutines)

        call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

        if (n_nodes == 3) n_points = 1
        if (n_nodes == 6) n_points = 9
        if (n_nodes == 4) n_points = 4
        !if (n_nodes == 10) n_points = 4
        if (n_nodes == 8) n_points = 9
        !if (n_nodes == 20) n_points = 27

    D = 0.d0
    E = element_properties(1)
    !write(*,*) 'E = ',E
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) !ENGINEERING STRAIN
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:2,1:2) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44
    !write(*,*) 'stiffness: ',D(1,1),',',D(1,2),',',D(2,1),',',D(2,2),',',D(3,3)
    !write(*,*) 'stupid!'
    !write(*,*) 'be happy'
    r0 = 6.D-4
    !write(*,*) 'god!'

        do i = 1, n_nodes
            iof = 2*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
            call extract_node_data(node_list(i),node_identifier,n_coords,x(1:2,i),n_dof, &
                                                 dof_increment(iof:iof+1),dof_total(iof:iof+1))
        end do
        !write(*,*) 'hey!'
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        dof_new = dof_total + dof_increment
        new_strain = strain + dstrain
        stress = matmul(D,new_strain)
        !write(*,*) 'real'
        strain_energy = 0.5*new_strain(1)*stress(1)+0.5*new_strain(2)*stress(2)+0.5*new_strain(3)*stress(3)

!        x_int = matmul(x(1:2,1:n_nodes),N(1:n_nodes))

        !write(*,*) 'code'
        du1dx2 = 0.D0
        du2dx2 = 0.D0
        x_int = 0.D0
        do i = 1 , n_nodes
            x_int(1) = x_int(1) + N(i)*x(1,i)
            x_int(2) = x_int(2) + N(i)*x(2,i)
            du1dx2 = du1dx2 + dNdx(i,2) * dof_new(2*i - 1)!matmul(dNdx(1:n_nodes,2),dof_total(1:2*n_nodes-1:2))
            du2dx2 = du2dx2 + dNdx(i,2) * dof_new(2*i)!matmul(dNdx(1:n_nodes,2),dof_total(2:2*n_nodes:2))
        end do
        r = dsqrt(x_int(1)**2+x_int(2)**2)
        dqdx1 = -x_int(1)/r/r0
        dqdx2 = -x_int(2)/r/r0
        J_integral_value = J_integral_value + w(kint) * determinant * (dqdx1*(stress(1)*du1dx2 + stress(3)*du2dx2) &
                        + dqdx2 * (stress(3) * du1dx2 + stress(2) * du2dx2) - strain_energy * dqdx2)
        !write(*,*) 'x1 = ',x_int(1), '; x2 = ', x_int(2)
        !write(*,*) 's11 = ', stress(1)
        !write(*,*) 'dqdx1 = ', dqdx1
        !write(*,*) 'strain energy = ', strain_energy
        !write(*,*) 'J-integral = ', J_integral_value

        !write(*,*) 'stark'

    end do
    !write(*,*) 'element--------------------------------------------------------'
    end do


    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)


    return




end subroutine compute_J_integral

