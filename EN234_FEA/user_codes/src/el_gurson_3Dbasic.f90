!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_gurson_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbdx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dB(length_dof_array,6)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12, Vel          ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    Vel = 0.d0
    dNbdx = 0.d0
    if ( element_identifier == 1) then ! B-bar element
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            dNbdx(1:n_nodes,1:3) = dNbdx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3) * w(kint) * determinant
            Vel = Vel + w(kint) * determinant
        end do
        dNbdx = dNbdx / Vel
    end if

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        if(element_identifier==1) then
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)+1/3*(dNbdx(1:n_nodes,1)-dNdx(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3)=1/3*(dNbdx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        B(1,3:3*n_nodes:3)=1/3*(dNbdx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3)=1/3*(dNbdx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)+1/3*(dNbdx(1:n_nodes,2) -dNdx(1:n_nodes,2))
        B(2,3:3*n_nodes:3)=1/3*(dNbdx(1:n_nodes,3) -dNdx(1:n_nodes,3))
        B(3,1:3*n_nodes-2:3)=1/3*(dNbdx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        B(3,2:3*n_nodes-1:3)=1/3*(dNbdx(1:n_nodes,2) -dNdx(1:n_nodes,2))
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)+1/3*(dNbdx(1:n_nodes,3) -dNdx(1:n_nodes,3))
        end if
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
  
    return
end subroutine el_gurson_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_gurson_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint,i,j,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    !real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E,xnu,Y,dedt0,m,q1,q2,q3,fn,en,sn,fc,ff              ! Material properties
    real (prec)  ::  eta,deta,JJ,F(3,3),dF(3,3),F_mid(3,3),dL(3,3),dLb(3,3),Vel
    real (prec)  ::  dNdy(n_nodes,3),dNbdy(n_nodes,3),dR(3,3),dW(3,3),deps(3,3),I3(3,3)
    real (prec)  ::  dof_new(length_dof_array),inv_F(3,3),inv_F_mid(3,3),temp,temp3(3,3)
    real (prec)  ::  s0(6),eps_b_mat0,Vf0,s_new(3,3),kstress(6)
    real (prec)  ::  eps_b_mat,Vf
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
!    D = 0.d0
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
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
!    write(*,*) 'dof_increment =',dof_increment
    eta = 0.d0
    deta = 0.d0
    dNbdy = 0.d0
    !     --  Loop over integration points
    element_deleted = .true.
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
        !if (Vf > ff) then
        !    element_deleted = .true.
        !    return
        !end if
        !write(*,*) 'kstress =', kstress
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = (dNbdy(1:n_nodes,1) - dNdy(1:n_nodes,1))/3.d0
        B(1,2:3*n_nodes-1:3) = (dNbdy(1:n_nodes,2) - dNdy(1:n_nodes,2))/3.d0
        B(1,3:3*n_nodes:3)   = (dNbdy(1:n_nodes,3) - dNdy(1:n_nodes,3))/3.d0
        B(2,1:3*n_nodes) = B(1,1:n_nodes)
        B(3,1:3*n_nodes) = B(1,1:n_nodes)
        B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = B(3,3:3*n_nodes:3) + dNdy(1:n_nodes,3)

        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        !write(*,*) 'B(1,1) =', B(1,1)
        !strain = matmul(B,dof_total)
        !dstrain = matmul(B,dof_increment)
!        kstress(1) = s_new(1,1)
!        kstress(2) = s_new(2,2)
!        kstress(3) = s_new(3,3)
!        kstress(4) = s_new(1,2)
!        kstress(5) = s_new(1,3)
!        kstress(6) = s_new(2,3)
!        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),kstress)*w(kint)*determinant
!        write(*,*) 'residual =', element_residual
!        write(*,*) 'kstress =', kstress
        !read(*,*)
        updated_state_variables(8*(kint-1)+1:8*(kint-1)+6) = kstress(1:6)
        updated_state_variables(8*(kint-1)+7) = eps_b_mat
        updated_state_variables(8*kint) = Vf
        if (Vf <= fn) element_deleted = .false.
    end do

    !write(*,*) 'element stiffness finished'

    return
end subroutine el_gurson_3dbasic_dynamic


subroutine get_stress_gurson(kstress,eps_b_mat,Vf,& ! Output
                            deps,s0,eps_b_mat0,Vf0,dR, & ! deformation of last step
                            E,xnu,Y, dedt0,m,q1,q2,q3,fn,en,sn,fc,ff) ! material properties
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : rotatesymvec
    implicit none

    real( prec ), intent(in)    :: deps(3,3)
    real( prec ), intent(in)    :: s0(6)
    real( prec ), intent(in)    :: eps_b_mat0
    real( prec ), intent(in)    :: Vf0
    real( prec ), intent(in)    :: dR(3,3)
    real( prec ), intent(in)    :: E,xnu,Y,dedt0,m,q1,q2,q3,fn,en,sn,fc,ff

    real( prec ), intent(out)    :: kstress(6)
    real( prec ), intent(out)    :: eps_b_mat
    real( prec ), intent(out)    :: Vf

    integer      ::  i,j
    real (prec)  ::  stress(3,3),sdev(3,3),dedev(3,3),p_star,S_star(3,3),devol,svol,phi
    real (prec)  ::  se_star,f_star,ff_b,dphidse,dphidp,dee,dev,test,deps_b_mat,se,p,temp
    real (prec)  ::  dF1ddee,dF1ddev,dF2ddee,dF2ddev,F1,F2,dseddee,dpddev,dphiddee,dphiddev
    real (prec)  ::  dtempddee,dtempddev
    real (prec)  ::  KK(2,2), inv_KK(2,2), RR(2), error, temp2,dde_vec(2), sdev_vec(6)

    !write(*,*) 'stress computing'
!    write(*,*) 'deps =', deps
    devol = (deps(1,1) + deps(2,2) + deps(3,3))/3.d0
    dedev(1:3,1:3) = deps(1:3,1:3)
    dedev(1,1) = dedev(1,1) - devol
    dedev(2,2) = dedev(2,2) - devol
    dedev(3,3) = dedev(3,3) - devol

    svol = (s0(1)+s0(2)+s0(3))/3.d0
    sdev(1,1) = s0(1) - svol
    sdev(2,2) = s0(2) - svol
    sdev(3,3) = s0(3) - svol
    sdev(1,2) = s0(4)
    sdev(2,1) = s0(4)
    sdev(1,3) = s0(5)
    sdev(3,1) = s0(5)
    sdev(2,3) = s0(6)
    sdev(3,2) = s0(6)
        if (.not. isnan(s0(1))) then
!        write(*,*) 's0 =',s0
        end if

    p_star = svol + E/(1-2.d0*xnu)*devol
!    write(*,*) 'p =',p_star
    sdev_vec(1:6) = s0(1:6)
    sdev_vec(1) = sdev_vec(1) - svol
    sdev_vec(2) = sdev_vec(2) - svol
    sdev_vec(3) = sdev_vec(3) - svol
    sdev_vec(1:6) = rotatesymvec(sdev_vec,dR)
    S_star(1:3,1:3) = E/(1.d0+xnu)*dedev(1:3,1:3) + &
!                            rotatesymvec(sdev_vec,dR)
                            matmul(dR(1:3,1:3),matmul(sdev(1:3,1:3),transpose(dR(1:3,1:3))))
    !write(*,*) 'S_star=',S_star
!    S_star(1,1) = S_star(1,1) + sdev_vec(1)
!    S_star(2,2) = S_star(2,2) + sdev_vec(2)
!    S_star(3,3) = S_star(3,3) + sdev_vec(3)
!    S_star(1,2) = S_star(1,2) + sdev_vec(4)
!    S_star(2,1) = S_star(2,1) + sdev_vec(4)
!    S_star(1,3) = S_star(1,3) + sdev_vec(5)
!    S_star(3,1) = S_star(3,1) + sdev_vec(5)
!    S_star(2,3) = S_star(2,3) + sdev_vec(6)
!    S_star(3,2) = S_star(3,2) + sdev_vec(6)
    se_star = 0.d0
    do i = 1,3
        do j = 1,3
            se_star = se_star + S_star(i,j)**2.d0
        end do
    end do
    se_star = dsqrt(1.5d0 * se_star)
    !write(*,*) 'se_star =',se_star
    ff_b = (q1+dsqrt(q1**2.d0 - q3))/q3
    if (Vf0 < fc) then
        f_star = Vf
    else if (Vf0 < ff) then
        f_star = fc + (ff_b - fc)/(ff - fc)* (Vf - fc)
    else
        f_star = ff_b
    end if

    test = se_star**2.d0/Y**2.d0 + 2.d0*q1*f_star*cosh(1.5d0*q2*p_star/Y) &
             - (1.d0+q3*f_star**2.d0)

    dseddee = -1.5d0 * E /(1.D0 + xnu)
    dpddev = - E/3.D0/(1-2.D0*xnu)
    !write(*,*) 'y =',y
    !write(*,*) 'dpddev =', dpddev
    if (test .le. 1.d-8) then ! elastic step
        deps_b_mat = 0.d0
        dev = 0.d0
        dee = 0.d0
    else
        error = 1.d0
        dee = 0.d0
        dev = 0.d0
        !write(*,*) 'iteration starts'
        do while (error >= 1.d-10)
            !write(*,*) error
            se = se_star - 1.5*E/(1.d0+xnu)*dee
            p = p_star - E/3.d0/(1.d0-2.d0*xnu)*dev
            !if (.not. isnan(se)) then
            !write(*,*) 'se =',se
            !end if
            !write(*,*) 'fstar =' , f_star
            phi = dsqrt(se**2.d0/Y**2.d0 + 2.d0*q1*f_star*cosh(1.5d0*q2*p/Y) &
                - (1.d0+q3*f_star**2.d0))
            dphidse = se / phi / Y**2.d0
            dphidp = 1.5d0 * q1 * q2 * f_star / phi / Y * sinh(1.5d0 * q2 * p / Y)
            dphiddee = dphidse * dseddee
            dphiddev = dphidp * dpddev
            temp = dsqrt(dphidse**2.d0 + 2.d0/9.d0*dphidp**2.d0)
            !write(*,*) 'phi =', phi
            !write(*,*) 'temp =', temp
            dtempddee = (2.d0*se*dseddee/phi/phi/Y**4.d0 - se**2.d0/Y**4.d0*2.d0/phi**3.d0*dphiddee &
                + 2.d0/9.d0*(dphidp**2.d0)*(-2.d0)/phi*dphiddee)/2.d0/temp
            dtempddev = (q1**2.d0*q2**2.d0*f_star**2/phi**2.d0/Y**2.d0*sinh(1.5*q2/Y*p)*cosh(1.5*q2/Y*p) &
                * 1.5*q2/Y*dpddev + 2.d0/9.d0*dphidp**2.d0*(-2.d0)/phi*dphiddev - se**2.d0/Y**4.d0*2.d0/ &
                phi**3.d0*dphiddev)/2.d0/temp
            F1 = temp * dee /dtime/dedt0 - dphidse * phi**m
            F2 = temp * dev /dtime/dedt0 - dphidp * phi**m

            dF1ddee = temp/dtime/dedt0 - dseddee / Y**2.d0 * phi**(m-1.d0) - &
                se/Y**2.d0*(m-1.d0)*phi**(m-2.d0)*dphiddee + dee/dtime/dedt0*dtempddee
            dF1ddev = - se/Y**2.d0*(m-1.d0)*phi**(m-2.d0)*dphiddev + dee/dtime/dedt0*dtempddev
            dF2ddee = dev/dtime/dedt0 * dtempddee - dphidp*(phi**m)*(m-1.d0)/(phi)*dphiddee
            dF2ddev = temp/dtime/dedt0 - dphidp*(phi**m)*(m-1.d0)/(phi)*dphiddev - 1.5d0*q1*q2*f_star/Y*phi**(m-1.d0)* &
                cosh(1.5*q2*p/Y) * dpddev*1.5d0*q2/Y + dev/dtime/dedt0 * dtempddev

            ! now start iteration with N-R method to solve for dee and dev
            ! assume each iteration takes the form: K*de = R
            RR = [-F1,-F2]
            KK(1,1) = dF1ddee
            KK(1,2) = dF1ddev
            KK(2,1) = dF2ddee
            KK(2,2) = dF2ddev
            call invert_small(KK,inv_KK,temp2)
            dde_vec = matmul(inv_KK,RR)
            error = dsqrt((dde_vec(1))**2.d0 + (dde_vec(2))**2.d0)
           ! write(*,*) 'K = ', KK
            !write(*,*) 'R = ', RR
            dee = dde_vec(1) + dee
            dev = dde_vec(2) + dev
!        if (isnan(dee)) then
!        write(*,*) 'dee =',dee
!        write(*,*) 'dev =',dev
!        stop
!        else
!
!        write(*,*) 'dee =',dee
!        write(*,*) 'dev =',dev
!        end if
        end do
        !write(*,*) 'iteration finished'
!        if (isnan(dee)) then
!        write(*,*) 'dee =',dee
!        write(*,*) 'dev =',dev
!        stop
!        else
!
!        write(*,*) 'dee =',dee
!        write(*,*) 'dev =',dev
!        end if
        se = se_star - 1.5*E/(1.d0+xnu)*dee
        p = p_star - E/3.d0/(1.d0-2.d0*xnu)*dev

        phi = dsqrt(se**2.d0/Y**2.d0 + 2.d0*q1*f_star*cosh(1.5d0*q2*p/Y) &
            - (1.d0+q3*f_star**2.d0))
        dphidse = se / phi / Y**2.d0
        dphidp = 1.5d0 * q1 * q2 * f_star / phi / Y * sinh(1.5d0 * q2 * p / Y)
        !dphiddee = dphidse * dseddee
        !dphiddev = dphidp * dpddev
        temp = dsqrt(dphidse**2.d0 + 2.d0/9.d0*dphidp**2.d0)

        deps_b_mat = dedt0*dtime/(1.d0-Vf0) *phi**m /temp *(dphidse*se + dphidp * p/3.d0)

    end if

    eps_b_mat = eps_b_mat0 + deps_b_mat
    Vf = 1.d0 + (Vf0 - 1.d0)*exp(-dev) + fn/sn/dsqrt(2.d0*3.14159265)*deps_b_mat* &
            exp(-0.5d0*((eps_b_mat0-en)/sn)**2.d0)
    !write(*,*) 'Vf =',Vf
    if (se_star /= 0.d0) then
    stress(1:3,1:3) = (1.d0 - dee * E/(1.d0+xnu)*1.5d0/se_star) * S_star(1:3,1:3)
    else
    stress(1:3,1:3) = S_star(1:3,1:3)
    end if
    kstress(1) = stress(1,1) + p_star - E/3.d0/(1.d0-2.d0*xnu)*dev
    kstress(2) = stress(2,2) + p_star - E/3.d0/(1.d0-2.d0*xnu)*dev
    kstress(3) = stress(3,3) + p_star - E/3.d0/(1.d0-2.d0*xnu)*dev
    kstress(4) = stress(1,2)
    kstress(5) = stress(1,3)
    kstress(6) = stress(2,3)
!    write(*,*) 'dev =',dev
!    write(*,*) 'dee =',dee
end subroutine get_stress_gurson


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_gurson_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNbdx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( inout )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,i,j,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dB(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E,xnu,Y,dedt0,m,q1,q2,q3,fn,en,sn,fc,ff              ! Material properties
    real (prec)  ::  eta,deta,JJ,F(3,3),dF(3,3),F_mid(3,3),dL(3,3),dLb(3,3),Vel
    real (prec)  ::  dNdy(n_nodes,3),dNbdy(n_nodes,3),dR(3,3),dW(3,3),deps(3,3),I3(3,3)
    real (prec)  ::  dof_new(length_dof_array),inv_F(3,3),inv_F_mid(3,3),temp,temp3(3,3)
    real (prec)  ::  s0(6),eps_b_mat0,Vf0,s_new(3,3),kstress(6)
    real (prec)  ::  eps_b_mat,Vf
    real (prec)  :: p, smises,s3                       ! Pressure and Mises stress

    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    !write(*,*) 'fieldvars called!'
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
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
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
!    write(*,*) 'dof_increment =',dof_increment
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
        !write(*,*) 'kstress =', kstress
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = (dNbdy(1:n_nodes,1) - dNdy(1:n_nodes,1))/3.d0
!        B(1,2:3*n_nodes-1:3) = (dNbdy(1:n_nodes,2) - dNdy(1:n_nodes,2))/3.d0
!        B(1,3:3*n_nodes:3)   = (dNbdy(1:n_nodes,3) - dNdy(1:n_nodes,3))/3.d0
!        B(2,1:3*n_nodes) = B(1,1:n_nodes)
!        B(3,1:3*n_nodes) = B(1,1:n_nodes)
!        B(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + dNdy(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + dNdy(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = B(3,3:3*n_nodes:3) + dNdy(1:n_nodes,3)
!
!        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),kstress)*w(kint)*determinant
!        write(*,*) 'residual =', element_residual
!        write(*,*) 'kstress =', kstress
        !read(*,*)
        updated_state_variables(8*(kint-1)+1:8*(kint-1)+6) = kstress(1:6)
        updated_state_variables(8*(kint-1)+7) = eps_b_mat
        updated_state_variables(8*kint) = Vf
        stress(1:6) = kstress(1:6) / JJ
!    write(*,*) 'kstress =', kstress
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
!        write(*,*) 'J =',JJ
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'VF',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
    !write(*,*) 'fieldvars'
  
    return
end subroutine fieldvars_gurson_3dbasic

