!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
!    use Element_Utilities, only:  dNbdx => vol_avg_shape_function_derivatives_3D
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
    integer      :: n_points,kint,i,j,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  kstress(6)
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
!    real (prec)  ::  dB(length_dof_array,6)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu, KK                        ! Material properties
    real (prec)  ::  dof_new(3*n_nodes)
    real (prec)  ::  GG(6,9),Bs(9,3*n_nodes),Iv(6),Bv(6), inv_Bv(6)
    real (prec)  ::  F(3,3), inv_F(3,3), BB(3,3), inv_BB(3,3), JJ
    real (prec)  ::  Iv_dyadic_invBv(6,6),Iv_dyadic_Iv(6,6),Bv_dyadic_invBv(6,6)
    real (prec)  ::  dNdy(n_nodes,3)
    real (prec)  ::  Sigma(3*n_nodes,3*n_nodes),SS(3,length_dof_array/3),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Svec(3*n_nodes), Smat(3*n_nodes,3*n_nodes)
!    real (prec)  :: temp1(3*n_nodes,6),temp2(3*n_nodes,9),temp3(3*n_nodes,3*n_nodes)
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
    mu = element_properties(1)
    KK = element_properties(2)
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = 1.d0
!    D(2,2) = 1.d0
!    D(3,3) = 1.d0
!    D(4,4) = 0.5d0
!    D(5,5) = 0.5d0
!    D(6,6) = 0.5d0
  
!    Vel = 0.d0
!    dNbdx = 0.d0
!    if ( element_identifier == 1) then ! B-bar element
!        do kint = 1, n_points
!            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!            call invert_small(dxdxi,dxidx,determinant)
!            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!
!            dNbdx(1:n_nodes,1:3) = dNbdx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3) * w(kint) * determinant
!            Vel = Vel + w(kint) * determinant
!        end do
!        dNbdx = dNbdx / Vel
!    end if

    !     --  Loop over integration points
    !write(*,*) 'n_points = ',n_points
    do kint = 1, n_points
!        write(*,*) 'kint = ',kint
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dof_new = 0.d0
        dof_new = dof_total + dof_increment

        F = 0.d0
        F(1,1) = 1.d0
        F(2,2) = 1.d0
        F(3,3) = 1.d0
        do i = 1,3 ! Compute deformation gradient tensor F
            do j = 1,3
                do a= 1,n_nodes
                    F(i,j) = F(i,j) +  dof_new(3*(a-1)+i)*dNdx(a,j)
                end do
            end do
        end do
        !WRITE(*,*) 'F'
        dNdy = 0.d0
        call invert_small(F,inv_F,JJ)
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F(1:3,1:3))
!        write(*,*) dNdy
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        !WRITE(*,*) 'B'
!        if(element_identifier==1) then
!        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+1/3*(dNbdy(1:n_nodes,1)-dNdy(1:n_nodes,1))
!        B(1,2:3*n_nodes-1:3)=1/3*(dNbdy(1:n_nodes,2)-dNdy(1:n_nodes,2))
!        B(1,3:3*n_nodes:3)=1/3*(dNbdy(1:n_nodes,3)-dNdy(1:n_nodes,3))
!        B(2,1:3*n_nodes-2:3)=1/3*(dNbdy(1:n_nodes,1) -dNdy(1:n_nodes,1))
!        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+1/3*(dNbdy(1:n_nodes,2) -dNdy(1:n_nodes,2))
!        B(2,3:3*n_nodes:3)=1/3*(dNbdy(1:n_nodes,3) -dNdy(1:n_nodes,3))
!        B(3,1:3*n_nodes-2:3)=1/3*(dNbdy(1:n_nodes,1) -dNdy(1:n_nodes,1))
!        B(3,2:3*n_nodes-1:3)=1/3*(dNbdy(1:n_nodes,2) -dNdy(1:n_nodes,2))
!        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)+1/3*(dNbdy(1:n_nodes,3) -dNdy(1:n_nodes,3))
!        end if
        BB = 0.D0
        BB = matmul(F,transpose(F))
!        write(*,*) BB
        Bv = 0.d0
        Bv = (/ BB(1,1), BB(2,2), BB(3,3), BB(1,2), BB(1,3), BB(2,3) /)
        call invert_small(BB,inv_BB,JJ)
        JJ = dsqrt(JJ)
!        write(*,*) JJ
        inv_Bv = 0.d0
        inv_Bv = (/ inv_BB(1,1), inv_BB(2,2), inv_BB(3,3), inv_BB(1,2), inv_BB(1,3), inv_BB(2,3) /)
!        JJ = dsqrt(BB(1,1)*(BB(2,2)*BB(3,3)-BB(2,3)**2)-BB(1,2)*(BB(2,1)*BB(3,3)-BB(2,3)*BB(3,1)) &
!                        +BB(1,3)*(BB(2,1)*BB(3,2)-BB(3,1)*BB(2,2)))
                        !WRITE(*,*) 'JJ'
        write(*,*) 'J=',JJ

        Iv = 0.d0
        Iv = (/ 1.d0,1.d0,1.d0,0.d0,0.d0,0.d0 /)
        Iv_dyadic_invBv = 0.d0
        Iv_dyadic_invBv = spread(Iv,dim=2,ncopies=6)*spread(inv_Bv,dim=1,ncopies=6)
        Iv_dyadic_Iv = 0.d0
        Iv_dyadic_Iv = spread(Iv,dim=2,ncopies=6)*spread(Iv,dim=1,ncopies=6)
        Bv_dyadic_invBv = 0.d0
        Bv_dyadic_invBv = spread(Bv,dim=2,ncopies=6)*spread(inv_Bv,dim=1,ncopies=6)

        GG = 0.D0
        GG(1,1:9)=(/2.d0*BB(1,1),0.d0,0.d0,2.d0*BB(1,2),0.d0,2.d0*BB(1,3),0.d0,0.d0,0.d0/)
        GG(2,1:9)=(/0.d0,2.d0*BB(2,2),0.d0,0.d0,2.d0*BB(1,2),0.d0,0.d0,2.d0*BB(2,3),0.d0/)
        GG(3,1:9)=(/0.d0,0.d0,2.d0*BB(3,3),0.d0,0.d0,0.d0,2.d0*BB(1,3),0.d0,2.d0*BB(1,3)/)
        GG(4,1:9)=(/2.d0*BB(1,2),2.d0*BB(1,2),0.d0,2.d0*BB(2,2),2.d0*BB(1,1),&
                2.d0*BB(2,3),0.d0,2.d0*BB(1,3),0.d0/)
        GG(5,1:9)=(/2.d0*BB(1,3),0.d0,2.d0*BB(1,3),2.d0*BB(2,3),0.d0,2.d0*BB(3,3)&
                ,2.d0*BB(1,1),0.d0,2.d0*BB(1,2)/)
        GG(6,1:9)=(/0.d0,2.d0*BB(2,3),2.d0*BB(2,3),0.d0,2.d0*BB(1,3),0.d0,2.d0*BB(1,2),&
                2.d0*BB(3,3),2.d0*BB(2,2)/)
!        write(*,*) GG(1,1)
        !WRITE(*,*) 'GG'
        Bs = 0.d0
        Bs(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bs(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bs(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        Bs(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bs(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        Bs(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        Bs(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        Bs(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        Bs(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!        write(*,*) Bs(1,1)
        !WRITE(*,*) 'BS'
        D = 0.d0
        D(1,1) = 1.d0
        D(2,2) = 1.d0
        D(3,3) = 1.d0
        D(4,4) = 0.5d0
        D(5,5) = 0.5d0
        D(6,6) = 0.5d0
        D = mu/JJ**(2.d0/3.d0) * D
        D(1:6,1:6) = D(1:6,1:6) + mu/3.d0/(JJ**(2.D0/3.D0))*((BB(1,1)+BB(2,2)+BB(3,3))/3.d0*Iv_dyadic_invBv - &
                                            Iv_dyadic_Iv - Bv_dyadic_invBv)
        D = D + KK*JJ*(JJ-0.5D0) * Iv_dyadic_invBv
        !WRITE(*,*) 'D'
!        write(*,*) D(1,1)
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        stress = 0.d0
        stress(1) = mu/(JJ**(5.d0/3.d0))*(BB(1,1)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1.d0)
        stress(2) = mu/(JJ**(5.d0/3.d0))*(BB(2,2)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1.d0)
        stress(3) = mu/(JJ**(5.d0/3.d0))*(BB(3,3)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1.d0)
        stress(4) = mu/(JJ**(5.d0/3.d0))*BB(1,2)
        stress(5) = mu/(JJ**(5.d0/3.d0))*BB(1,3)
        stress(6) = mu/(JJ**(5.d0/3.d0))*BB(2,3)
        kstress = 0.d0
        kstress = stress*JJ
        !write(*,*) 'kstress =',kstress
        SS = 0.d0
        SS = reshape(matmul(transpose(B),kstress),(/3,length_dof_array/3/))
        Pvec = 0.d0
        Pmat = 0.d0
        Svec = 0.d0
        Smat = 0.d0
        do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(SS(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do
        Sigma = Pmat*transpose(Smat)
!        write(*,*) Sigma(1,1)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),kstress) &
               *w(kint)*determinant
        !WRITE(*,*) 'RESIDUAL'
        !write(*,*) element_residual(1)
        !temp1 = matmul(transpose(B(1:6,1:3*n_nodes)),D)
!        write(*,*) temp1(1,1)
        !temp2 = matmul(temp1,GG)
!        write(*,*) temp2(1,1)
        !temp3  = matmul(temp2,Bs)
!        write(*,*) temp3(1,1)
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B),matmul(D,matmul(GG, &
                    Bs)))*w(kint)*determinant
!        write(*,*) element_stiffness(1,1)
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            - Sigma(1:3*n_nodes,1:3*n_nodes) *w(kint)*determinant
        !    WRITE(*,*) 'STIFFNESS'

    end do
   write(*,*) 'element stiffness computed'
    return
end subroutine el_hyperelast_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hyperelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
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

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hyperelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
!    use Element_Utilities, only: dNbdx => vol_avg_shape_function_derivatives_3D
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
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,i,j,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
!    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dB(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
!    real (prec)  :: E, xnu, D44, D11, D12, Vel         ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    real (prec)  ::  mu, KK                       ! Material properties
    real (prec)  ::  dof_new(length_dof_array)
    real (prec)  ::   Iv(6),Bv(6), inv_Bv(6)
    real (prec)  ::  F(3,3), inv_F(3,3), BB(3,3), inv_BB(3,3), JJ, temp
!    real (prec)  ::  Iv_dyadic_invBv(6,6),Iv_dyadic_Iv(6,6),Bv_dyadic_invBv(6,6)
    real (prec)  ::  dNdy(n_nodes,3)
!    real (prec)  ::  Sigma(3*n_nodes,3*n_nodes),SS(3,length_dof_array/3),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
!    real (prec)  ::  Svec(3*n_nodes), Smat(3*n_nodes,3*n_nodes)

    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
!    D = 0.d0
    mu = element_properties(1)
    KK = element_properties(2)

!    Vel = 0d0
!    if ( element_identifier == 1) then ! B-bar element
!        do kint = 1, n_points
!            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!            call invert_small(dxdxi,dxidx,determinant)
!            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!
!            dNbdx(1:n_nodes,1:3) = dNbdx(1:n_nodes,1:3) + dNdx(1:n_nodes,1:3) * w(kint) * determinant
!            Vel = Vel + w(kint) * determinant
!        end do
!        dNbdx(1:n_nodes,1:3) = dNbdx(1:n_nodes,1:3) / Vel
!    end if

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        dof_new = dof_total + dof_increment

        F = 0.d0
        F(1,1) = 1.d0
        F(2,2) = 1.d0
        F(3,3) = 1.d0
        do i = 1,3 ! Compute deformation gradient tensor F
            do j = 1,3
                do a= 1,n_nodes
                    F(i,j) = F(i,j) +  dof_new(3*(a-1)+i)*dNdx(a,j)
                end do
            end do
        end do
        !WRITE(*,*) 'F'
        dNdy = 0.d0
        call invert_small(F,inv_F,JJ)
        dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F(1:3,1:3))
!        write(*,*) dNdy
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        !WRITE(*,*) 'B'
!        if(element_identifier==1) then
!        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+1/3*(dNbdy(1:n_nodes,1)-dNdy(1:n_nodes,1))
!        B(1,2:3*n_nodes-1:3)=1/3*(dNbdy(1:n_nodes,2)-dNdy(1:n_nodes,2))
!        B(1,3:3*n_nodes:3)=1/3*(dNbdy(1:n_nodes,3)-dNdy(1:n_nodes,3))
!        B(2,1:3*n_nodes-2:3)=1/3*(dNbdy(1:n_nodes,1) -dNdy(1:n_nodes,1))
!        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+1/3*(dNbdy(1:n_nodes,2) -dNdy(1:n_nodes,2))
!        B(2,3:3*n_nodes:3)=1/3*(dNbdy(1:n_nodes,3) -dNdy(1:n_nodes,3))
!        B(3,1:3*n_nodes-2:3)=1/3*(dNbdy(1:n_nodes,1) -dNdy(1:n_nodes,1))
!        B(3,2:3*n_nodes-1:3)=1/3*(dNbdy(1:n_nodes,2) -dNdy(1:n_nodes,2))
!        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)+1/3*(dNbdy(1:n_nodes,3) -dNdy(1:n_nodes,3))
!        end if
        BB = matmul(F,transpose(F))
!        write(*,*) BB
        Bv = (/ BB(1,1), BB(2,2), BB(3,3), BB(1,2), BB(1,3), BB(2,3) /)
        call invert_small(BB,inv_BB,JJ)
        inv_Bv = (/ inv_BB(1,1), inv_BB(2,2), inv_BB(3,3), inv_BB(1,2), inv_BB(1,3), inv_BB(2,3) /)
        JJ = dsqrt(JJ)
        !JJ = dsqrt(BB(1,1)*(BB(2,2)*BB(3,3)-BB(2,3)**2)-BB(1,2)*(BB(2,1)*BB(3,3)-BB(2,3)*BB(3,1)) &
        !                +BB(1,3)*(BB(2,1)*BB(3,2)-BB(3,1)*BB(2,2)))
                        !WRITE(*,*) 'JJ'

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        stress(1) = mu/JJ**(5.d0/3.d0)*(BB(1,1)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1)
        stress(2) = mu/JJ**(5.d0/3.d0)*(BB(2,2)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1)
        stress(3) = mu/JJ**(5.d0/3.d0)*(BB(3,3)-(BB(1,1)+BB(2,2)+BB(3,3))/3.d0) + KK*(JJ-1)
        stress(4) = mu/JJ**(5.d0/3.d0)*BB(1,2)
        stress(5) = mu/JJ**(5.d0/3.d0)*BB(1,3)
        stress(6) = mu/JJ**(5.d0/3.d0)*BB(2,3)
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
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
            endif
        end do
 
    end do
   write(*,*) 'field vars finished'
    return
end subroutine fieldvars_hyperelast_3dbasic

