!
! This code is written for the peridynamic differential operator to solve the Timoshenko's two-dimensional cantilevered beam problem given in Section 5.2.
!
!
! Dorduncu M, Ren H, Zhuang X, Silling S, Madenci E, Rabczuk T. "A review of peridynamic theory and nonlocal operators 
! along with their computer implementations." Computers & Structures 2024;299:107395. https://doi.org/10.1016/j.compstruc.2024.107395.
!
! Input files are PDnav2D.inp and PDgeom2D.dat (subroutine 'input' is used to read the files)
!
! PDnav2D.inp - provide the details of
!      -   name of input file for geometry
!      -   Order of Taylor series expansion for PD differential operator used to evaluate deformation gradient
!      -   Specify the horizon size
!      -   specify Youngs Modulus , Poisson's ratio, and coefficient of thermal expansion
!      -   Number of pre-existing cracks
!      -   location of pre-existing cracks
!      -   Number  of boundaries of geometry
!      -   Min and Max coordinates of the boundary along with normals associated with the boundaries
!      -   Number of Dirichlet boundary conditions
!      -   Details of Dirichlet boundary conditions (coordinates, direction and value)
!      -   Number of stress boundary conditions
!      -   Details of stress boundary conditions (coordinates (x1, x2, y1, y2), stress values)
!   
!  The construction of stiffness matrix and load vector for body points is provided in suboutine 'SetupODEMatrixVector2D'.
!   
!  The construction of stiffness matrix and load vector for boundary points is provided in suboutine 'SetupODEBoundaryConditions2D'.
!
!  The program uses MKL library available in Visual Studio along with Inter Pardiso  for the solution of the system of equations.
!
!
!
module dbase
implicit none
integer totnode
!
real *8, allocatable :: coord(:,:)
integer, allocatable :: numfam(:)
!
    Type material_library
       real *8 dens
       real *8 emod1
       real *8 emod2
       real *8 emod3
       real *8 pratio12
       real *8 pratio21
       real *8 g12
       real *8 g13
    end type material_library
    Type(material_library) matlib(20)
!
    Type DiffOperator_
         integer n1, n2, n3
         real *8 val, coef
    end type DiffOperator_
!
    Type BoundaryCondition
         integer num_diff_ops_u, num_diff_ops_v, u_flag, v_flag, sxx_flag, syy_flag, sxy_flag
         integer bx_flag, by_flag
         Type(DiffOperator_), allocatable :: uDiffs(:), vDiffs(:)
         real *8 val, xmin, xmax, ymin, ymax, zmin, zmax
    end type BoundaryCondition
!    
    Type Equation
         Type(DiffOperator_), allocatable :: uDiffs(:), vDiffs(:)
         integer num_diff_ops_u, num_diff_ops_v
         real *8 RHSval
    end type Equation
!    
    Type PDoperator_
         Type(Equation) Eq1, Eq2
         Type(DiffOperator_), allocatable :: u_out(:), v_out(:)
         integer norder, n1order, n2order, n3order, nsize, morder, asymFlag
         integer num_bc, nteqs, ncons, atype, num_out_u, num_out_v, nwk, nslits
         Type(BoundaryCondition), allocatable :: bc(:)
         Type(DiffOperator_), allocatable :: order(:)
         character *80 fname
         real *8 lambda, mu, Emod, nu, Gmod, C11, C22, C12, C66
    end type PDoperator_
    Type(PDoperator_) PDoperator
!
    Type nodefam_
         integer , allocatable :: node(:)
         integer , allocatable :: nodeFlag(:)
    end type nodefam_    
!
    Type SlitLine
         real *8 x1, y1, x2, y2
    end type SlitLine
    Type(SlitLine) slits(10)     
!
integer atype
!
real *8 length
real *8 width
real *8 hole_radius
!
Type(nodefam_), allocatable :: family(:)
real *8, allocatable :: DiffAmat3D(:,:), DiffAmatInv3D(:,:), DiffBvec3D(:), DiffAvec(:), plist(:), blist(:), weight(:)
real *8, allocatable :: DiffAmat2D(:,:), DiffAmatInv2D(:,:), DiffBvec2D(:)
real *8, allocatable :: DiffAmat1D(:,:), DiffAmatInv1D(:,:), DiffBvec1D(:)
real *8, allocatable :: rcvec(:), fvec(:), dvolume(:), deltax(:), deltay(:), deltaz(:)
real *8, allocatable :: SpSysMat(:), SysMat(:,:), SysVec(:), Coefs(:), sysnorm(:), SpColMax(:), bforce_x(:), bforce_y(:)
integer, allocatable :: irow(:)
real *8, allocatable :: icol(:)

end module dbase
!
!
!
      module GlobalVariables
!
      integer, parameter :: node_above=1, node_below=2, ndof=3, nload=4, ndiv_impactor=30
!
      real *8 dx, dy, dz, area
      real *8, parameter :: pi=dacos(-1.0d0)
!
      integer ilay, ijply
      real *8 pratio12
!
      real *8 emod1, emod2, emod3, g12, g13, ply_theta, ply_thick, dens
      real *8 zloc
      integer ndvxy
      real *8 scale
      character *6 strNum
      character *80 restart_file_name
      character *50 cskip
      real *8 wall1, wall0, twall
      end module GlobalVariables
!  
!
!  
      program main
      use dbase
      use GlobalVariables
      implicit none
      integer i, j, k, ii, n_impactors, ij
      real *8 TIME
      real *8 OMP_GET_WTIME
      integer START_CLOCK,END_CLOCK,ClOCK_RATE
      real *8, allocatable :: matF(:,:,:)
!
!
      call system('rd results /s /q' )
      call system('mkdir results' )
!
!
      call input()
!
!
      call AllocateAndInitializeArrays()
!
      call GenerateNodeFamilies()
      call ApplySlit()
!
      ij = 0
      call SetupODEMatrixVector2D(ij)
      call SetupODEBoundaryConditions2D(ij)
      call SparseSolve()
      do i = 1 , 2*totnode
          fvec(i) = SysVec(i)
      enddo
      call CalcStressStrain()
!
      end program main
!
!
!
      subroutine input()
      use dbase
      use GlobalVariables
      implicit none
      character line
      integer i,j, k, ii, imat, idiff, ibc
      real *8 radius, psi1, psi2, theta1, theta2
      integer mR, mPsi, mTheta, nPsi, nTheta
      real *8 Lx, Ly, Lz
      integer ndivx, ndivy, mdivx, mdivy
      integer ndivR, ndivT, ndivZ, mdivR, mdivT, mdivZ
      real *8 critical_stretch_fiber, critical_stretch_matrix, critical_stretch_fiber_comp, critical_stretch_matrix_comp
      real *8 critical_stretch_peel, critical_stretch_shear
      integer itarget, iout, islit
      real *8 lambda, mu
      character *80 fname
!
      open(1,file='PDnav2D.inp')
!
! 
      read(1,'(A)') line
      read(1,*) PDoperator%fname
      open(2,file=trim(adjustl(PDoperator%fname)))
      read(1,'(A)') line
      read(1,*) PDoperator%n1order, PDoperator%n2order
      PDoperator%atype = 1
      PDoperator%morder = 2
      PDoperator%n3order = 0
      PDoperator%asymFlag = 0
      read(1,'(A)') line
      read(1,*) PDoperator%Emod, PDoperator%nu
      PDoperator%lambda = PDoperator%Emod*PDoperator%nu/(1.0d0+PDoperator%nu)/(1.0d0-2.0d0*PDoperator%nu)
      PDoperator%mu = PDoperator%Emod/2.0d0/(1.0d0+PDoperator%nu)
      lambda = PDoperator%lambda
      mu = PDoperator%mu
      PDoperator%C11 = PDoperator%Emod/(1-PDoperator%nu**2)
      PDoperator%C22 = PDoperator%Emod/(1-PDoperator%nu**2)
      PDoperator%C12 = PDoperator%nu*PDoperator%Emod/(1-PDoperator%nu**2)
      PDoperator%C66 = PDoperator%mu
!
      read(2,*) totnode
      allocate(fvec(2*totnode))
      allocate(bforce_x(totnode))
      allocate(bforce_y(totnode))
      allocate(coord(totnode,3))
      allocate(deltax(totnode))
      allocate(deltay(totnode))
      allocate(deltaz(totnode))
      allocate(dvolume(totnode))
      do i = 1 , totnode
         read(2,*) coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), bforce_x(i), bforce_y(i)
      enddo
      close(2)
!
      call getSize2D(PDoperator%n1order, PDoperator%n2order, PDoperator%nsize)
!        
      PDoperator%Eq1%num_diff_ops_u = 2
      PDoperator%Eq1%num_diff_ops_v = 1
      PDoperator%Eq2%num_diff_ops_u = 1
      PDoperator%Eq2%num_diff_ops_v = 2
      allocate(PDoperator%Eq1%uDiffs(PDoperator%Eq1%num_diff_ops_u))
      allocate(PDoperator%Eq1%vDiffs(PDoperator%Eq1%num_diff_ops_v))
      allocate(PDoperator%Eq2%uDiffs(PDoperator%Eq2%num_diff_ops_u))
      allocate(PDoperator%Eq2%vDiffs(PDoperator%Eq2%num_diff_ops_v))
!      
      PDoperator%Eq1%uDiffs(1)%n1 = 2
      PDoperator%Eq1%uDiffs(1)%n2 = 0
      PDoperator%Eq1%uDiffs(1)%n3 = 0
      PDoperator%Eq1%uDiffs(1)%coef = PDoperator%mu + PDoperator%Emod/2.0d0/(1.0d0-PDoperator%nu)
!          
      PDoperator%Eq1%uDiffs(2)%n1 = 0
      PDoperator%Eq1%uDiffs(2)%n2 = 2
      PDoperator%Eq1%uDiffs(2)%n3 = 0
      PDoperator%Eq1%uDiffs(2)%coef = PDoperator%mu
!      
      
      PDoperator%Eq1%vDiffs(1)%n1 = 1
      PDoperator%Eq1%vDiffs(1)%n2 = 1
      PDoperator%Eq1%vDiffs(1)%n3 = 0
      PDoperator%Eq1%vDiffs(1)%coef = PDoperator%Emod/2.0d0/(1.0d0-PDoperator%nu)
!      
      PDoperator%Eq2%uDiffs(1)%n1 = 1
      PDoperator%Eq2%uDiffs(1)%n2 = 1
      PDoperator%Eq2%uDiffs(1)%n3 = 0
      PDoperator%Eq2%uDiffs(1)%coef = PDoperator%Emod/2.0d0/(1.0d0-PDoperator%nu)

!      
      PDoperator%Eq2%vDiffs(1)%n1 = 2
      PDoperator%Eq2%vDiffs(1)%n2 = 0
      PDoperator%Eq2%vDiffs(1)%n3 = 0
      PDoperator%Eq2%vDiffs(1)%coef = PDoperator%mu
!
      PDoperator%Eq2%vDiffs(2)%n1 = 0
      PDoperator%Eq2%vDiffs(2)%n2 = 2
      PDoperator%Eq2%vDiffs(2)%n3 = 0
      PDoperator%Eq2%vDiffs(2)%coef = PDoperator%mu + PDoperator%Emod/2.0d0/(1.0d0-PDoperator%nu)

!
      read(1,'(A)') line
      read(1,*) PDoperator%nslits
      if(PDoperator%nslits>0) then
          read(1,'(A)') line
          do islit = 1 , PDoperator%nslits
             read(1,*) slits(islit)%x1, slits(islit)%y1, slits(islit)%x2, slits(islit)%y2
          enddo
      endif      
!      
      read(1,'(A)') line
      read(1,*) PDoperator%num_bc
      if( PDoperator%num_bc > 0) then
          allocate(PDoperator%bc(PDoperator%num_bc))
          do ibc = 1 , PDoperator%num_bc
             read(1,'(A)') line
             read(1,*) PDoperator%bc(ibc)%xmin, PDoperator%bc(ibc)%xmax, PDoperator%bc(ibc)%ymin, PDoperator%bc(ibc)%ymax, PDoperator%bc(ibc)%u_flag, PDoperator%bc(ibc)%v_flag, &
                       PDoperator%bc(ibc)%sxx_flag, PDoperator%bc(ibc)%syy_flag, PDoperator%bc(ibc)%sxy_flag, PDoperator%bc(ibc)%bx_flag, PDoperator%bc(ibc)%by_flag, PDoperator%bc(ibc)%val
             if( PDoperator%bc(ibc)%u_flag > 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 0
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = 1.0d0
             elseif( PDoperator%bc(ibc)%u_flag < 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = PDoperator%bc(ibc)%u_flag
                 PDoperator%bc(ibc)%num_diff_ops_v = 0
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = 1.0d0
             else if( PDoperator%bc(ibc)%v_flag > 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 0
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%vDiffs(1).coef = 1.0d0
             else if( PDoperator%bc(ibc)%v_flag < 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 0
                 PDoperator%bc(ibc)%num_diff_ops_v = PDoperator%bc(ibc)%v_flag
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%vDiffs(1).coef = 1.0d0
             else if( PDoperator%bc(ibc)%sxx_flag > 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C11
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C12
             else if( PDoperator%bc(ibc)%syy_flag > 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C12
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C22
             else if( PDoperator%bc(ibc)%sxy_flag > 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C66
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C66
             else if( PDoperator%bc(ibc)%sxx_flag < 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C11
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C12
             else if( PDoperator%bc(ibc)%syy_flag < 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C12
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C22
             else if( PDoperator%bc(ibc)%sxy_flag < 0 ) then
                 PDoperator%bc(ibc)%num_diff_ops_u = 1
                 PDoperator%bc(ibc)%num_diff_ops_v = 1
                 allocate(PDoperator%bc(ibc)%uDiffs(1))
                 allocate(PDoperator%bc(ibc)%vDiffs(1))
                 PDoperator%bc(ibc)%uDiffs(1).n1 = 0
                 PDoperator%bc(ibc)%uDiffs(1).n2 = 1
                 PDoperator%bc(ibc)%uDiffs(1).coef = PDoperator%C66
                 PDoperator%bc(ibc)%vDiffs(1).n1 = 1
                 PDoperator%bc(ibc)%vDiffs(1).n2 = 0
                 PDoperator%bc(ibc)%vDiffs(1).coef = PDoperator%C66
             endif
          enddo
      endif
!
      read(1,'(A)') line
      read(1,*) PDoperator%num_out_u
      if(PDoperator%num_out_u > 0 ) then
         allocate(PDoperator%u_out(PDoperator%num_out_u))
         read(1,'(A)') line
         read(1,*) (PDoperator%u_out(iout)%n1, PDoperator%u_out(iout)%n2, iout = 1 , PDoperator%num_out_u)
      endif
      read(1,'(A)') line
      read(1,*) PDoperator%num_out_v
      if( PDoperator%num_out_v > 0 ) then
         allocate(PDoperator%v_out(PDoperator%num_out_v))
         read(1,'(A)') line
         read(1,*) (PDoperator%v_out(iout)%n1, PDoperator%v_out(iout)%n2, iout = 1 , PDoperator%num_out_v) 
      endif
!
      close(1)
!
      return
      end
!
!
!
      subroutine AllocateAndInitializeArrays()
      use dbase
      use GlobalVariables
      implicit none
      integer i, j, k, ii
!
      PRINT *,'AllocateAndInitializeArrays'
!
      allocate(numfam(totnode))
      allocate(family(totnode))
      !
!
      do i = 1, totnode
         numfam(i) = 0
      enddo
!
      return
      end subroutine AllocateAndInitializeArrays
!
!
!
      subroutine GenerateNodeFamilies()
      use dbase
      use GlobalVariables
      implicit none
      integer, parameter :: maxFamNodes=10000
      real *8, parameter :: tol=1.0d-5
      integer i, j, k, ii
      real *8 ijdist, xdist, ydist, zdist, delta_mag, xyz_mag
      integer :: nodefam(maxFamNodes)
!
      PRINT *,'GenerateNodeFamilies'
!
!
      do i = 1,totnode
          delta_mag = dsqrt(deltax(i)**2+deltay(i)**2+deltaz(i)**2)
          numfam(i) = 1
          nodefam(1) = i
	      do j = 1,totnode
	         if(j==i) cycle
	         xdist = dabs(coord(j,1) - coord(i,1)) 
	         ydist = dabs(coord(j,2) - coord(i,2)) 
	         zdist = dabs(coord(j,3) - coord(i,3)) 
	         xyz_mag = dsqrt(xdist**2+ydist**2+zdist**2)
	         if( xdist <= deltax(i) .and. ydist <= deltay(i) .and. zdist <= deltaz(i)) then
	             if( (PDoperator.asymFlag == 0) .or. &
	                 ( PDoperator.asymFlag == 1 .and. coord(j,1) <= coord(i,1)+tol ) .or. &
	                 ( PDoperator.asymFlag == 2 .and. coord(j,2) <= coord(i,2)+tol ) .or. &
	                 ( PDoperator.asymFlag == 3 .and. coord(j,3) <= coord(i,3)+tol ) ) then
	                   numfam(i) = numfam(i) + 1
	                   if( numfam(i) > maxFamNodes )then
	                       print *, 'numfam(i) > maxFamNodes'
	                       stop
	                   else
	                       nodefam(numfam(i)) = j
	                   endif
	              endif
	         endif  
	      enddo
          !print *, 'allocate node family arrays'
	      allocate( family(i)%node(numfam(i)))
	      allocate( family(i)%nodeFlag(numfam(i)))
	      do j = 1 , numfam(i)
	         family(i)%node(j) = nodefam(j)
	         family(i)%nodeFlag(j) = 1
	      enddo
	  enddo
!
	  return
!
	  end subroutine GenerateNodeFamilies
!
!
!
    subroutine SetupODEMatrixVector2D(ij)
    use dbase
    use GlobalVariables
    implicit none
    integer i, j, k, l, jmem, morder, n1order, n2order,n3order,nsize, n1, n2, m, idiff
    integer ibc, icons, ij, ij_rel
    real *8 coef, gfunval, xsi1, xsi2, RHSVal, val1, val2, val3, val4, val5, val6, val7
    real *8 delta_mag, rcond, max_rcond
    integer ival1, jflag
!    
!    
    n1order = PDoperator%n1order
    n2order = PDoperator%n2order
    n3order = 0
    nsize = PDoperator%nsize
    morder = PDoperator%morder
!
    call getSysMatIndex()  
    call SetDiffOperators2D(n1order, n2order, nsize)
    do i = 1 , PDoperator%nwk
       SpSysMat(i) = 0.0d0
    enddo
    do i = 1 , PDoperator%nteqs
       SysVec(i) = 0.0d0
    enddo
!
    allocate(DiffAmat2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffAmatInv2D(PDoperator%nsize,PDoperator%nsize))
    allocate(DiffBvec2D(PDoperator%nsize))
    allocate(DiffAvec(PDoperator%nsize))
    allocate(rcvec(PDoperator%nsize))
    allocate(plist(PDoperator%nsize))
    allocate(blist(PDoperator%nsize))
    allocate(weight(PDoperator%nsize))
!    
    open(12,file=trim(adjustl(PDoperator%fname)))
    read(12,*) ival1
!    
    ij_rel = 0
    do k = 1 , totnode
       call FormDiffAmat2D(n1order, n2order, nsize, k)
       call recondition(DiffAmat2D,rcvec,nsize)
       call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)

       delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
       PDoperator%eq1%RHSval = bforce_x(k)
       PDoperator%eq2%RHSval = bforce_y(k)

!      1st equation
       do idiff = 1 , PDoperator%eq1%num_diff_ops_u
          n1 = PDoperator%eq1%uDiffs(idiff)%n1
          n2 = PDoperator%eq1%uDiffs(idiff)%n2
          coef = PDoperator%eq1%uDiffs(idiff)%coef
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             SpSysMat(ij_rel + jmem) = SpSysMat(ij_rel + jmem) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
             SpSysMat(ij_rel + 1)    = SpSysMat(ij_rel + 1)    - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
          enddo
       enddo
       ij_rel = ij_rel + numfam(k)
       do idiff = 1 , PDoperator%eq1%num_diff_ops_v
          n1 = PDoperator%eq1%vDiffs(idiff)%n1
          n2 = PDoperator%eq1%vDiffs(idiff)%n2
          coef = PDoperator%eq1%vDiffs(idiff)%coef
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             SpSysMat(ij_rel + jmem) = SpSysMat(ij_rel + jmem) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
             SpSysMat(ij_rel + 1)    = SpSysMat(ij_rel + 1)    - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
          enddo
       enddo
       SysVec(2*k-1) = PDoperator%eq1%RHSval
       ij_rel = ij_rel + numfam(k)
!
!      2nd equation
       do idiff = 1 , PDoperator%eq2%num_diff_ops_u
          n1 = PDoperator%eq2%uDiffs(idiff)%n1
          n2 = PDoperator%eq2%uDiffs(idiff)%n2
          coef = PDoperator%eq2%uDiffs(idiff)%coef
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             SpSysMat(ij_rel + jmem) = SpSysMat(ij_rel + jmem) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
             SpSysMat(ij_rel + 1)    = SpSysMat(ij_rel + 1)    - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
          enddo
       enddo
       ij_rel = ij_rel + numfam(k)
       do idiff = 1 , PDoperator%eq2%num_diff_ops_v
          n1 = PDoperator%eq2%vDiffs(idiff)%n1
          n2 = PDoperator%eq2%vDiffs(idiff)%n2
          coef = PDoperator%eq2%vDiffs(idiff)%coef
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          !call ApplyConstraintOnBvec(DiffBvec2D, PDoperator%nsize, k)
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             SpSysMat(ij_rel + jmem) = SpSysMat(ij_rel + jmem) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
             SpSysMat(ij_rel + 1)    = SpSysMat(ij_rel + 1)    - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
          enddo
       enddo
       SysVec(2*k) = PDoperator%eq2%RHSval
       ij_rel = ij_rel + numfam(k)
    enddo
    ij = ij_rel
    close(12)
    return
    end
!
!
!
    subroutine inverse2(nz, mat, invmat, i, rcond)
    implicit none
    real *8 mat(nz,nz), invmat(nz,nz)
    integer, allocatable :: ipiv(:), iwork(:)
    real *8, allocatable :: work(:)
    integer info, kz, k2, k1, k, i, nz
    real *8 anorm, a, rcond
!    
    allocate (ipiv(nz), work(4*nz), iwork(nz))
    call dgetrf( nz, nz, mat, nz, ipiv, info )   
    anorm=0.0d0
    do k2=1, nz
       a=0.0d0
       do k1=1, nz
          a=a+dabs(mat(k1,k2))
       enddo
       anorm=max(anorm,a)
    enddo
    call dgecon( '1', nz, mat, nz, anorm, rcond, work, iwork, info )
    rcond=1.0d0/rcond
    if( rcond > 1.0d8 .or. rcond .ne. rcond )then
        print *, 'E R R O R: ill-conditioned for k=', i, 'rcond=', rcond
        !read(*,*)
        !stop
    endif
    call dgetri( nz, mat, nz, ipiv, work, 4*nz, info ) 
    do k=1, nz
       !invmat(k,:) = fact(k) * ds**sum(order(k,:)) * mat(k,:)
       invmat(k,:) = mat(k,:)
    enddo
    deallocate (ipiv)
    deallocate(work)
    deallocate(iwork)
    return
    end
!
!
!
    subroutine getSysMatIndex()
    use dbase
    implicit none
    integer nwk, i, k, j, ibc, jmem, icons, row_count, nmask
    integer, allocatable :: mask(:)
!    
    allocate(mask(2*totnode))
    nwk = 0
    do k = 1 , totnode
       nwk = nwk + 2*numfam(k)
       nwk = nwk + 2*numfam(k)
    enddo
    icons = 0
    do ibc = 1 , PDoperator%num_bc
       if(PDoperator%bc(ibc)%u_flag > 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             nwk = nwk + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%u_flag < 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             nwk = nwk + 4
          enddo
       elseif(PDoperator%bc(ibc)%v_flag > 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             nwk = nwk + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%v_flag < 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             nwk = nwk + 4
          enddo
       elseif(PDoperator%bc(ibc)%sxx_flag > 0 .or. PDoperator%bc(ibc)%syy_flag > 0 .or. PDoperator%bc(ibc)%sxy_flag > 0  ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             nwk = nwk + 4*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%sxx_flag < 0 .or. PDoperator%bc(ibc)%syy_flag < 0 .or. PDoperator%bc(ibc)%sxy_flag < 0  ) then
          icons = icons + 1
          do k = 1 , 2*totnode
             mask(k) = 0
          enddo
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                mask(2*j-1) = 1
                mask(2*j) = 1
             enddo
          enddo
          nmask = 0
          do k = 1 , 2*totnode
             nmask = nmask + mask(k)
          enddo
          nwk = nwk + 2*nmask
       endif
    enddo
!        
    PDoperator%ncons = icons
    PDoperator%nteqs = 2*totnode + PDoperator%ncons
    PDoperator%nwk = nwk
    print *, 'before:', nwk
!
    allocate(icol(PDoperator%nwk))
    allocate(irow(PDoperator%nteqs+1))
    allocate(SpSysMat(PDoperator%nwk))
    allocate(SysVec(PDoperator%nteqs))
    allocate(SpColMax(PDoperator%nteqs))
!
    do i = 1 , PDoperator%nteqs
       irow(i) = 0
    enddo
    do i = 1 , PDoperator%nwk
       icol(i) = 0
    enddo
!
    nwk = 0
    do k = 1 , totnode
       do jmem = 1 , numfam(k)
          j = family(k).node(jmem)
          irow(2*k-1) = irow(2*k-1) + 1
          icol(nwk + jmem) = PDoperator%nteqs*((2*k-1)-1) + 2*j - 1
       enddo
       nwk=nwk+numfam(k)
       do jmem = 1 , numfam(k)
          j = family(k).node(jmem)
          irow(2*k-1) = irow(2*k-1) + 1
          icol(nwk+jmem) = PDoperator%nteqs*((2*k-1)-1) + 2*j
       enddo
       nwk=nwk+numfam(k)
       do jmem = 1 , numfam(k)
          j = family(k).node(jmem)
          irow(2*k) = irow(2*k) + 1
          icol(nwk+jmem) = PDoperator%nteqs*((2*k)-1) + 2*j - 1
       enddo
       nwk=nwk+numfam(k)
       do jmem = 1 , numfam(k)
          j = family(k).node(jmem)
          irow(2*k) = irow(2*k) + 1
          icol(nwk+jmem) = PDoperator%nteqs*((2*k)-1) + 2*j
       enddo
       nwk = nwk + numfam(k)
    enddo
    icons = 0
    do ibc = 1 , PDoperator%num_bc
       if(PDoperator%bc(ibc)%u_flag > 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                irow(2*totnode+icons) = irow(2*totnode+icons) + 1
                icol(nwk+2*jmem-1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*j - 1
                irow(2*j-1) = irow(2*j-1) + 1
                icol(nwk+2*jmem) = PDoperator%nteqs*((2*j-1)-1) + 2*totnode + icons
             enddo
             nwk=nwk + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%u_flag < 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             irow(2*totnode+icons) =  irow(2*totnode+icons) + 2
             icol(nwk+1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*k-1
             icol(nwk+2) = PDoperator%nteqs*(2*totnode+icons-1) + 2*(-PDoperator%bc(ibc)%num_diff_ops_u) - 1
             irow(2*k-1) =  irow(2*k-1) + 1
             icol(nwk+3) = PDoperator%nteqs*(2*k-2) + 2*totnode+icons
             irow(2*(-PDoperator%bc(ibc)%num_diff_ops_u)-1) =  irow(2*(-PDoperator%bc(ibc)%num_diff_ops_u)-1) + 1
             icol(nwk+4) = PDoperator%nteqs*(2*(-PDoperator%bc(ibc)%num_diff_ops_u)-2) + 2*totnode + icons
             nwk=nwk+4
          enddo
       elseif(PDoperator%bc(ibc)%v_flag > 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                irow(2*totnode+icons) = irow(2*totnode+icons) + 1
                icol(nwk+2*jmem-1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*j
                irow(2*j) = irow(2*j) + 1
                icol(nwk+2*jmem) = PDoperator%nteqs*((2*j)-1) + 2*totnode + icons
             enddo
             nwk=nwk + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%v_flag < 0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             irow(2*totnode+icons) =  irow(2*totnode+icons) + 2
             icol(nwk+1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*k
             icol(nwk+2) = PDoperator%nteqs*(2*totnode+icons-1) + 2*(-PDoperator%bc(ibc)%num_diff_ops_v)
             irow(2*k) = irow(2*k) + 1
             icol(nwk+3) = PDoperator%nteqs*(2*k-1) + 2*totnode+icons
             irow(2*(-PDoperator%bc(ibc)%num_diff_ops_v)) = irow(2*(-PDoperator%bc(ibc)%num_diff_ops_v)) + 1
             icol(nwk+4) = PDoperator%nteqs*(2*(-PDoperator%bc(ibc)%num_diff_ops_v)-1) + 2*totnode+icons
             nwk=nwk+4
          enddo
       elseif(PDoperator%bc(ibc)%sxx_flag >0 .or. PDoperator%bc(ibc)%syy_flag > 0 .or. PDoperator%bc(ibc)%sxy_flag >0 ) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                irow(2*totnode+icons) = irow(2*totnode+icons) + 1
                icol(nwk+2*jmem-1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*j - 1
                irow(2*j-1) = irow(2*j-1) + 1
                icol(nwk+2*jmem) = PDoperator%nteqs*((2*j-1)-1) + 2*totnode + icons
             enddo
             nwk=nwk + 2*numfam(k)
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                irow(2*totnode+icons) = irow(2*totnode+icons) + 1
                icol(nwk+2*jmem-1) = PDoperator%nteqs*(2*totnode+icons-1) + 2*j
                irow(2*j) = irow(2*j) + 1
                icol(nwk+2*jmem) = PDoperator%nteqs*((2*j)-1) + 2*totnode + icons
             enddo
             nwk=nwk + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%sxx_flag <0 .or. PDoperator%bc(ibc)%syy_flag <0 .or. PDoperator%bc(ibc)%sxy_flag <0 ) then
          icons = icons + 1
          do k = 1 , 2*totnode
             mask(k) = 0
          enddo
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                mask(2*j-1) = 1
                mask(2*j) = 1
             enddo
          enddo
          nmask = 0
          do k = 1 , 2*totnode
             if( mask(k) > 0 ) then
                 nmask = nmask + 1
                 mask(k) = nmask
             endif
          enddo
          irow(2*totnode+icons) = nmask
          do k = 1 , 2*totnode
             if(mask(k)>0) then
                icol(nwk+mask(k)) = PDoperator%nteqs*(2*totnode+icons-1) + k
                irow(k)=irow(k) + 1
                icol(nwk+nmask+mask(k)) = PDoperator%nteqs*(k-1) + 2*totnode+icons
             endif                
          enddo
          nwk = nwk + 2*nmask
       endif
    enddo
!      
    nwk = 0  
    do i = 1 , PDoperator%nteqs
       row_count = irow(i)
       irow(i) = nwk + 1
       nwk = nwk + row_count
    enddo
    irow(PDoperator%nteqs+1) = nwk+1
    print *, 'after:',nwk
    deallocate(mask)
!           
    return
    end    
!
!
!
    subroutine SparseSolve()
    use dbase
    include 'mkl_rci.fi'
    integer *8, allocatable :: key(:)
    real *8, allocatable :: mat2(:), dvec(:), key2(:)
    
    integer *8, allocatable :: col2(:)
    integer , allocatable :: col(:)
    integer *8 j

    allocate(dvec(PDoperator%nteqs))
    allocate(key(PDoperator%nwk),mat2(PDoperator%nwk),col2(PDoperator%nwk),col(PDoperator%nwk))
    do ii = 1 , PDoperator%nwk
       mat2(ii)=SpSysMat(ii)
!       col2(ii)=real(icol(ii),8)
    enddo
    do ii=1, PDoperator%nwk
       key(ii)=ii
    enddo

    PRINT *,'Sorting begins'
!    call dlasrt2('i', PDoperator%nwk, col2, key, info)
!    call quick_sort(PDoperator%nwk, col2, key)
!    call dsort(PDoperator%nwk,PDoperator%nwk,key,col2)
!    call sort(PDoperator%nwk,col2,key2)
    call quicksort(icol,key,1,PDoperator%nwk)
    col2(:)= int(icol(:),8)

!    open(14,file='col2.dat')
!    do j = 1,PDoperator%nwk
!    write(14,"(i15,2x,i15)") col2(j), key(j)
!    enddo
!    close(14)

!    key(:)=int(key2(:))
    PRINT *,'Sorting - Done!'
    do ii=1, PDoperator%nwk
       SpSysMat(ii)=mat2(key(ii))
    enddo
!    deallocate(key,mat2,col2)
!
!print *, 'row'     
    j=0
    do ii=1, PDoperator%nwk
       if( col2(ii) > j )then
           j=j+PDoperator%nteqs
       endif
       col2(ii)=col2(ii)-(j-PDoperator%nteqs)
    enddo
!
    col(:) = int(col2(:))
    
!    open(14,file='col2.dat')
!    do j = 1,PDoperator%nwk
!    write(14,"(i15,2x,i15)") col(j), col2(j)
!    enddo
!    close(14)

    print *, 'Solve started'
    dvec(:)=0.0d0
    call intel_direct_solver(PDoperator%nteqs,PDoperator%nwk,irow,col,SpSysMat,dvec,SysVec)
    print *, 'Solve ended'
    do i = 1 , PDoperator%nteqs
       SysVec(i) = dvec(i)
!       SysVec(i) = dvec(i)
    enddo
    deallocate(dvec)
!
    return
    end
!
!
!
    subroutine SparseRecondition()
    use dbase
    implicit none
    integer i, j, k
!    
    do i = 1 , PDoperator%nteqs
       SpColMax(i) = 0.0d0
    enddo
    i = 1
    do k = 1 , PDoperator.nwk
       j = icol(k)
       if( dabs(SpSysMat(k)) > SpColMax(j) ) SpColmax(j) = dabs(SpSysMat(k))
    enddo
    do k = 1 , PDoperator.nwk
       j = icol(k)
       SpSysMat(k) = SpSysMat(k)/SpColMax(j) 
    enddo
    
    return
    end       
!
!
!
    subroutine SetupODEBoundaryConditions2D(ij)
    use dbase
    use GlobalVariables
    implicit none
    integer i, j, ij, k, l, jmem, morder, n1order, n2order, nsize, n1, n2, m, idiff
    integer icons, ibc, ic, ij_rel, jflag
    real *8 coef, gfunval, xsi1, xsi2, delta_mag, rcond, max_rcond, vol, darea
    logical is_found
    integer, allocatable :: mask(:)
    integer nmask
    real *8 emod, pp, nu, ii, dd, ll
!    
    
    emod = 30.0d9
    nu = 0.3d0
    pp = -1000.0d0
    ll = 8.0d0
    dd = 3.0d0
    ii = dd**3 / 12.0d0
    
    allocate(mask(2*totnode))
    n1order = PDoperator%n1order
    n2order = PDoperator%n2order
    nsize = PDoperator%nsize
    morder = PDoperator%morder
!
    icons = 0
    ij_rel = ij
    do ibc = 1 , PDoperator%num_bc
       if(PDoperator%bc(ibc)%u_flag > 0 ) then
          do k = 1 , totnode
             delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             !SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val

             SysVec(2*totnode+icons) = pp*coord(k,2)/(6.0d0 * emod * ii) * ( (6.0d0 * ll - 3.0d0 * coord(k,1)) * coord(k,1) + (2.0d0 + nu) * (coord(k,2)**2- dd**2 / 4.0d0))

             
             
             call FormDiffAmat2D(n1order, n2order, nsize, k)
             call recondition(DiffAmat2D,rcvec,nsize)
             call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)             
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_u
                n1 = PDoperator%bc(ibc)%uDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%uDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%uDiffs(idiff)%coef
                if( n1 > 0 .or. n2 > 0 ) then
                   call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                   do i = 1 , nsize
                      DiffAvec(i) = 0.0d0
                      do j = 1 , nsize
                         DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                      enddo
                   enddo
                   call recover(DiffAvec,rcvec,nsize)
                   do jmem = 2 , numfam(k)
                      j = family(k)%node(jmem)
                      jflag = family(k)%nodeFlag(jmem)
                      xsi1 = coord(j,1) - coord(k,1)
                      xsi2 = coord(j,2) - coord(k,2)
                      call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                      call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                      gfunval = 0.0d0
                      do l = 1 , nsize 
                         gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                      enddo
!
                      SpSysMat(ij_rel + 2*jmem - 1) = SpSysMat(ij_rel + 2*jmem - 1) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
                      SpSysMat(ij_rel + 1         ) = SpSysMat(ij_rel + 1)          - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
!
                      SpSysMat(ij_rel + 2*jmem) = SpSysMat(ij_rel + 2*jmem-1)
                      SpSysMat(ij_rel + 2     ) = SpSysMat(ij_rel + 1)
                   enddo
                else
                   SpSysMat(ij_rel + 1) = coef
                   SpSysMat(ij_rel + 2) = coef
                endif
             enddo
             ij_rel = ij_rel + 2*numfam(k)
          enddo
       elseif( PDoperator%bc(ibc)%u_flag < 0) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val
             SpSysMat(ij_rel + 1) = coef
             SpSysMat(ij_rel + 2) =-coef
             SpSysMat(ij_rel + 3) = coef
             SpSysMat(ij_rel + 4) =-coef
             ij_rel = ij_rel + 4
          enddo
       elseif(PDoperator%bc(ibc)%v_flag > 0 ) then
          do k = 1 , totnode
             delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             !SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val
             
              SysVec(2*totnode+icons) = -pp/(6.0d0 * emod * ii) * (3.0d0 * nu * coord(k,2)**2 * (ll - coord(k,1)) + (4.0d0 + 5.0d0 * nu) * dd**2 * coord(k,1)/4.0d0 + (3.0d0*ll-coord(k,1)) * coord(k,1)**2 )
             
             
             call FormDiffAmat2D(n1order, n2order, nsize, k)
             call recondition(DiffAmat2D,rcvec,nsize)
             call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)             
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_v
                n1 = PDoperator%bc(ibc)%vDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%vDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%vDiffs(idiff)%coef
                if( n1 > 0 .or. n2 > 0 ) then
                   call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                   do i = 1 , nsize
                      DiffAvec(i) = 0.0d0
                      do j = 1 , nsize
                         DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                      enddo
                   enddo
                   call recover(DiffAvec,rcvec,nsize)
                   do jmem = 2 , numfam(k)
                      j = family(k)%node(jmem)
                      jflag = family(k)%nodeFlag(jmem)
                      xsi1 = coord(j,1) - coord(k,1)
                      xsi2 = coord(j,2) - coord(k,2)
                      call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                      call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                      gfunval = 0.0d0
                      do l = 1 , nsize 
                         gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                      enddo
!
                      SpSysMat(ij_rel + 2*jmem - 1) = SpSysMat(ij_rel + 2*jmem - 1) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
                      SpSysMat(ij_rel + 1)          = SpSysMat(ij_rel + 1)          - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
!                      
                      SpSysMat(ij_rel + 2*jmem) = SpSysMat(ij_rel + 2*jmem-1)
                      SpSysMat(ij_rel + 2)      = SpSysMat(ij_rel + 1)
                   enddo
                else
                   SpSysMat(ij_rel + 1) = coef
                   SpSysMat(ij_rel + 2) = coef
                endif
             enddo
             ij_rel = ij_rel + 2*numfam(k)
          enddo
       elseif( PDoperator%bc(ibc)%v_flag < 0) then
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val
             SpSysMat(ij_rel + 1) = coef
             SpSysMat(ij_rel + 2) =-coef
             SpSysMat(ij_rel + 3) = coef
             SpSysMat(ij_rel + 4) =-coef
             ij_rel = ij_rel + 4
          enddo          
       elseif(PDoperator%bc(ibc)%sxx_flag > 0 .or. PDoperator%bc(ibc)%syy_flag > 0 .or. PDoperator%bc(ibc)%sxy_flag > 0) then
          do k = 1 , totnode
             delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             icons = icons + 1
             
             if (ibc==4) then
             SysVec(2*totnode+icons) = -pp/(2.0d0 * ii) * (dd**2 / 4.0d0 - coord(k,2)**2)
             else
             SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val
             endif
             
             call FormDiffAmat2D(n1order, n2order, nsize, k)
             call recondition(DiffAmat2D,rcvec,nsize)
             call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)             
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_u
                n1 = PDoperator%bc(ibc)%uDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%uDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%uDiffs(idiff)%coef
                call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                do i = 1 , nsize
                   DiffAvec(i) = 0.0d0
                   do j = 1 , nsize
                      DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                   enddo
                enddo
                call recover(DiffAvec,rcvec,nsize)
                do jmem = 2 , numfam(k)
                   j = family(k)%node(jmem)
                   jflag = family(k)%nodeFlag(jmem)
                   xsi1 = coord(j,1) - coord(k,1)
                   xsi2 = coord(j,2) - coord(k,2)
                   call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                   call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                   gfunval = 0.0d0
                   do l = 1 , nsize 
                      gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                   enddo
!
                   SpSysMat(ij_rel + 2*jmem - 1) = SpSysMat(ij_rel + 2*jmem - 1) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
                   SpSysMat(ij_rel + 1         ) = SpSysMat(ij_rel + 1)          - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
!
                   SpSysMat(ij_rel + 2*jmem) = SpSysMat(ij_rel + 2*jmem-1)
                   SpSysMat(ij_rel + 2     ) = SpSysMat(ij_rel + 1)
                enddo
             enddo
             ij_rel = ij_rel + 2*numfam(k)
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_v
                n1 = PDoperator%bc(ibc)%vDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%vDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%vDiffs(idiff)%coef
                call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                do i = 1 , nsize
                   DiffAvec(i) = 0.0d0
                   do j = 1 , nsize
                      DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                   enddo
                enddo
                call recover(DiffAvec,rcvec,nsize)
                do jmem = 2 , numfam(k)
                   j = family(k)%node(jmem)
                   jflag = family(k)%nodeFlag(jmem)
                   xsi1 = coord(j,1) - coord(k,1)
                   xsi2 = coord(j,2) - coord(k,2)
                   call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                   call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                   gfunval = 0.0d0
                   do l = 1 , nsize 
                      gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                   enddo
!
                   SpSysMat(ij_rel + 2*jmem - 1) = SpSysMat(ij_rel + 2*jmem - 1) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
                   SpSysMat(ij_rel + 1)          = SpSysMat(ij_rel + 1)          - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))
!                      
                   SpSysMat(ij_rel + 2*jmem) = SpSysMat(ij_rel + 2*jmem-1)
                   SpSysMat(ij_rel + 2)      = SpSysMat(ij_rel + 1)
                enddo
             enddo
             ij_rel = ij_rel + 2*numfam(k)
          enddo
       elseif(PDoperator%bc(ibc)%sxx_flag < 0 .or. PDoperator%bc(ibc)%syy_flag <0 .or. PDoperator%bc(ibc)%sxy_flag <0 ) then
          icons = icons + 1
          SysVec(2*totnode+icons) = PDoperator%bc(ibc)%val
          do k = 1 , 2*totnode
             mask(k) = 0
          enddo          
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             do jmem = 1 , numfam(k)
                j = family(k).node(jmem)
                mask(2*j-1) = 1
                mask(2*j) = 1
             enddo
          enddo
          nmask = 0
          do k = 1 , 2*totnode
             if(mask(k) > 0 ) then
                nmask = nmask + 1
                mask(k) = nmask
             endif
          enddo
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             darea = 1.0d0
             delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
             call FormDiffAmat2D(n1order, n2order, nsize, k)
             call recondition(DiffAmat2D,rcvec,nsize)
             call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)             
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_u
                n1 = PDoperator%bc(ibc)%uDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%uDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%uDiffs(idiff)%coef
                call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                do i = 1 , nsize
                   DiffAvec(i) = 0.0d0
                   do j = 1 , nsize
                      DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                   enddo
                enddo
                call recover(DiffAvec,rcvec,nsize)
                do jmem = 2 , numfam(k)
                   j = family(k)%node(jmem)
                   jflag = family(k)%nodeFlag(jmem)
                   xsi1 = coord(j,1) - coord(k,1)
                   xsi2 = coord(j,2) - coord(k,2)
                   call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                   call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                   gfunval = 0.0d0
                   do l = 1 , nsize 
                      gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                   enddo
!
                   SpSysMat(ij_rel + mask(2*j-1)) = SpSysMat(ij_rel + mask(2*j-1)) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))*darea
                   SpSysMat(ij_rel + mask(2*k-1)) = SpSysMat(ij_rel + mask(2*k-1)) - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))*darea
!
                   SpSysMat(ij_rel + nmask + mask(2*j-1)) = SpSysMat(ij_rel + mask(2*j-1))
                   SpSysMat(ij_rel + nmask + mask(2*k-1)) = SpSysMat(ij_rel + mask(2*k-1))
                enddo
             enddo
             do idiff = 1 , PDoperator%bc(ibc)%num_diff_ops_v
                n1 = PDoperator%bc(ibc)%vDiffs(idiff)%n1
                n2 = PDoperator%bc(ibc)%vDiffs(idiff)%n2
                coef = PDoperator%bc(ibc)%vDiffs(idiff)%coef
                call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
                do i = 1 , nsize
                   DiffAvec(i) = 0.0d0
                   do j = 1 , nsize
                      DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
                   enddo
                enddo
                call recover(DiffAvec,rcvec,nsize)
                do jmem = 2 , numfam(k)
                   j = family(k)%node(jmem)
                   jflag = family(k)%nodeFlag(jmem)
                   xsi1 = coord(j,1) - coord(k,1)
                   xsi2 = coord(j,2) - coord(k,2)
                   call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
                   call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
                   gfunval = 0.0d0
                   do l = 1 , nsize 
                      gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
                   enddo
!
                   SpSysMat(ij_rel + mask(2*j)) = SpSysMat(ij_rel + mask(2*j)) + coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))*darea
                   SpSysMat(ij_rel + mask(2*k)) = SpSysMat(ij_rel + mask(2*k)) - coef*gfunval*dvolume(j)*jflag/(delta_mag**(n1+n2))*darea
!                      
                   SpSysMat(ij_rel + nmask + mask(2*j)) = SpSysMat(ij_rel + mask(2*j))
                   SpSysMat(ij_rel + nmask + mask(2*k)) = SpSysMat(ij_rel + mask(2*k))
                enddo
             enddo
          enddo
          ij_rel = ij_rel + 2*nmask
       elseif( PDoperator%bc(ibc)%bx_flag > 0) then
          vol = 0
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             vol = vol + dvolume(k)
          enddo          
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             SysVec(2*k-1) =-PDoperator%bc(ibc)%val/vol
          enddo          
       elseif( PDoperator%bc(ibc)%by_flag > 0) then
          vol = 0
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             vol = vol + dvolume(k)
          enddo          
          do k = 1 , totnode
             if( coord(k,1) < PDoperator%bc(ibc)%xmin .or. coord(k,1) > PDoperator%bc(ibc)%xmax ) cycle
             if( coord(k,2) < PDoperator%bc(ibc)%ymin .or. coord(k,2) > PDoperator%bc(ibc)%ymax ) cycle
             SysVec(2*k) =-PDoperator%bc(ibc)%val/vol
          enddo          
       endif
    enddo
!
    deallocate(mask)
!    
    return
    end
!
!
!
!
    subroutine CalcStressStrain()
    use dbase
    implicit none
    real *8 uval_k, vval_k, uval_j, vval_j
    integer k, j, i, l, jflag
    integer n1, n2, nsize, n1order, n2order, m, jmem, morder
    real *8 rcond, delta_mag, dudx, dudy, dvdx, dvdy
    real *8 epsxx, epsyy, gamxy, sigxx, sigyy, sigxy, xsi1, xsi2, gfunval
    real *8 emod, pp, nu, ii, dd, ll, exactu, exactv, exactsxx, exactsxy, errx
!    
    
    emod = 30.0d9
    nu = 0.3d0
    pp = -1000.0d0
    ll = 8.0d0
    dd = 3.0d0
    ii = dd**3 / 12.0d0
!    
    if( PDoperator%atype == 0 ) then
        allocate(DiffAmat2D(PDoperator%nsize,PDoperator%nsize))
        allocate(DiffAmatInv2D(PDoperator%nsize,PDoperator%nsize))
        allocate(DiffBvec2D(PDoperator%nsize))
        allocate(DiffAvec(PDoperator%nsize))
        allocate(plist(PDoperator%nsize))
        allocate(blist(PDoperator%nsize))
        allocate(weight(PDoperator%nsize))
        allocate(rcvec(PDoperator%nsize))
    endif
    n1order = PDoperator%n1order
    n2order = PDoperator%n2order
    nsize = PDoperator%nsize
    morder = 2
    open(28,file='./results/derivatives.dat')
    do k = 1 , totnode
       uval_k = fvec(2*k-1)
       vval_k = fvec(2*k)
       call FormDiffAmat2D(n1order, n2order, nsize, k)
       call recondition(DiffAmat2D,rcvec,nsize)
       call inverse2(nsize, DiffAmat2D, DiffAmatInv2D, k, rcond)
       delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
          n1 = 1
          n2 = 0          
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          dudx = 0.0d0
          dvdx = 0.0d0
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             uval_j = fvec(2*j-1)
             vval_j = fvec(2*j)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             dudx = dudx + uval_j*gfunval*dvolume(j)*jflag
             dudx = dudx - uval_k*gfunval*dvolume(j)*jflag
             dvdx = dvdx + vval_j*gfunval*dvolume(j)*jflag
             dvdx = dvdx - vval_k*gfunval*dvolume(j)*jflag
          enddo
          dudx = (1.0d0/delta_mag**(n1+n2))*dudx
          dvdx = (1.0d0/delta_mag**(n1+n2))*dvdx
!
          n1 = 0
          n2 = 1          
          call FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
          do i = 1 , nsize
             DiffAvec(i) = 0.0d0
             do j = 1 , nsize
                DiffAvec(i) = DiffAvec(i) + DiffAmatInv2D(i,j)*DiffBvec2D(j)
             enddo
          enddo
          call recover(DiffAvec,rcvec,nsize)
          dudy = 0.0d0
          dvdy = 0.0d0
          do jmem = 2 , numfam(k)
             j = family(k)%node(jmem)
             jflag = family(k)%nodeFlag(jmem)
             uval_j = fvec(2*j-1)
             vval_j = fvec(2*j)
             xsi1 = coord(j,1) - coord(k,1)
             xsi2 = coord(j,2) - coord(k,2)
             call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag ) 
             call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
             gfunval = 0.0d0
             do l = 1 , nsize 
                gfunval = gfunval + DiffAvec(l)*plist(l)*weight(l)
             enddo
             dudy = dudy + uval_j*gfunval*dvolume(j)*jflag
             dudy = dudy - uval_k*gfunval*dvolume(j)*jflag
             dvdy = dvdy + vval_j*gfunval*dvolume(j)*jflag
             dvdy = dvdy - vval_k*gfunval*dvolume(j)*jflag
          enddo
          dudy = (1.0d0/delta_mag**(n1+n2))*dudy
          dvdy = (1.0d0/delta_mag**(n1+n2))*dvdy
          epsxx = dudx
          epsyy = dvdy
          gamxy = dudy + dvdx
          sigxx = PDoperator%C11*epsxx + PDoperator%C12*epsyy 
          sigyy = PDoperator%C12*epsxx + PDoperator%C22*epsyy 
          sigxy = PDoperator%C66*gamxy
          
          exactu   = pp*coord(k,2)/(6.0d0 * emod * ii) * ( (6.0d0 * ll - 3.0d0 * coord(k,1)) * coord(k,1) + (2.0d0 + nu) * (coord(k,2)**2- dd**2 / 4.0d0))
          exactv   =  -pp/(6.0d0 * emod * ii) * (3.0d0 * nu * coord(k,2)**2 * (ll - coord(k,1)) + (4.0d0 + 5.0d0 * nu) * dd**2 * coord(k,1)/4.0d0 + (3.0d0*ll-coord(k,1)) * coord(k,1)**2 )
          exactsxx =  pp * (ll - coord(k,1)) * coord(k,2) / ii
          exactsxy =  -pp/(2.0d0 * ii) * (dd**2 / 4.0d0 - coord(k,2)**2)
          
          errx = uval_k-exactu

!      write(28,"(11(e16.9,2x))") (coord(k,j), j=1,2), uval_k, vval_k, epsxx, epsyy, gamxy, sigxx, sigyy, sigxy
       write(28,"(11(e16.9,2x))") (coord(k,j), j=1,2), uval_k, vval_k, sigxx, sigxy, exactu, exactv, exactsxx, exactsxy
    enddo
    close(28)
    deallocate(DiffAmat2D)
    deallocate(DiffBvec2D)
    deallocate(DiffAvec)
    deallocate(plist)
    deallocate(blist)
    deallocate(weight)
    return
    end

!
!
!
    subroutine FormDiffBvec2D(n1order, n2order, nsize, n1, n2, m )
    use dbase
    implicit none
    integer morder, i, nsize, n1order, n2order, n1, n2, m
!    
    morder = 2
    call b_operator_2d( n1order, n2order, nsize, n1, n2, m )    
    do i = 1 , nsize
       DiffBvec2D(i) = blist(i)
    enddo
    return
    end
!
!
!    
    subroutine FormDiffAmat2D(n1order, n2order, nsize, k )
    use dbase
    use GlobalVariables
    implicit none
    integer morder, n1order, n2order, nsize, k, jmem, j, ii, jj, jflag
    real *8 xsi1, xsi2, delta_mag, tol
!    
    morder = 2
    delta_mag = dsqrt(deltax(k)*deltax(k)+deltay(k)*deltay(k))
!
    do ii = 1 , nsize
       do jj = 1 , nsize
          DiffAmat2D(ii,jj) = 0.0d0
       enddo
    enddo
    
    do jmem = 1 , numfam(k)
       j = family(k)%node(jmem)
       jflag = family(k)%nodeFlag(jmem)
       xsi1 = coord(j,1) - coord(k,1)
       xsi2 = coord(j,2) - coord(k,2)
       call p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
       call weights_2d( n1order, n2order, nsize, xsi1, xsi2, dsqrt(deltax(k)**2+deltay(k)**2) )
       do ii = 1 , nsize 
          do jj = 1 , nsize
             DiffAmat2D(ii,jj) = DiffAmat2D(ii,jj) + weight(jj)*plist(ii)*plist(jj)*dvolume(j)*jflag
          enddo
       enddo
    enddo
    return
    end
!
!
!
    subroutine getSize2D(n1order, n2order, nsize)
    implicit none
    integer n1order, n2order, nsize, iterm
!        
        iterm = 0
        !if( n1order >=0 .and. n2order >= 0 ) then
        !    iterm = iterm + 1
        !endif
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
            iterm = iterm + 1
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
        endif
        nsize = iterm
    return
    end        
!
!
!
    subroutine SetDiffOperators2D(n1order, n2order, nsize)
    use dbase
    implicit none
    integer n1order, n2order, nsize, iterm
!        
    allocate(PDoperator.order(nsize))
        iterm = 0
        !if( n1order >=0 .and. n2order >= 0 ) then
        !    iterm = iterm + 1
        !    PDoperator.order(iterm).n1=0
        !    PDoperator.order(iterm).n2=0
        !    PDoperator.order(iterm).n3=0
        !endif
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
           PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=5
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=5
            PDoperator.order(iterm).n3=0
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=6
            PDoperator.order(iterm).n2=0
            PDoperator.order(iterm).n3=0
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=5
            PDoperator.order(iterm).n2=1
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=4
            PDoperator.order(iterm).n2=2
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=3
            PDoperator.order(iterm).n2=3
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=2
            PDoperator.order(iterm).n2=4
            PDoperator.order(iterm).n3=0
            iterm = iterm + 1
            PDoperator.order(iterm).n1=1
            PDoperator.order(iterm).n2=5
            PDoperator.order(iterm).n3=0
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
            PDoperator.order(iterm).n1=0
            PDoperator.order(iterm).n2=6
            PDoperator.order(iterm).n3=0
        endif
    return
    end        
!
!
!
    subroutine p_operator_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
    use dbase
    implicit none
    integer morder, n1order, n2order, nsize, iterm
    real *8 xsi1, xsi2, xsi1p, xsi2p, delta_mag
!    
    xsi1p = xsi1/delta_mag
    xsi2p = xsi2/delta_mag
    iterm = 0
    !if( n1order >= 0 .and. n2order >= 0 ) then
    !    iterm = iterm + 1
    !    plist(iterm) = 1
    !endif
    if( n1order >= 1 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p
    endif
    if( n2order >= 1 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p
    endif
!    
    if( n1order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p
    endif
    if( n1order >= 2 .and. n2order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi2p
    endif
    if( n2order >= 2 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi2p*xsi2p
    endif
!    
    if( n1order >= 3 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 3 .and. n2order >= 3 ) then
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm)  = xsi1p*xsi2p*xsi2p
    endif
    if( n2order >= 3 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p*xsi2p*xsi2p
    endif
!
    if( n1order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 4 .and. n2order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p
    endif
    if( n2order >= 4 ) then
        iterm = iterm + 1
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p
    endif
!
    if( n1order >= 5 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p
    endif
    if( n1order >= 5 .and. n2order >= 5 ) then
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
    endif
    if( n2order >= 5 ) then
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
    endif
!    
    if( n1order >= 6 ) then
        iterm = iterm + 1
    endif
    if( n1order >= 6 .and. n2order >= 6 ) then
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p*xsi1p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi1p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi1p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi1p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi1p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
        plist(iterm) = xsi1p*xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
        iterm = iterm + 1
    endif
    if( n1order >= 6 ) then
        plist(iterm) = xsi2p*xsi2p*xsi2p*xsi2p*xsi2p*xsi2p
    endif
!    
    return
    end
!
!
!
    subroutine weights_2d( n1order, n2order, nsize, xsi1, xsi2, delta_mag )
    use dbase
    implicit none
    integer morder, n1order, n2order, nsize, iterm
    real *8 xsi1, xsi2, xsi_mag, wt, delta_mag
!    
    morder = 2
    xsi_mag = dsqrt(xsi1*xsi1+xsi2*xsi2)
!
    wt = exp(-4*(xsi_mag/delta_mag)**2)
!
!    
    iterm = 0

!
    if( n1order >= 1 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 1 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 2 ) then
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
    if( n1order >= 2 .and. n2order >= 2 ) then
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
    if( n2order >= 2 ) then
        iterm = iterm + 1
        weight(iterm)  = wt
    endif
!
    if( n1order >= 3 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 3 .and. n2order >= 3 ) then
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 3 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 4 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 4 .and. n2order >= 4) then
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 4 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!
    if( n1order >= 5 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 5 .and. n2order >= 5 ) then
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 5 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    if( n1order >= 6 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n1order >= 6 .and. n2order >= 6 ) then
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
        iterm = iterm + 1
        weight(iterm) = wt
    endif
    if( n2order >= 6 ) then
        iterm = iterm + 1
        weight(iterm) = wt
    endif
!    
    return
    end
!
!
!    
    subroutine b_operator_2d( n1order, n2order, nsize, n1, n2, m )
    use dbase
    implicit none
    integer n1, n2, m, morder, n1order, n2order, nsize, iterm, i
    real *8 fn1, fn2, coef
!    
    morder = 2
    if( n1 > n1order .or. n2 > n2order ) stop
    do i = 1 , nsize
       blist(i) = 0.0d0
    enddo
    fn1=1
    fn2=1
    do i = 1 , n1
       fn1 = fn1*i
    enddo
    do i = 1 , n2
       fn2 = fn2*i
    enddo
    coef = fn1*fn2
        iterm = 0
!        
        if( n1order >= 1 ) then
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n2order >= 1 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 1 ) then 
                m = iterm
            endif
        endif
!        
        if( n1order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 2 .and. n2order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 1 ) then
                m = iterm
            endif
        endif
        if( n2order >= 2 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 2 ) then
                m = iterm
            endif
        endif
!        
        if( n1order >= 3 ) then
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 3 .and. n2order >= 3) then
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 2 ) then
                m = iterm
            endif
        endif
        if( n2order >= 3 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 3 ) then
                m = iterm
            endif
        endif
!
        if( n1order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 4 .and. n2order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 3 ) then
                m = iterm
            endif
        endif
        if( n2order >= 4 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 4 ) then
                m = 15
            endif
        endif
!
        if( n1order >= 5 ) then
            iterm = iterm + 1
            if( n1 == 5 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 5 .and. n1order >= 5) then
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 3 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 4 ) then
                m = iterm
            endif
        endif
        if( n2order >= 5 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 5 ) then
                m = iterm
            endif
        endif
!
        if( n1order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 6 .and. n2 == 0 ) then
                m = iterm
            endif
        endif
        if( n1order >= 6 .and. n2order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 5 .and. n2 == 1 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 4 .and. n2 == 2 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 3 .and. n2 == 3 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 2 .and. n2 == 4 ) then
                m = iterm
            endif
            iterm = iterm + 1
            if( n1 == 1 .and. n2 == 5 ) then
                m = iterm
            endif
        endif
        if( n2order >= 6 ) then
            iterm = iterm + 1
            if( n1 == 0 .and. n2 == 6 ) then
               m = iterm
            endif
        endif
    blist(m) = coef                                       
    return                                      
    end
!
!
!
!
!////////////////////////////////UTILITY SUBROUTINES /////////////////////////
!
!
!
      subroutine ApplySlit()
      use dbase
      implicit none
      integer k, j, islit, jmem
      real *8 xavrg, yavrg
!      
      do islit = 1 , PDoperator%nslits
         if(dabs(slits(islit)%x2-slits(islit)%x1)>dabs(slits(islit)%y2-slits(islit)%y1)) then
            yavrg = (slits(islit)%y2+slits(islit)%y1)/2.0d0
            do k = 1 , totnode
               if( (coord(k,1) < slits(islit)%x1) .or. (coord(k,1) > slits(islit)%x2) ) cycle
               do jmem = 2 , numfam(k)
                  j = family(k)%node(jmem)
                  if( (coord(j,1) < slits(islit)%x1) .and. (coord(j,1) > slits(islit)%x2) ) cycle
                  if( (coord(k,2)-yavrg)*(coord(j,2)-yavrg) < 0.0d0 ) then
                      family(k)%nodeFlag(jmem) = 0
                  endif
               enddo
            enddo
         else
            xavrg = (slits(islit)%x2+slits(islit)%x1)/2.0d0
            do k = 1 , totnode
               if( (coord(k,2) < slits(islit)%y1) .or. (coord(k,2) > slits(islit)%y2) ) cycle
               do jmem = 2 , numfam(k)
                  j = family(k)%node(jmem)
                  if( (coord(j,2) < slits(islit)%y1) .and. (coord(j,2) > slits(islit)%y2) ) cycle
                  if( (coord(k,1)-xavrg)*(coord(j,1)-xavrg) < 0.0d0 ) then
                      family(k)%nodeFlag(jmem) = 0
                  endif
               enddo
            enddo
         endif
      enddo
      return
      end         
!
!
!
      subroutine inverse( n, a, y, np )
      implicit none
      integer nmax
      parameter (nmax=5000)
      integer np, indx(NMAX), n, i, j
      real *8 a(np,np), y(np,np), d
      do i = 1 , n
         do j = 1 , n
            y(i,j) = 0.
         enddo
         y(i,i) = 1.
      enddo
      call ludcmp(a,n,np,indx,d)
      do j = 1 , n
         call lubksb(a,n,np,indx,y(1,j))
      enddo
      return
      end subroutine inverse
!
!
!
      subroutine solve( n, a, y, np )      
      implicit none
      integer nmax
      parameter (nmax=5000)
      integer np, indx(NMAX), n, i, j
      real *8 a(np,np), y(np), d
      call ludcmp(a,n,np,indx,d)
      call lubksb(a,n,np,indx,y)
      return
      end subroutine solve
!
!
!
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none
      integer nmax
      parameter(nmax=5000)
      INTEGER n,np,indx(n)
      REAL *8 d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL *8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.) then
            print *, 'singular matrix in ludcmp'
            read(*,*)
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE ludcmp
!
!
!
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      REAL *8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL *8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END SUBROUTINE lubksb
!
!
!
      subroutine recondition( amat, rcvec, n )
      implicit none
      integer n, i, j
      real *8 amat(n,n), rcvec(n), cmax
      do j = 1 , n
         cmax = dabs(amat(1,j))
         do i = 1 , n
            if(dabs(amat(i,j)) > cmax) cmax = dabs(amat(i,j))
         enddo
         rcvec(j) = cmax
         do i = 1 , n
            amat(i,j) = amat(i,j)/cmax
         enddo
      enddo
!         
      return
      end
!
!
!
      subroutine recover( avec, rcvec, n )
      implicit none
      integer n, i
      real *8 avec(n), rcvec(n)
      do i = 1 , n
         avec(i) = avec(i)/rcvec(i)
      enddo
!         
      return
    end
!
    recursive subroutine quicksort(a, order, first, last)
  implicit none
  real *8  a(*), x, t
  integer first, last
  integer *8 order(*)
  integer *8 i, j, itemp

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     itemp = order(i); order(i) = order(j); order(j) = itemp
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, order , first, i-1)
  if (j+1 < last)  call quicksort(a, order , j+1, last)
end subroutine quicksort        
