 !
!   Peridynamics
!
module dbase
!
implicit none
!
integer totnode
real *8, allocatable :: coord(:,:)
integer morder
real *8 length
real *8 width
real *8 height
real *8, allocatable :: c0vec(:), c1vec(:), c2vec(:), c3vec(:), c4vec(:), c5vec(:)
real *8, allocatable :: fvec(:), dvolume(:), deltax(:), deltay(:), deltaz(:)
integer ndivx,ndivy,ndivz
real *8 dx, dy, dz, deltaxk, deltayk, deltazk, xmin, xmax, ymin, ymax, zmin, zmax
real *8, parameter :: pi=dacos(-1.0d0)

end module dbase
!
!
program main
use dbase
implicit none
!
call input()
!
call ProblemDimensions()
!
call AllocateAndInitializeArrays()
!
call GridGeneration()
if(morder == 1 ) then
   call getFunVal1D()
elseif(morder == 2 ) then
   call getFunVal2D()
elseif(morder == 3 ) then
   call getFunVal3D()
endif
!
end program main
!
!
!
      subroutine input()
      use dbase
      implicit none
      character line
!
      open(1,file='grid.inp')
!
      read(1,'(A)') line
      read(1,*) morder
      if(morder==1) then
         read(1,'(A)') line
         read(1,*) xmin, xmax
         length = xmax - xmin
         ymin = 0.0d0
         ymax = 0.0d0
         zmin = 0.0d0
         zmax = 0.0d0
         width = 0.0d0
         height = 0.0d0
         read(1,'(A)') line
         read(1,*) dx, deltaxk
         dy = 1.0d0
         dz = 1.0d0
         deltayk = 0.0d0
         deltazk = 0.0d0
      elseif(morder==2) then
         read(1,'(A)') line
         read(1,*) xmin, xmax, ymin, ymax
         length = xmax - xmin
         width = ymax - ymin
         zmin = 0.0d0
         zmax = 0.0d0
         height = 0.0d0
         read(1,'(A)') line
         read(1,*) dx, dy, deltaxk, deltayk
         dz = 1.0d0
         deltazk = 0.0d0
      elseif(morder==3) then
         read(1,'(A)') line
         read(1,*) xmin, xmax, ymin, ymax, zmin, zmax
         length = xmax - xmin
         width = ymax - ymin
         height = zmax - zmin
         read(1,'(A)') line
         read(1,*) dx, dy, dz, deltaxk, deltayk, deltazk
      endif
!
      close(1)
!
      return
!
      end subroutine input
!
!
!
      subroutine getFunVal1D()
      use dbase
      implicit none
      integer i
      real *8 epsilon
!
      allocate(fvec(totnode))
      allocate(c0vec(totnode))
      allocate(c1vec(totnode))
      allocate(c2vec(totnode))
!
      open(230,file='PDgeom1D.dat')
      write(230,"(i6)") totnode
      epsilon = 0.0001d0
      do i = 1 , totnode
         c2vec(i) = coord(i,1)
         c1vec(i) = 2.0d0
         c0vec(i) = coord(i,1)
         fvec(i) =0.0d0
         write(230,"(11(2x,e16.9))") coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), c2vec(i), c1vec(i), c0vec(i), fvec(i)
      enddo
      close(230)
      return
      end
!
!
!
      subroutine getFunVal1Da()
      use dbase
      implicit none
      integer i
      real *8 epsilon
!
      allocate(fvec(totnode))
      allocate(c0vec(totnode))
      allocate(c1vec(totnode))
      allocate(c2vec(totnode))
!
      open(230,file='PDgeom1D.dat')
      write(230,"(i6)") totnode
      do i = 1 , totnode
         c2vec(i) = 1.0d0
         c1vec(i) = 2.0d0*coord(i,1)
         c0vec(i) =-4.0d0
         fvec(i) = 4.0d0*(derf(coord(i,1))**2-2.0d0/dsqrt(pi)*coord(i,1)*dexp(-coord(i,1)**2)*derf(coord(i,1))-2.0d0/pi*dexp(-2.0d0*coord(i,1)**2)+2.0d0/pi*dexp(-coord(i,1)**2)-1.0d0)
         write(230,"(12(2x,e16.9))") coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), c2vec(i), c1vec(i), c0vec(i), fvec(i)
      enddo
      close(230)
      return
      end
!
!
!
      subroutine getFunVal2Da()
      use dbase
      implicit none
      integer i
      real *8, allocatable :: c20vec(:), c02vec(:), c00vec(:), c11vec(:), c10vec(:), c01vec(:)
      real *8 rho, c_sh, k_tc
!      
      allocate(c20vec(totnode))      
      allocate(c11vec(totnode))      
      allocate(c02vec(totnode))      
      allocate(c01vec(totnode))     
      allocate(c10vec(totnode))     
      allocate(c00vec(totnode))     
      allocate(fvec(totnode))
!
      k_tc = 50.0d0
      c_sh = 500.0d0
      rho = 7850.0d0
!
      open(230,file='PDgeom2D.dat')
      write(230,"(i6)") totnode
      do i = 1 , totnode
         c20vec(i) = 1.0d0
         c01vec(i) =-rho*c_sh/600.0d0/25.0/k_tc
         fvec(i) =0.0d0
         write(230,"(10(2x,e16.9))") coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), c20vec(i), c01vec(i), fvec(i)
      enddo
      close(230)
      return
      end
!
!
!
      subroutine getFunVal2D()
      use dbase
      implicit none
      integer i
      real *8, allocatable :: c20vec(:), c02vec(:), c00vec(:), c11vec(:), c10vec(:), c01vec(:), fvec_x(:), fvec_y(:)
      real *8 rho, c_sh, k_tc
!      
      allocate(c20vec(totnode))      
      allocate(c11vec(totnode))      
      allocate(c02vec(totnode))      
      allocate(c01vec(totnode))     
      allocate(c10vec(totnode))     
      allocate(c00vec(totnode))     
      allocate(fvec_x(totnode))
      allocate(fvec_y(totnode))
!
      k_tc = 50.0d0
      c_sh = 500.0d0
      rho = 7850.0d0
!
      open(230,file='PDgeom2D.dat')
      write(230,"(i6)") totnode
      do i = 1 , totnode
         fvec_x(i) =0.0d0
         fvec_y(i) =0.0d0
         write(230,"(10(2x,e16.9))") coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), fvec_x(i), fvec_y(i)
      enddo
      close(230)
      return
      end
!
!
!
      subroutine getFunVal3D()
      use dbase
      implicit none
      integer i
!
      open(230,file='PDgeom3D.dat')
      write(230,"(i6)") totnode
      do i = 1 , totnode
         fvec(i) = 3.0d0*dcos(pi*coord(i,1)/length)*dcos(pi*coord(i,2)/width)
         write(230,"(10(2x,e16.9))") coord(i,1), coord(i,2), coord(i,3), dvolume(i), deltax(i), deltay(i), deltaz(i), fvec(i)
      enddo
      close(230)
      return
      end
!
!
!
      subroutine ProblemDimensions()
      use dbase
      implicit none
!
      PRINT *, 'ProblemDimensions'
!
      ndivx = idnint(length/dx)
      if(ndivx==0) ndivx = 1
      ndivy = idnint(width/dy)
      if(ndivy==0) ndivy = 1
      ndivz = idnint(height/dz)
      if(ndivz==0) ndivz = 1
      totnode = totnode + ndivx*ndivy*ndivz
      print *, 'totnode =', totnode
      return
      end subroutine ProblemDimensions
!
!
!
      subroutine AllocateAndInitializeArrays()
      use dbase
      implicit none
!
      PRINT *,'AllocateAndInitializeArrays'
!
      allocate(coord(totnode,3))
      allocate(dvolume(totnode))
      allocate(deltax(totnode))
      allocate(deltay(totnode))
      allocate(deltaz(totnode))
!
      return
      end subroutine AllocateAndInitializeArrays
!
!
!
      subroutine GridGeneration()
      use dbase
      implicit none
      integer i, j, k
      integer nnum
!
      PRINT *,'GridGeneration'
!
!
      nnum = 0
      do k = 1 , ndivz
         do i = 1, ndivy
            do j = 1, ndivx
               nnum = nnum + 1
               coord(nnum,1) = xmin + (dx / 2.0d0) + (j-1) * dx
               coord(nnum,2) = ymin + (dy / 2.0d0) + (i-1) * dy
               coord(nnum,3) = zmin + (dz / 2.0d0) + (k-1) * dz
               dvolume(nnum) = dx*dy*dz
               deltax(nnum) = deltaxk
               deltay(nnum) = deltayk
               deltaz(nnum) = deltazk
             enddo
          enddo
      enddo
!
      return
      end subroutine GridGeneration
!
!
!
