 PROGRAM relax

  IMPLICIT NONE
  INTEGER  :: natms, ndim, i, iargc, neigv
  INTEGER  :: nthird
  INTEGER  :: idnode, mxnode
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: atm1,atm2,atm3,c1,c2

  REAL(8)                              :: eig, twopi, temp, domega
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dynmat, v3ijk
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: eigenvalues,eigvec,wgt
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: occupation, occupaclas
  LOGICAL :: linterac, ltube = .false.
  CHARACTER(LEN=100) :: buffer

  include "mpif.h"
 
  call initmpi(idnode,mxnode)

  if(iargc()<2) stop "Usage: lifetime.x <#atoms> <temperature> <delta(meV)>"
  call getarg(1,buffer)
  read(buffer,*) natms
  call getarg(2,buffer)
  read(buffer,*) temp
  call getarg(3,buffer)
  read(buffer,*) domega
  if (iargc()>3) ltube =.true.

  ndim = natms*3
 
  call readers

  call occupa

  call relaxation
  
  call exitmpi
  stop

 CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE occupa
 
 INTEGER i
 REAL(8), PARAMETER                   :: omegaK  = 1.4387752      !1/cm --> K
 
 ALLOCATE(occupation(ndim))
 ALLOCATE(occupaclas(ndim))

 !print *,'compute occupation'
 do i = 1, ndim
   occupation(i) = 1.d0/(exp(eigenvalues(i)*omegaK/temp)-1.d0)
   occupaclas(i) = temp/eigenvalues(i)/omegaK
 enddo
 return
 
 END SUBROUTINE occupa
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE relaxation
  
 INTEGER i, j, k, i3rd, ati, atj, atk, inph, ngarbage
 INTEGER, ALLOCATABLE, DIMENSION(:)   :: nintplus, nintminus
 REAL(8), PARAMETER     :: invcm2mev  = 1.239841875d-1 !1/cm --> meV
 REAL(8), PARAMETER     :: invcm2thz =  2.997924580d-2 !1/cm --> THz
 REAL(8), PARAMETER     :: hbar = 6.35075751 ! 10 j/mol ps
 REAL(8), PARAMETER     :: pi = 3.1415926535d0
 REAL(8), ALLOCATABLE, DIMENSION(:)   :: tau, eigvi, eigvj, eigvk, eigthz, vvec
 REAL(8), ALLOCATABLE, DIMENSION(:)   :: gamma2, gammaCL
 REAL(8), ALLOCATABLE, DIMENSION(:,:)   :: wrk3ijk
 REAL(8)  :: deltaa, deltab, occfactor, occ_classic, vijksq, vijksqCL
 CHARACTER(20) :: fname,frate

 INTERFACE 
   FUNCTION ddot(ndim,vvec,i,eigvk,j)
     REAL(8) ddot
     INTEGER, INTENT(IN) :: ndim, i, j
     REAL(8), INTENT(IN) :: vvec(ndim),eigvk(ndim)
   END FUNCTION ddot
 END INTERFACE 

 ngarbage = 0           ! don't use the last ngarbage phonons 
 domega = domega/4.135  ! meV to THz

 if(ltube) then 
   inph = 5
 else
   inph = 4
 endif

 ALLOCATE( eigvi(ndim),eigvj(ndim),eigvk(ndim),eigthz(ndim) )
 ALLOCATE( gamma2(ndim), gammaCL(ndim) )
 ALLOCATE( nintplus(ndim), nintminus(ndim) )

 eigthz = eigenvalues*invcm2thz

 if (idnode<10) then
   write(frate,'("Decay.00",i1)') idnode
 else if (idnode<100) then
   write(frate,'("Decay.0",i2)') idnode
 else if (idnode<1000) then
   write(frate,'("Decay.",i3)') idnode
 endif
 open(unit=44,file=frate)
 if(idnode==0) & 
   write(44,'("#a triangular delta, width (THz, meV): ",2f14.4)') domega, domega*4.135
 if(idnode==0) & 
   write(44,'("#b frequency (THz)   bandwidth QM (meV)   bandwidth CL (meV) 2->1  1->2")')
 call flush(44)
! write(fname,'("COUPLING",i1)') idnode
! open(unit=100,file=fname,form='formatted')

! rescale the force constants with atom masses

 do i3rd = 1, nthird
   v3ijk(:,i3rd) = v3ijk(:,i3rd)  & 
           /wgt(atm1(i3rd))/wgt(atm2(i3rd))/wgt(atm3(i3rd))
 enddo

! write(*,*) 'delta width (THz): ', domega

 ALLOCATE( wrk3ijk(ndim,ndim), vvec(ndim) )

 do i =  idnode+inph, ndim-ngarbage, mxnode  ! 1st loop over phonons
! do i =  idnode+inph, 7, mxnode  ! 1st loop over phonons
   print *, "phonon", i,"on node", idnode,mxnode
!  if( (eigthz(i)-eigthz(i-1))<1.d-9 ) then
!    gamma2(i) = gamma2(i-1)
!    gammaCL(i) = gammaCL(i-1)
!    write(44,*) eigthz(i), gamma2(i), gammaCL(i) 
!    cycle
!  endif
   nintplus(i)  = 0
   nintminus(i)  = 0
   eigvi = dynmat(:,i)
   gamma2(i) = 0.d0
   gammaCL(i) = 0.d0
   wrk3ijk = 0.d0

   do i3rd = 1, nthird
     ati = 3*(atm3(i3rd)-1)
     atj = 3*(atm2(i3rd)-1) + c2(i3rd)
     atk = 3*(atm1(i3rd)-1) + c1(i3rd)
     wrk3ijk(atj,atk) = wrk3ijk(atj,atk) + v3ijk(1,i3rd)*eigvi(ati+1) +      &
                                           v3ijk(2,i3rd)*eigvi(ati+2) +      &
                                           v3ijk(3,i3rd)*eigvi(ati+3) 
   enddo
   
   do j = inph, ndim-ngarbage  ! 2nd loop over phonons

!    Check whether there is interaction between i and j
     linterac=.false.
     do k = inph, ndim
       deltaa = abs(eigthz(i)-eigthz(j)-eigthz(k))
       deltab = abs(eigthz(i)+eigthz(j)-eigthz(k))
       if( deltaa<domega .or. deltab<domega ) then
         linterac=.true.
         exit
       endif
     enddo

     if(.not.linterac) cycle

     eigvj = dynmat(:,j)
!!!  vvec = MATMUL(eigvj,wrk3ijk)
     CALL DGEMM('N','N',ndim,1,ndim, 1.d0, wrk3ijk, ndim, eigvj, ndim, 0.d0, vvec, ndim )      
!!!  CALL DSYMM('L','U',ndim,1,1.d0,wrk3ijk,ndim,eigvj,ndim, 0.d0, vvec, ndim )

     do k = inph, ndim-ngarbage  ! 3rd loop over phonons
       
!      energy conservation

       deltaa = abs(eigthz(i)-eigthz(j)-eigthz(k))
       deltab = abs(eigthz(i)+eigthz(j)-eigthz(k))
       if (deltaa<domega) then
         deltaa = 1.d0/domega * (1.d0 - deltaa/domega) ! delta is approximated as a triangle
         !deltaa = 1.d0/domega  ! delta is approximated as a rectangle
         nintminus(i) = nintminus(i) + 1
       else 
         deltaa = 0.d0
       endif
       if (deltab<domega) then
         deltab = 1.d0/domega * (1.d0 - deltab/domega) ! delta is approximated as a triangle
         !deltab = 1.d0/domega
         nintplus(i) = nintplus(i) + 1
       else
         deltab = 0.d0
       endif
       if( (deltaa+deltab)<0.001 ) cycle
       eigvk = dynmat(:,k)
     

       vijksq = DOT_PRODUCT(vvec,eigvk)
!      vijksq = DDOT(ndim,vvec,1,eigvk,1)

       if (abs(vijksq)>1.d-6) &

!!     Write COUPLING files
!      write(100,'(3i6,4g16.8)') i, j, k, vijksq 

!      Calculate the phonon occupation factors (classical and QM)
       
       occfactor = 0.5d0*deltaa*(1.d0+occupation(j)+occupation(k)) +  &
                         deltab*(occupation(j)-occupation(k) )
!      occ_classic = 0.5d0*deltaa*(1.d0+occupaclas(j)+occupaclas(k)) +  &
!                        deltab*(occupaclas(j)-occupaclas(k) )
       occ_classic = (0.5d0*deltaa + deltab)*occupaclas(j)*occupaclas(k)
       vijksq = vijksq**2/eigthz(j)/eigthz(k)/16.d0/pi**4   ! 16 pi**4 converts frequencies in omega

       gamma2(i) = gamma2(i) + vijksq*occfactor
       gammaCL(i)= gammaCL(i) + vijksq*occ_classic

!!$ DEBUG      write(777,*) occ_classic, occfactor
     enddo   ! close 3rd phonon loop
   enddo     ! close 2nd phonon loop
   gamma2(i)  = hbar**2*pi*gamma2(i) /eigthz(i)/4.d0/9.648538 
   gammaCL(i) = hbar**2*pi*gammaCL(i)/eigthz(i)/4.d0/9.648538 / occupaclas(i) 
!   gammaCL(i) = hbar**2*pi*gammaCL(i)/eigthz(i)/4.d0/9.648538 
!     note: 1 meV = 9.648538 * 10 j/mol
   write(44,'(3g14.6,2x,2i8)') eigthz(i), gamma2(i), gammaCL(i), nintplus(i), nintminus(i)  ! sum(eigvi)
 enddo

 END SUBROUTINE relaxation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE readers

 INTEGER i, idx
 REAL(8) :: x,y,z

! read the dynamical matrix

  ALLOCATE(wgt(natms))
  ALLOCATE(dynmat(ndim,ndim))
  ALLOCATE(eigenvalues(ndim))
  ALLOCATE(eigvec(ndim))

  OPEN (unit=10,file="EIGVEC",status="old",form="unformatted")
  OPEN (unit=11,file="THIRD",status="old",form="formatted")
  OPEN (unit=12,file="REFERENCE",status="old",form="formatted")
  OPEN (unit=13,file="EIGENVALUES",status="old",form="formatted")

  DO i = 1, ndim
     READ(10) dynmat(:,i)
  ENDDO
  REWIND(10)

  DO i = 1,ndim
    READ(13,*) idx, eigenvalues(i)
  ENDDO
  
  DO i = 1, natms
    READ(12,*) x,y,z, wgt(i)
  ENDDO
  wgt = sqrt(wgt)

! read the third order force constants
  
  nthird = 0
  DO WHILE (.true.)
    READ(11,*,END=101) 
    nthird = nthird + 1
  ENDDO
  101 CONTINUE
  REWIND(11)
  
  print *, '# of anharmonic force constants:', nthird,' *3'
   
  ALLOCATE(atm1(nthird),atm2(nthird),atm3(nthird),c1(nthird),c2(nthird))
  ALLOCATE(v3ijk(3,nthird))

  DO i = 1, nthird
    READ(11,*) atm1(i),c1(i),atm2(i),c2(i),atm3(i),v3ijk(:,i)
  ENDDO
  print *, 'read THIRD'
  RETURN

 END SUBROUTINE readers   

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END PROGRAM relax
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE initmpi(idnode,mxnode)

 include "mpif.h"

 INTEGER idnode,mxnode
 INTEGER ierr

 call MPI_init(ierr)
 idnode = 0
 call MPI_COMM_RANK(MPI_COMM_WORLD, idnode ,ierr)
 mxnode = 1
 call MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode, ierr)

 return

 END SUBROUTINE initmpi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE exitmpi

 include "mpif.h"
 INTEGER ierr

 call MPI_FINALIZE(ierr)
 return
 
 END SUBROUTINE exitmpi
