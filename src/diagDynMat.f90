  PROGRAM diagmat

  INTEGER  :: natms, ndim, i, iargc, neigv
  REAL(8)                              :: eig, twopi
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dynmat
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: eigenvalues,eigvec,weight
  CHARACTER(LEN=100) :: buffer,feigv  
  LOGICAL ::  leig = .false., lacou ! acoustic sum rule

  if(iargc()<1) stop "Usage: diagDynMat <#atoms>"
  CALL GETARG(1,buffer)
  READ(buffer,*) natms
  
!  if (iargc()>1) then
!    CALL GETARG(2,buffer)
!    READ(buffer,*) neigv
!    leig = .true.
!  endif
  lacou = .false.
  if (iargc()>1) lacou =.true.
  twopi = 2.d0*dacos(-1.d0)
  ndim = 3*natms

  ALLOCATE(dynmat(ndim,ndim))
  ALLOCATE(eigenvalues(ndim))
  ALLOCATE(eigvec(ndim))
  ALLOCATE(weight(ndim))

!  OPEN (unit=12,file="DYNMAT",form="unformatted")
  OPEN (unit=12,file="Dyn.form",form="formatted")

  DO i = 1, ndim
!      READ(12) dynmat(:,i)
     READ(12,*) dynmat(:,i)
!     print *,i
  ENDDO
  
  CALL readweight
  if (lacou)  CALL acoustic(dynmat,ndim,weight)
  CALL diag(dynmat,eigenvalues,ndim)
  
  WRITE(*,*) 'Eigenvalues are in cm^-1'
  OPEN (unit=13,file="EIGENVALUES",form="formatted")
  DO i = 1,ndim
    eig = eigenvalues(i)
    eigenvalues(i) = sqrt(ABS(eig))*eig/ABS(eig)*33.356/twopi
    eig = sqrt(ABS(eig))*eig/ABS(eig)/twopi
    WRITE(13,'(i6,3x,2g16.8)') i, eigenvalues(i), eig
  ENDDO

! ... write out a selected eigenvector

  OPEN(unit=14,file='EIGVEC',form='unformatted')
  DO i = 1, ndim
      WRITE(14) dynmat(1:ndim,i)
  ENDDO
  CLOSE(14) 

! OPEN(unit=54,file='EIGVEC.form',form='formatted')
! DO i = 1, ndim
!     WRITE(54,'(3f12.8)') dynmat(1:ndim,i)
! ENDDO
! CLOSE(54)

  IF(leig) THEN
    eigvec = dynmat(:,neigv)
    WRITE(feigv,'(a5,i5)') "EIGV.",neigv
    !print *,feigv
    OPEN (unit=24,file=feigv,form="formatted")
    DO i = 1,ndim,3
      WRITE(24,'(i6,3f12.5)') i/3+1, eigvec(i:i+2)
    ENDDO
  ENDIF

! CALL mldgau(natms,eigenvalues,dynmat,ndim)

  STOP

  CONTAINS

  subroutine readweight

  integer i
  real(8) :: x,y,z, ww

  open (unit=17,file="REFERENCE")
  do i =1,natms
    read(17,*) x,y,z,ww
    weight(3*(i-1)+1:3*i) = ww
  enddo 
  close(17)

  end subroutine readweight

  END PROGRAM diagmat

!-------------------------------------------------------------------
!  Implementing the acoustic sum rule
!-------------------------------------------------------------------
  SUBROUTINE acoustic(dynmat,ndim,weight)
!-------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: ndim
  REAL(8), DIMENSION(ndim)  :: weight
  REAL(8), DIMENSION(ndim,ndim)  :: dynmat

  INTEGER :: i, j, ix, jx, nat
  REAL(8), ALLOCATABLE, DIMENSION(:) :: acsum
  REAL(8) :: mxdif, sumrulecorr, offdiagsum
 
  nat = ndim/3
  ALLOCATE(acsum(nat))

  mxdif = 0.d0
  print *,'Symmetrize the Dyn Mat'
  do i = 1, ndim
    do j = i+1, ndim
      mxdif=MAX(mxdif, abs(dynmat(j,i)-dynmat(i,j)) )
      dynmat(i,j) = 0.5d0*(dynmat(i,j)+dynmat(j,i))
      dynmat(j,i) = dynmat(i,j)
    enddo  
  enddo
  print *,'Maximum symm correction: ', mxdif

! apply the acoustic sum rule

      sumrulecorr = 0.d0
      DO i = 1, ndim
        offdiagsum = 0.d0
        ix = mod(i,3)
        if (ix==0) ix = 3
        do j = ix, ndim, 3
          offdiagsum = offdiagsum+dynmat(j,i)*sqrt(weight(i)*weight(j)) 
        enddo
        dynmat(i,i) = dynmat(i,i) - offdiagsum/weight(i)
        sumrulecorr = max(sumrulecorr,offdiagsum)
      ENDDO

      WRITE(*,'("Max. sum rule correction:",g12.5)') sumrulecorr

  RETURN
  END SUBROUTINE acoustic
!-------------------------------------------------------------------
!  Driver to call blas-lapack diagonalisation libraries
!-------------------------------------------------------------------
  SUBROUTINE diag(itmp, d, n)

! itmp : IN  -> input tensor  OUT -> eigenvectors
! d    : OUT -> eigenvalues
! n    : IN  -> size of the matrix

  INTEGER    :: N, lwork
  REAL(8), DIMENSION(N,N)  :: itmp

  INTEGER                  :: INFO
  REAL(8)                  :: D(n), E(n-1), TAU(n-1)
  REAL(8), ALLOCATABLE, DIMENSION(:) :: work
  CHARACTER(1)             :: uplo = 'U'

  lwork = 3*n-1

  print *,'Entering lapack routines for diagonalization'
  ALLOCATE (work(lwork))

  CALL DSYEV('V','U',N,itmp,N,D,work,lwork,INFO)

! CALL DSYTRD(uplo, N, itmp, N, D, E, TAU, WORK, LWORK, INFO )
! CALL DORGTR(uplo, N, itmp, N, TAU, WORK, LWORK,INFO)
! print *, '\2'
! CALL DSTEQR('V', N, D, E, itmp, N, WORK, INFO )
! print *, '\3'

  IF (INFO.NE.0) print*, "Error!"
  RETURN
  END SUBROUTINE diag
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE MLDGAU(NAT,VIBE,SDER,NDIM)
!     ==--------------------------------------------------------------==
!     write gaussianformated output into VIB1.log and VIB2.log files,
!     which are readable in molden/molekel to visualise the vibrations.
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
!     Arguments
      INTEGER   NDIM, NAT
      REAL*8    VIBE(NDIM),SDER(NDIM,NDIM)
!     Variables
      LOGICAL   LW
      CHARACTER FILE1*100,FILE2*100,typ*2,a*8,ck(3)*7,fr*17,rm(5)*17
      INTEGER   i,j,k,l,n,t(NAT),is,ia,at,he,IUNIT,m
      REAL*8    gz,b,dbltk
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xc,yc,zc
!     ==--------------------------------------------------------------==
 
      ALLOCATE(xc(nat),yc(nat),zc(nat))

      FILE1='VIB1.log'

      fr=' Frequencies --  '
      rm(1)=' Red. masses --  '
      rm(2)=' Frc consts  --  '
      rm(3)=' IR Inten    --  '
      rm(4)=' Raman Activ --  '
      rm(5)=' Depolar     --  '
      typ='?A'
      a=' Atom AN'
      ck(1)='      X'
      ck(2)='      Y'
      ck(3)='      Z'
      at=0
      b=0.529177249
      gz=0.
      LW=.false.
! ---- open files ---------------------------
      OPEN (UNIT=115, FILE=FILE1, ERR=999)
      OPEN (UNIT=17, FILE='REFERENCE',status='old', ERR=998)
! ---- read coordinates and atomtyps ------------
      DO k=1,NAT
           READ(17,*) xc(k), yc(k), zc(k), dbltk
           t(k) = INT(dbltk/2)
           if(t(k)==0) t(k) =1
      ENDDO
! ---- write some lines needed for molden
      IUNIT=115
      WRITE (IUNIT,*)'Entering Gaussian System'
      WRITE (IUNIT,*)'this file is generated from the MLDGAU', &
                     ' subroutine in the file secder.F'
      WRITE (IUNIT,*)'Standard orientation:'
      WRITE (IUNIT,*)'---------------------------------------', &
                     '------------------------------'
      WRITE (IUNIT,'(A,2(5X,A),14X,A)') &
                'Center','Atomic','Atomic','Coordinates (Angstroms)'
      WRITE (IUNIT,'(2(A,5X),1X,A,3X,3(11X,A))') &
                    'Number','Number','Type','X','Y','Z'
      WRITE (IUNIT,*)'---------------------------------------', &
                    '------------------------------'
! ---- make sure that the center of mass is the origin or move the atoms
! ---- from the center of the box to the origin and printout
      DO  i=1,NAT
         WRITE (115,22)  i, t(i), at, xc(i), yc(i), zc(i)
      ENDDO
! ---- write some lines for molden -----------------------------
      
      WRITE(IUNIT,*)'--------------------------------------------', &
                    '-------------------------'
      WRITE(IUNIT,*)'      basis functions          primitive ', &
                    'gaussians'
      WRITE(IUNIT,*)'      alpha electrons          beta electrons' 
      WRITE(IUNIT,*)'********************************************', &
                    '**************************'
      WRITE(IUNIT,*)
      WRITE(IUNIT,*)'Harmonic frequencies (cm**-1), IR intensities ', &
                    '(KM/Mole),'
      WRITE(IUNIT,*)'Raman scattering activities (A**4/AMU), Raman ', &
                    'depolarization ratios,'
      WRITE(IUNIT,*)'reduced masses (AMU), force constants ', &
                    '(mDyne/A) and normal coordinates:'
! ---- write eigenvalues and eigenvectors in both files
      DO 15 i=1,(NDIM),3
        WRITE(115,23) i-3, i-2, i-1
        WRITE(115,24) typ, typ, typ
        WRITE(115,25) fr, (VIBE(l),l=i,i+2)
        DO n=1,5
           WRITE(115,25) rm(n), gz, gz, gz
        ENDDO
        WRITE(115,26) a,(ck(n),n=1,3),(ck(n),n=1,3),(ck(n),n=1,3)
        DO 16 j=1,NDIM,3
          he=(j-1)/3+1
          WRITE(115,27) he,t(he), &
                       (SDER(j,m),SDER(j+1,m),SDER(j+2,m),m=i,i+2)
 16     CONTINUE
 15   CONTINUE
      WRITE(115,*) 'Normal termination of Gaussian 98.'
      CLOSE (115)
 22   FORMAT(i5,i11,i14,4x,3(3x,f10.6))
 23   FORMAT(i22,2i23)
 24   FORMAT(20x,a2,2(21x,a2))
 25   FORMAT(a17,f9.4,2f23.4)
 26   FORMAT(a8,3a7,2(2x,3a7))
 27   FORMAT(2i4,3(f9.2,2f7.2))
      goto 888
 999  CONTINUE
       WRITE(6,*) 'COULD NOT OPEN FILE VIB[1,2].log!'
 888  CONTINUE
!     ==--------------------------------------------------------------==
      RETURN
 998  stop "REFERENCE FILE MISSING"
      END
!     ==================================================================


 
