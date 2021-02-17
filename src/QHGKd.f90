  PROGRAM QHGK

  IMPLICIT NONE

  INTEGER  :: natms, ndim, i, j, k, iargc, idx,ntype, jj, iphskip
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: natt
  REAL(8)                              :: part, partot, pvec(3), pi
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: eigenvalues,eigvec,linewidth
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: partype
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dynmat, diff
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: xxx,yyy,zzz, weight
  REAL(8), PARAMETER  :: invcm = 33.35641d0     ! 1 THz = 33.35641  1/cm
  REAL(8)     :: cell(9), delta, volume, temperature, radius ! it is the diameter
  CHARACTER(LEN=100) :: buffer,feigv  
  LOGICAL :: lsijx = .false., lphskip = .false., lpr =.false.
  LOGICAL :: lsijy = .false.
  LOGICAL :: lsijz = .false.

  pi = dacos(-1.d0)
  temperature = 300.d0 ! Default temperature

  if(iargc()<1) stop "Usage: QHGKd.x <#atoms> <#nspecies> [-P] -bX(Y,Z) -T <temperature> -skip <nph>"
  CALL GETARG(1,buffer)
  READ(buffer,*) natms
  CALL GETARG(2,buffer)
  READ(buffer,*) ntype

  IF (iargc()>2) then 
    do i = 3, iargc()
      CALL GETARG(i,buffer)
      if (INDEX(buffer,'-P' ).NE.0) lpr   = .true. ! compute participation ratio
      if (INDEX(buffer,'-bX').NE.0) lsijx = .true.
      if (INDEX(buffer,'-bY').NE.0) lsijy = .true.
      if (INDEX(buffer,'-bZ').NE.0) lsijz = .true.
      if (INDEX(buffer,'-T').NE.0) then
         CALL GETARG(i+1,buffer) 
         READ(buffer,*) temperature
      endif        
      if (INDEX(buffer,'-skip').NE.0) then
        lphskip = .true.
        CALL GETARG(i+1,buffer)
        READ(buffer,*) iphskip
      endif
    enddo
  endif 

  ALLOCATE (natt(ntype))
  ALLOCATE (partype(ntype))

!!!  if (ntype>1) then
!!!    do i = 3, 2+ntype
!!!      CALL GETARG(i, buffer)
!!!      READ(buffer,*) natt(i-2)
!!!    enddo  
!!!  endif

  ndim = 3*natms

  ALLOCATE(xxx(natms),yyy(natms),zzz(natms),weight(natms))
  ALLOCATE(eigenvalues(ndim))
  ALLOCATE(eigvec(ndim))
  ALLOCATE(dynmat(ndim,ndim))
  ALLOCATE (diff(3, ndim))
  ALLOCATE (linewidth(ndim))

  OPEN (unit=11,file="REFERENCE",form="formatted")
  OPEN (unit=12,file="EIGVEC",form="unformatted")
  OPEN (unit=13,file="EIGENVALUES",form="formatted")
  OPEN (unit=15,file="DECAY",form="formatted")
  
  OPEN (unit=14,file='partratio.dat')    ! participation ratio

  DO i = 1,ndim
    READ(13,*) idx, eigenvalues(i)       ! frequencies are given in 1/cm
  ENDDO

! compute the participation ratio

  if (lpr) then
    DO i = 1, ndim
      READ(12) eigvec(1:ndim)
  
      do jj = 0, ntype
        CALL partratio(jj,part)
        if(jj==0) partot = part
        if(jj>0) partype(jj) = part!*partot
      enddo
      partype = partype/sqrt(sum(partype**2))
      if (ntype>1) then
        write(14,'(g14.6,2x,6g12.4)') eigenvalues(i), partot, partype(1:ntype)
      else
        write(14,'(g14.6,2x,g12.4)') eigenvalues(i), partot
      endif
  
  !   CALL polar(pvec)
  !   write(15,'(g14.6,2x,3g12.4)') eigenvalues(i), pvec
  !   CALL phaseq(phase)
    ENDDO
  endif  

  eigenvalues = eigenvalues*2.d0*pi/invcm ! converts frequencies in rad/ps
  CALL readref
  volume = cell(1)*cell(5)*cell(9)  
  print *, "!!! WARNING: This code works only with orthorhombic cells!"
! CALL fconst(ndim)

  if (lsijx .or. lsijy .or. lsijz ) then
    write(*,*) "Build heatflux matrix elements"
    if(lsijx) CALL buildSij(ndim,dynmat,1)  
    if(lsijy) CALL buildSij(ndim,dynmat,2)  
    if(lsijz) CALL buildSij(ndim,dynmat,3)  
  else
          print *,"!!! Use S_ij from previous calculations"
  endif
 
  CALL direcdiff(ndim, eigenvalues, linewidth)
  
  print *,"!!! Computing QHGK conductivity"
  CALL AFcond(ndim, eigenvalues, temperature)
  
  STOP
  CONTAINS 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE fconst(ndim)

     IMPLICIT NONE

     INTEGER :: ndim
     REAL(8), DIMENSION(ndim,ndim) :: dynmat
     REAL(8) :: dx, dy, dz, dr

! read dynmat
     OPEN (unit=10,file="Dyn.form",form="formatted")
     DO i = 1, ndim
        READ(10,*) dynmat(:,i)
     ENDDO
     CLOSE(10)

     open(10,file='force-const.dat')

     do i =1, ndim
       do j = i + 1, ndim
         dx = xxx((i-1)/3 + 1)-xxx((j-1)/3 + 1)
         dy = yyy((i-1)/3 + 1)-yyy((j-1)/3 + 1)
         dz = zzz((i-1)/3 + 1)-zzz((j-1)/3 + 1)
         dx = dx - NINT(dx/cell(1)) * cell(1)
         dy = dy - NINT(dy/cell(5)) * cell(5)
         dz = dz - NINT(dz/cell(9)) * cell(9)
         dr = sqrt(dx**2+dy**2+dz**2)
         write(10,'(2g14.4,2i5)') dr, abs(dynmat(i,j)), i, j
         
       enddo
     enddo
     close(10)

     END SUBROUTINE fconst
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE AFcond(ndim, omega, temp)
     
     IMPLICIT NONE

     INTEGER :: ndim
     REAL(8) :: omega(ndim)

     REAL(8) kappaQM(3), kappaCL(3)
     REAL(8) x, kboltz, expx, cvx, temp, dtemp, cvqm, cvcl
     REAL(8), PARAMETER :: cnv = 47.992374d0 !conversion THz --> Kelvin

     print *,"Computing thermal conductivity at T (K) =",temp
     kboltz = 1.3806504d1
     dtemp = 1.d0
     omega = omega/2.d0/pi * cnv

!    open (unit=100,file='kappa.dat')
     open (unit=103,file='kappa.dat')   

!    do i = 5, 1700 ! Temperature scale
!      temp = DBLE(i)*dtemp
     kappaCL = 0.d0
     kappaQM = 0.d0
     cvqm = 0.d0
     do j = 4+iphskip, ndim
       x = omega(j)/temp
       expx = exp(x)
       cvx   = x*x*expx/(expx-1.0)**2
       cvqm = cvqm + cvx
       kappaCL(:) = kappaCL(:) + kboltz*diff(:,j)/volume
       kappaQM(:) = kappaQM(:) + kboltz*cvx*diff(:,j)/volume
       write(103,'(10g13.3)') omega(j)/cnv, kappaQM(:), kappaCL(:), kboltz*cvx*diff(:,j)/volume
     enddo 
     write(*,*) temp, kappaQM, kappaCL,  cvqm/volume
!    enddo

     return
     END SUBROUTINE AFcond
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$$     SUBROUTINE diffusivity(ndim, omega, linewidth)
!!$$     
!!$$     IMPLICIT NONE
!!$$    
!!$$     INTEGER :: ndim, ix
!!$$     REAL(8) :: omega(ndim), linewidth(ndim)
!!$$    
!!$$     REAL(8), ALLOCATABLE :: ework(:,:)
!!$$     REAL(8) :: lorentz, omthz, delta, hbar
!!$$     CHARACTER(LEN=5) :: filename
!!$$  
!!$$
!!$$     ALLOCATE (ework(ndim,ndim))    
!!$$     ework = 0.d0
!!$$
!!$$     linewidth(1:3) = 0.d0
!!$$     do i = 4, ndim
!!$$        READ(15,*) omthz, linewidth(i)
!!$$     enddo 
!!$$
!!$$     hbar = 0.65821195d0  ! converts meV to rad/ps
!!$$     linewidth = linewidth/hbar
!!$$
!!$$     do ix = 1, 3
!!$$        write(filename,'(a4,i1)') "sij.", ix
!!$$        open(unit=10,file=filename,form='unformatted')
!!$$        
!!$$        print *,"!!! Reading Sij.",ix
!!$$        DO i = 1, ndim
!!$$           READ(10) dynmat(:,i)
!!$$        ENDDO
!!$$        ! dynmat is now sij!!!     
!!$$     
!!$$        do i = 4, ndim
!!$$          do j = 4,ndim
!!$$            ework(j,i) = ework(j,i) + dynmat(j,i)**2 
!!$$          enddo
!!$$        enddo     
!!$$        close(10)
!!$$     enddo
!!$$     do i = 4, ndim
!!$$        do j = 4,ndim
!!$$           ework(j,i) = ework(j,i)*(omega(i)+omega(j))**2/16.d0/omega(i)/omega(j)
!!$$        enddo
!!$$     enddo
!!$$    
!!$$     diff = 0.d0
!!$$     do i = 4, ndim
!!$$       do j = 4, ndim 
!!$$         delta = linewidth(i)+linewidth(j)
!!$$         lorentz = delta/((omega(i)-omega(j))**2 + delta**2)/pi
!!$$         diff(i) = diff(i) + ework(j,i)*lorentz
!!$$       enddo
!!$$       diff(i) = pi*diff(i)/3.d0/omega(i)**2/100.d0 ! convert to mm^2/s
!!$$     enddo
!!$$   
!!$$     OPEN (unit=88,file='diffusivity')
!!$$     do i = 4, ndim
!!$$        write(88,'(6g14.6)') omega(i)/2.d0/pi, diff(i), linewidth(i)
!!$$     enddo
!!$$     close(88)
!!$$     deallocate(ework)
!!$$    
!!$$     return
!!$$     END SUBROUTINE diffusivity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE buildSij(ndim,dynmat,ix)

     IMPLICIT NONE

     INTEGER :: ndim, ix
     REAL(8), DIMENSION(ndim,ndim) :: dynmat

     INTEGER i, j, k1, k2
     REAL(8) :: drij, dx, dy, dz
     REAL(8), ALLOCATABLE, DIMENSION(:,:) :: ework, ework2
     REAL(8), ALLOCATABLE, DIMENSION(:) :: sij, ewvec
     CHARACTER(LEN=5) :: filename

     write(filename,'(a4,i1)') "sij.", ix
     open(unit=22,file=filename,form='unformatted')
     ALLOCATE(ework(ndim,ndim))
     ALLOCATE(ework2(ndim,ndim))
     ALLOCATE(sij(ndim), ewvec(ndim))

! read dynmat
     OPEN (unit=10,file="Dyn.form",form="formatted")
     DO i = 1, ndim
        READ(10,*) dynmat(:,i)
     ENDDO
     CLOSE(10)
     print*,'read dynmat', ix

! compute dz * dynmat

     do i =1, ndim 
       dynmat(i,i) = 0.d0
       do j = i + 1, ndim
         dx = xxx((i-1)/3 + 1)-xxx((j-1)/3 + 1) 
         dy = yyy((i-1)/3 + 1)-yyy((j-1)/3 + 1) 
         dz = zzz((i-1)/3 + 1)-zzz((j-1)/3 + 1) 
         dx = dx - NINT(dx/cell(1)) * cell(1)
         dy = dy - NINT(dy/cell(5)) * cell(5)
         dz = dz - NINT(dz/cell(9)) * cell(9)
         select case (ix)
         case(1)
           drij = dx
         case(2)
           drij = dy
         case(3)
           drij = dz
         end select
         dynmat(i,j) = (dynmat(i,j)+dynmat(j,i))/2.d0 * drij
         dynmat(j,i) = - dynmat(i,j)
       enddo
     enddo

! read the eigenvectors and store them

     REWIND(12)
     do i =1, ndim
       READ(12) ework(1:ndim,i)
     enddo

! Fast way with blas libraries

    CALL DGEMM('T','N',ndim,ndim,ndim, 1.d0, dynmat , ndim, ework, ndim, 0.d0, ework2, ndim )
    
    
!   ework2 = MATMUL(TRANSPOSE(ework),dynmat) ! this is correct but crashes with big systems

     CALL DGEMM('T','N',ndim,ndim,ndim, 1.d0, ework , ndim, ework2, ndim, 0.d0, dynmat, ndim )
     sij(:) = dynmat(:,1)
!    write(89,*) sij(1:5)

     do i =1, ndim
       
       sij(:) = dynmat(:,i)

!      do j = 1, ndim
!           sij(j) = DOT_PRODUCT( ework2(:,i), ework(:,j) )
!      enddo
!      ewvec = ework2(:,i)
!      sij = matmul(ewvec,ework)
!      CALL DGEMV('T',ndim,ndim,1.d0,ework,ndim,ewvec,1,0.d0,sij,1)

       write(22) sij(:)
     enddo

     close(22)
     DEALLOCATE(ework)
     DEALLOCATE(ework2)
     DEALLOCATE(sij)

     return
     end SUBROUTINE buildSij
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE readref
     
     IMPLICIT NONE

     INTEGER     :: k, j
 
! read reference
  
     DO k = 1, natms
       read(11,*) xxx(k), yyy(k), zzz(k), weight(k)
     ENDDO
     read(11,*) cell
     
     END SUBROUTINE readref
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE phaseq(phase)

     IMPLICIT NONE

     INTEGER     :: k, j, idx
     REAL(8)     :: phase
     
     RETURN    
     END SUBROUTINE phaseq
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE polar(pvec)

     IMPLICIT NONE

     INTEGER     :: k, j, idx
     REAL(8)     :: pvec(3), pjvec(3), norm
     
     pvec = 0.d0
     idx = 0
     DO j = 1, natms
       DO k = 1, 3
         idx = idx + 1
         pjvec(k) = eigvec(idx)
       ENDDO
       norm = sum(pjvec**2)
       if ((i>5).and.(i<201)) write(i+20,'(3f12.5)') pjvec
       pvec = pvec + abs(pjvec)
     ENDDO
     
     END SUBROUTINE polar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE partratio(itype, part)
     
     IMPLICIT NONE

     INTEGER idx, j, k, itype, n1, n2
     REAL(8)     :: ppj, invpart, part



     if (jj == 1) then
        n1 = 0
        n2 = natt(1)
     else
        n1 = SUM(natt(1:jj-1))
        n2 = SUM(natt(1:jj))
     endif

     invpart = 0.d0
     idx = 0
     part = 0.d0

!    total participation ratio

     if (jj == 0) then
       do j = 1, natms
         ppj  = 0.d0
         do k = 1, 3
           idx = idx + 1
           ppj = ppj + eigvec(idx)*eigvec(idx)
         enddo
         invpart = invpart + ppj**2
       enddo
       part = 1.d0/DBLE(natms)/invpart

!    partial participation ratio

     else
       do j = n1+1, n2
         ppj  = 0.d0
         do k = 1, 3
           idx = 3*(j-1)+k
           ppj = ppj + eigvec(idx)*eigvec(idx)
         enddo
         part = part + ppj
       enddo
       part = part*DBLE(natms)/(n2-n1)/2.d0
     endif
    
     RETURN
     END SUBROUTINE partratio

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE direcdiff(ndim, omega, linewidth)

     IMPLICIT NONE

     INTEGER :: ndim, ix
     REAL(8) :: omega(ndim), linewidth(ndim)

     REAL(8), ALLOCATABLE :: ework(:,:)
     REAL(8) :: lorentz, omthz, delta, hbar
     CHARACTER(LEN=5) :: filename


     ALLOCATE (ework(ndim,ndim))

     linewidth(1:3) = 0.d0
     do i = 4, ndim
        READ(15,*) omthz, linewidth(i)
        print *, i, omthz, ndim
     enddo

     hbar = 0.65821195d0  ! converts meV to rad/ps
     linewidth = linewidth/hbar
     diff = 0.d0

     do ix = 1, 3
!       ework = 0.d0
        write(filename,'(a4,i1)') "sij.", ix
        open(unit=10,file=filename,form='unformatted')

        print *,"!!! Reading Sij.",ix
        DO i = 1, ndim
           READ(10) dynmat(:,i)
        ENDDO
        close(10)
        ! dynmat is now sij!!!     

        do i = 4+iphskip, ndim
          do j = 4+iphskip,ndim
            ework(j,i) = dynmat(j,i)**2 *(omega(i)+omega(j))**2/16.d0/omega(i)/omega(j)
          enddo
        enddo

        do i = 4+iphskip, ndim
          do j = 4+iphskip, ndim
            delta = linewidth(i)+linewidth(j)
            lorentz = delta/((omega(i)-omega(j))**2 + delta**2)/pi
            diff(ix,i) = diff(ix,i) + ework(j,i)*lorentz
          enddo
          diff(ix,i) = pi*diff(ix,i)/omega(i)**2/100.d0 ! convert to mm^2/s
        enddo

        write(filename,'(a4,i1)') "dif.", ix
        OPEN (unit=88,file=filename)
        do i = 4+iphskip, ndim
           write(88,'(6g14.6)') omega(i)/2.d0/pi, diff(ix,i), linewidth(i)
        enddo
        close(88)
     enddo
     deallocate(ework)

     return
     END SUBROUTINE direcdiff   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  END PROGRAM QHGK

