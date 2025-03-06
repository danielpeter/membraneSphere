! numerical recipies functions

!-----------------------------------------------------------------------
! http://www.library.cornell.edu/nr/cbookf90pdf.html, see chapter B13
! or
! http://hebb.mit.edu/courses/9.29/2002/readings/c13-2.pdf
! or
! http://www.phys.uu.nl/DU/num_recipes/fortran.208/f90/recipes/


! correl.f90
! new single and double precision routines (daniel peter,22.5.2005)
!-----------------------------------------------------------------------
      function correl_sp(data1,data2)       
!-----------------------------------------------------------------------
! Computes the correlation of two real datasets data1and data2 of length N (includ- 
! ing any user-supplied zeropadding). N must be an integer power of 2. The answer is 
! returned as the function correl, an array of length N. The answer is stored in wrap- 
! around order, i.e., correlations at increasingly negative lags are in correl(N) on down to 
! correl(N/2+1), while correlations at increasingly positive lags are in correl(1) (zero 
! lag) on up to correl(N/2). Sign convention of this routine: if data1 lags data2, i.e., 
! is shifted to the right of it, then correl will show a peak at positive lags. 
      USE nrtype; USE nrutil, ONLY : assert,assert_eq
      USE nr, ONLY : realft
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: data1,data2
      REAL(SP), DIMENSION(size(data1)) :: correl_sp
      COMPLEX(SPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
      INTEGER(I4B) :: no2,n

      !print*,'processing correlation...'      
      n=assert_eq(size(data1),size(data2),'correl_sp')
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl_sp')
      no2=n/2
      
      call realft(data1,1,cdat1)
      call realft(data2,1,cdat2)
      
      cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/no2, &
        aimag(cdat1(1))*aimag(cdat2(1))/no2, kind=spc)
      cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2      
      call realft(correl_sp,-1,cdat1)
      
      END FUNCTION correl_sp


!-----------------------------------------------------------------------
      function correl_dp(data1,data2)
!-----------------------------------------------------------------------
! Computes the correlation of two real datasets data1and data2 of length N (includ- 
! ing any user-supplied zeropadding). N must be an integer power of 2. The answer is 
! returned as the function correl, an array of length N. The answer is stored in wrap- 
! around order, i.e., correlations at increasingly negative lags are in correl(N) on down to 
! correl(N/2+1), while correlations at increasingly positive lags are in correl(1) (zero 
! lag) on up to correl(N/2). Sign convention of this routine: if data1 lags data2, i.e., 
! is shifted to the right of it, then correl will show a peak at positive lags.       
      USE nrtype; USE nrutil, ONLY : assert,assert_eq
      USE nr, ONLY : realft
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: data1,data2
      REAL(DP), DIMENSION(size(data1)) :: correl_dp
      COMPLEX(DPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
      INTEGER(I4B) :: no2,n
      
      !print*,'checking...'
      n=assert_eq(size(data1),size(data2),'correl_dp')
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl_dp')
      no2=n/2
      
      !print*,'fourier transform...'
      call realft(data1,1,cdat1)
      call realft(data2,1,cdat2)
      
      cdat1(1)=cmplx(dreal(cdat1(1))*dreal(cdat2(1))/no2,dimag(cdat1(1))*dimag(cdat2(1))/no2, kind=dpc)
      cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2
      
      !print*,'back fourier transform...'
      call realft(correl_dp,-1,cdat1)
      
      END FUNCTION correl_dp

!-----------------------------------------------------------------------
! http://www.phys.uu.nl/DU/num_recipes/fortran.208/f90/recipes/
!
!-----------------------------------------------------------------------
      FUNCTION convlv(data,respns,isign)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,nrerror
      USE nr, ONLY : realft
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
      REAL(SP), DIMENSION(:), INTENT(IN) :: respns
      INTEGER(I4B), INTENT(IN) :: isign
      REAL(SP), DIMENSION(size(data)) :: convlv
      INTEGER(I4B) :: no2,n,m
      COMPLEX(SPC), DIMENSION(size(data)/2) :: tmpd,tmpr
      n=size(data)
      m=size(respns)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in convlv')
      call assert(mod(m,2)==1, 'm must be odd in convlv')
      convlv(1:m)=respns(:)
      convlv(n-(m-3)/2:n)=convlv((m+3)/2:m)
      convlv((m+3)/2:n-(m-1)/2)=0.0
      no2=n/2
      call realft(data,1,tmpd)
      call realft(convlv,1,tmpr)
      if (isign == 1) then
        tmpr(1)=cmplx(real(tmpd(1))*real(tmpr(1))/no2, &
          aimag(tmpd(1))*aimag(tmpr(1))/no2, kind=spc)
        tmpr(2:)=tmpd(2:)*tmpr(2:)/no2
      else if (isign == -1) then
        if (any(abs(tmpr(2:)) == 0.0) .or. real(tmpr(1)) == 0.0 &
          .or. aimag(tmpr(1)) == 0.0) call nrerror &
          ('deconvolving at response zero in convlv')
        tmpr(1)=cmplx(real(tmpd(1))/real(tmpr(1))/no2, &
          aimag(tmpd(1))/aimag(tmpr(1))/no2, kind=spc)
        tmpr(2:)=tmpd(2:)/tmpr(2:)/no2
      else
        call nrerror('no meaning for isign in convlv')
      end if
      call realft(convlv,-1,tmpr)
      END FUNCTION convlv


!-----------------------------------------------------------------------
! four1.f90

!-----------------------------------------------------------------------
      SUBROUTINE four1_sp(data,isign)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : arth,assert
      USE nr, ONLY : fourrow
      IMPLICIT NONE
      COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
! Replaces a complex array data by its discrete Fourier transform, if isign is input as 1; 
! or replaces data by its inverse discrete Fourier transform times the size of data, if isign 
! is input as −1. The size of data must be an integer power of 2. Parallelism is achieved 
! byinternallyreshapingtheinput arraytotwodimensions. (Usethisversionif fourrowis 
! faster thanfourcolonyour machine.)      
      COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
      COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
      REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
      INTEGER(I4B) :: n,m1,m2,j
      n=size(data)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp')
      m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
      m2=n/m1
      allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
      dat=reshape(data,shape(dat))
      call fourrow(dat,isign)
      theta=arth(0,isign,m1)*TWOPI_D/n
      wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
      w=cmplx(1.0_dp,0.0_dp,kind=dpc)
      do j=2,m2
        w=w*wp+w
        dat(:,j)=dat(:,j)*w
      end do
      temp=transpose(dat)
      call fourrow(temp,isign)
      data=reshape(temp,shape(data))
      deallocate(dat,w,wp,theta,temp)
      END SUBROUTINE four1_sp

!-----------------------------------------------------------------------
      subroutine four1_dp(data,isign)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : arth,assert
      USE nr, ONLY : fourrow
      IMPLICIT NONE
      COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
      COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
      COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
      REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
      INTEGER(I4B) :: n,m1,m2,j
      n=size(data)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
      m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
      m2=n/m1
      allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
      dat=reshape(data,shape(dat))
      call fourrow(dat,isign)
      theta=arth(0,isign,m1)*TWOPI_D/n
      wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
      w=cmplx(1.0_dp,0.0_dp,kind=dpc)
      do j=2,m2
        w=w*wp+w
        dat(:,j)=dat(:,j)*w
      end do
      temp=transpose(dat)
      call fourrow(temp,isign)
      data=reshape(temp,shape(data))
      deallocate(dat,w,wp,theta,temp)
      END SUBROUTINE four1_dp

!-----------------------------------------------------------------------
! fourrow.f90

!-----------------------------------------------------------------------
      SUBROUTINE fourrow_sp(data,isign)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,swap
      IMPLICIT NONE
      COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
! Replaceseachrow(constantfirstindex)ofdata(1:M,1:N)byitsdiscreteFouriertrans- 
! form(transformonsecondindex), if isignis input as 1; or replaces eachrowof data 
! byNtimesitsinversediscreteFourier transform, if isignisinput as −1. Nmust bean 
! integer power of 2. ParallelismisM-foldonthefirst indexof data.       
      INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
      REAL(DP) :: theta
      COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
      COMPLEX(DPC) :: w,wp
      COMPLEX(SPC) :: ws
      n=size(data,2)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
      n2=n/2
      j=n2
      do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
        end do
        j=j+m
      end do
      mmax=1
      do
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI_D/(isign*mmax)
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
        w=cmplx(1.0_dp,0.0_dp,kind=dpc)
        do m=1,mmax
          ws=w
          do i=m,n,istep
            j=i+mmax
            temp=ws*data(:,j)
            data(:,j)=data(:,i)-temp
            data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
        end do
        mmax=istep
      end do
      END SUBROUTINE fourrow_sp

!-----------------------------------------------------------------------
      SUBROUTINE fourrow_dp(data,isign)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,swap
      IMPLICIT NONE
      COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
      INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
      REAL(DP) :: theta
      COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
      COMPLEX(DPC) :: w,wp
      COMPLEX(DPC) :: ws
      n=size(data,2)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
      n2=n/2
      j=n2
      do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
        end do
        j=j+m
      end do
      mmax=1
      do
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI_D/(isign*mmax)
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
        w=cmplx(1.0_dp,0.0_dp,kind=dpc)
        do m=1,mmax
          ws=w
          do i=m,n,istep
            j=i+mmax
            temp=ws*data(:,j)
            data(:,j)=data(:,i)-temp
            data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
        end do
        mmax=istep
      end do
      END SUBROUTINE fourrow_dp

!-----------------------------------------------------------------------
! realft.f90

!-----------------------------------------------------------------------
      SUBROUTINE realft_sp(data,isign,zdata)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,assert_eq,zroots_unity
      USE nr, ONLY : four1
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
      COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
! Whenisign=1, calculates theFourier transformof aset of Nreal-valueddatapoints, 
! inputinthearraydata. Iftheoptional argumentzdataisnotpresent, thedataarereplaced 
! bythepositivefrequencyhalf of its complexFourier transform. The real-valued first and 
! last components of the complex transform are returned as elements data(1)and data(2), 
! respectively. If thecomplexarrayzdataof lengthN/2ispresent, dataisunchangedand 
! thetransformisreturnedinzdata. Nmustbeapowerof 2. If isign= −1, thisroutine 
! calculatestheinversetransformof acomplexdataarrayif itisthetransformof real data. 
! (Resultinthiscasemustbemultipliedby2/N.) Thedatacanbesuppliedeither indata, 
! withzdataabsent, or inzdata. 
      INTEGER(I4B) :: n,ndum,nh,nq
      COMPLEX(SPC), DIMENSION(size(data)/4) :: w
      COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
      COMPLEX(SPC), DIMENSION(:), POINTER :: cdata
      COMPLEX(SPC) :: z
      REAL(SP) :: c1=0.5_sp,c2
      n=size(data)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
      nh=n/2
      nq=n/4
      if (present(zdata)) then
        ndum=assert_eq(n/2,size(zdata),'realft_sp')
        cdata=>zdata
        if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      else
        allocate(cdata(n/2))
        cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      end if
      if (isign == 1) then
        c2=-0.5_sp
        call four1(cdata,+1)
      else
        c2=0.5_sp
      end if
      w=zroots_unity(sign(n,isign),n/4)
      w=cmplx(-aimag(w),real(w),kind=spc)
      h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
      h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
      cdata(2:nq)=h1+w(2:nq)*h2
      cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
      z=cdata(1)
      if (isign == 1) then
        cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
      else
        cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
        call four1(cdata,-1)
      end if
      if (present(zdata)) then
        if (isign /= 1) then
          data(1:n-1:2)=real(cdata)
          data(2:n:2)=aimag(cdata)
        end if
      else
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
        deallocate(cdata)
      end if
      END SUBROUTINE realft_sp


!-----------------------------------------------------------------------
      SUBROUTINE realft_dp(data,isign,zdata)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,assert_eq,zroots_unity
      USE nr, ONLY : four1
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
      INTEGER(I4B), INTENT(IN) :: isign
      COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
      INTEGER(I4B) :: n,ndum,nh,nq
      COMPLEX(DPC), DIMENSION(size(data)/4) :: w
      COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
      COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
      COMPLEX(DPC) :: z
      REAL(DP) :: c1=0.5_dp,c2
      n=size(data)
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
      nh=n/2
      nq=n/4
      if (present(zdata)) then
        ndum=assert_eq(n/2,size(zdata),'realft_dp')
        cdata=>zdata
        if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      else
        allocate(cdata(n/2))
        cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
      end if
      if (isign == 1) then
        c2=-0.5_dp
        call four1(cdata,+1)
      else
        c2=0.5_dp
      end if
      w=zroots_unity(sign(n,isign),n/4)
      w=cmplx(-aimag(w),real(w),kind=dpc)
      h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
      h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
      cdata(2:nq)=h1+w(2:nq)*h2
      cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
      z=cdata(1)
      if (isign == 1) then
        cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
      else
        cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
        call four1(cdata,-1)
      end if
      if (present(zdata)) then
        if (isign /= 1) then
          data(1:n-1:2)=real(cdata)
          data(2:n:2)=aimag(cdata)
        end if
      else
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
        deallocate(cdata)
      end if
      END SUBROUTINE realft_dp

!-----------------------------------------------------------------------
! twofft.f90

!-----------------------------------------------------------------------
      SUBROUTINE twofft(data1,data2,fft1,fft2)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert,assert_eq
      USE nr, ONLY : four1
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
      COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
      INTEGER(I4B) :: n,n2
      COMPLEX(SPC), PARAMETER :: C1=(0.5_sp,0.0_sp), C2=(0.0_sp,-0.5_sp)
      COMPLEX, DIMENSION(size(data1)/2+1) :: h1,h2
      n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft')
      call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft')
      fft1=cmplx(data1,data2,kind=spc)
      call four1(fft1,1)
      fft2(1)=cmplx(aimag(fft1(1)),0.0_sp,kind=spc)
      fft1(1)=cmplx(real(fft1(1)),0.0_sp,kind=spc)
      n2=n/2+1
      h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1)))
      h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
      fft1(2:n2)=h1(2:n2)
      fft1(n:n2:-1)=conjg(h1(2:n2))
      fft2(2:n2)=h2(2:n2)
      fft2(n:n2:-1)=conjg(h2(2:n2))
      END SUBROUTINE twofft


!-----------------------------------------------------------------------
! numerical recipies
! dfridr.f90      
      FUNCTION dfridr(func,x,h,err)
!-----------------------------------------------------------------------
! Returns the derivative of a function func at a point x by Ridders' method of polynomial
! extrapolation. The value h is input as an estimated initial stepsize; it need not be small,
! but rather should be an increment in x over which func changes substantially. An estimate
! of the error in the derivative is returned as err.
! Parameters: Stepsize is decreased by CON at each iteration. Max size of tableau is set by
! NTAB. Return when error is SAFE worse than the best so far.      
      USE nrtype; USE nrutil, ONLY : assert,geop,iminloc
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: x,h
      REAL(DP), INTENT(OUT) :: err
      REAL(DP) :: dfridr
      INTERFACE
        FUNCTION func(x)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: x
        REAL(DP) :: func
        END FUNCTION func
      END INTERFACE
      INTEGER(I4B),PARAMETER :: NTAB=10
      REAL(SP), PARAMETER :: CON=1.4_sp
      REAL(SP), PARAMETER :: CON2=CON*CON
      REAL(SP), PARAMETER :: BIG= 1.e+30 !huge(real(x)) ! (e.g. XLF or SUN): BIG=huge(real(x))
      REAL(SP), PARAMETER :: SAFE=2.0
      INTEGER(I4B) :: ierrmin,i,j
      REAL(SP) :: hh
      REAL(SP), DIMENSION(NTAB-1) :: errt,fac
      REAL(SP), DIMENSION(NTAB,NTAB) :: a
            
      call assert(h /= 0.0, 'dfridr arg')
      hh=h      
      a(1,1)=(func(dble(x+hh))-func(dble(x-hh)))/(2.0_sp*hh)
      err=BIG
      fac(1:NTAB-1)=geop(CON2,CON2,NTAB-1)
      do i=2,NTAB
        hh=hh/CON
        a(1,i)=(func(dble(x+hh))-func(dble(x-hh)))/(2.0_sp*hh)
        do j=2,i
          a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.0_sp)
        end do
        errt(1:i-1)=max(abs(a(2:i,i)-a(1:i-1,i)),abs(a(2:i,i)-a(1:i-1,i-1)))
        ierrmin=iminloc(errt(1:i-1))
        if (errt(ierrmin) <= err) then
          err=errt(ierrmin)
          dfridr=dble(a(1+ierrmin,i))
        end if
        if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) RETURN
      end do
      END FUNCTION dfridr

!-----------------------------------------------------------------------
      double precision FUNCTION dfridr_spline(xin,h,err)
!-----------------------------------------------------------------------
! calculates the first derivative of the function which is stored in splineFunction module
! at the specified location 
! 
! input:
!     xin     - location to caluculate value for
!     h       - distance
!     err     - error
!
! returns: first derivated function value at xin, error estimation in err
      use splineFunction, only: i1,i2,X,Y,Q,ilength
      USE nrtype; USE nrutil, ONLY : assert,geop,iminloc
      use parallel
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: xin,h
      REAL(DP), INTENT(OUT) :: err

      !double precision,external:: splineRepresentation
      double precision,external:: drsple
      
      ! parameter
      INTEGER(I4B),PARAMETER :: NTAB=10
      REAL(SP), PARAMETER :: CON=1.4_sp
      REAL(SP), PARAMETER :: CON2=CON*CON
      REAL(SP), PARAMETER :: BIG= 1.e+30 !huge(real(xin)) ! (e.g. XLF or SUN): BIG=huge(real(xin))
      REAL(SP), PARAMETER :: SAFE=2.0
      INTEGER(I4B) :: ierrmin,i,j
      REAL(SP) :: hh
      REAL(SP), DIMENSION(NTAB-1) :: errt,fac
      REAL(SP), DIMENSION(NTAB,NTAB) :: a

      call assert(h /= 0.0, 'dfridr_spline arg')
      hh=h
      err = 0.0
            
      !a(1,1)=(splineRepresentation(dble(xin+hh))-splineRepresentation(dble(xin-hh)))/(2.0_sp*hh)
      a(1,1)=(drsple(i1,i2,X,Y,Q,ilength,dble(xin+hh))-drsple(i1,i2,X,Y,Q,ilength,dble(xin-hh)))/(2.0_sp*hh)

      err=BIG
      fac(1:NTAB-1)=geop(CON2,CON2,NTAB-1)
      do i=2,NTAB
        hh=hh/CON
        
        !a(1,i)=(splineRepresentation(dble(xin+hh))-splineRepresentation(dble(xin-hh)))/(2.0_sp*hh)
        a(1,i)=(drsple(i1,i2,X,Y,Q,ilength,dble(xin+hh))-drsple(i1,i2,X,Y,Q,ilength,dble(xin-hh)))/(2.0_sp*hh)
        
        do j=2,i
          a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.0_sp)
        end do
        errt(1:i-1)=max(abs(a(2:i,i)-a(1:i-1,i)),abs(a(2:i,i)-a(1:i-1,i-1)))
        ierrmin=iminloc(errt(1:i-1))
        if (errt(ierrmin) <= err) then
          err=errt(ierrmin)
          dfridr_spline = dble(a(1+ierrmin,i))
        end if
        if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) then
         RETURN
        endif
      end do
      
      END FUNCTION dfridr_spline


!-----------------------------------------------------------------------
	    double precision function splineRepresentation(location)
!-----------------------------------------------------------------------
      !use nrtype
      use splineFunction,only: i1,i2,X,Y,Q,ilength
      implicit none
      double precision, INTENT(IN) :: location
      double precision,external:: drsple
      ! calculate value for this representation
      splineRepresentation = drsple(i1,i2,X,Y,Q,ilength,location)      
      return
      end function splineRepresentation

!-----------------------------------------------------------------------
      FUNCTION spline_func(trace)
!-----------------------------------------------------------------------
      use splineFunction,only: i1,i2,X,Y,Q,ilength
      !USE nrtype 
      IMPLICIT NONE
      double precision, DIMENSION(:),INTENT(IN) :: trace
      double precision, DIMENSION(size(trace)) :: spline_func
      integer:: i
      double precision,external:: drsple
      !interface
      !  function splineRepresentation(location)
      !  USE nrtype
      !  IMPLICIT NONE
      !  REAL(DP),INTENT(IN) :: location
      !  REAL(DP):: splineRepresentation
      !  END FUNCTION splineRepresentation
      !end interface
      
      do i=1,size(trace)
        !spline_func(i) = splineRepresentation(x(i))        
        spline_func(i) = drsple(i1,i2,X,Y,Q,ilength,trace(i))
      enddo
      END FUNCTION spline_func




!-----------------------------------------------------------------------
      FUNCTION qsimp(func,a,b)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : nrerror
      USE nr, ONLY : trapzd
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: a,b
      REAL(DP) :: qsimp
      INTERFACE
        FUNCTION func(x)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
      END INTERFACE
      INTEGER(I4B), PARAMETER :: JMAX=20
      REAL(DP), PARAMETER :: EPS=1.0e-6_dp
      INTEGER(I4B) :: j
      REAL(DP) :: os,ost,st
      ost=0.0
      os= 0.0    
      do j=1,JMAX
        call trapzd(func,a,b,st,j)
        qsimp=(4.0_dp*st-ost)/3.0_sp
        if (j > 5) then
          if (abs(qsimp-os) < EPS*abs(os) .or. &
            (qsimp == 0.0 .and. os == 0.0)) then 
              !print*,'qsimp:',qsimp,a,b              
              RETURN
          end if
        end if
        os=qsimp
        ost=st
      end do      
      call nrerror('qsimp: too many steps')
      END FUNCTION qsimp
          
!-----------------------------------------------------------------------
      SUBROUTINE trapzd(func,a,b,s,n)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : arth
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: a,b
      REAL(DP), INTENT(INOUT) :: s
      INTEGER(I4B), INTENT(IN) :: n
      INTERFACE
        FUNCTION func(x)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
      END INTERFACE
      REAL(DP) :: del,fsum
      INTEGER(I4B) :: it
      if (n == 1) then
        s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
      else
        it=2**(n-2)
        del=(b-a)/it
        fsum=sum(func(arth(a+0.5_dp*del,del,it)))
        s=0.5_dp*(s+del*fsum)
      end if
      END SUBROUTINE trapzd


!-----------------------------------------------------------------------
      SUBROUTINE polint(xa,ya,x,y,dy)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(OUT) :: y,dy
      INTEGER(I4B) :: m,n,ns
      REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho

      n=assert_eq(size(xa),size(ya),'polint')
      c=ya
      d=ya
      ho=xa-x
      ns=iminloc(sngl(abs(x-xa)))
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)
        if (any(den(1:n-m) == 0.0)) &
          call nrerror('polint: calculation failure')
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m)
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        end if
        y=y+dy
      end do
      END SUBROUTINE polint

!-----------------------------------------------------------------------
      FUNCTION qromb(func,a,b)
!-----------------------------------------------------------------------
      USE nrtype; USE nrutil, ONLY : nrerror
      USE nr, ONLY : polint,trapzd
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: a,b
      REAL(DP) :: qromb
      INTERFACE
        FUNCTION func(x)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
      END INTERFACE
      INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
      REAL(DP), PARAMETER :: EPS=1.0e-6_dp
      REAL(DP), DIMENSION(JMAXP) :: h,s
      REAL(DP) :: dqromb
      INTEGER(I4B) :: j
      h(1)=1.0
      do j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j >= K) then
          call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
          if (abs(dqromb) <= EPS*abs(qromb)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j)
      end do
      call nrerror('qromb: too many steps')
      END FUNCTION qromb




!-----------------------------------------------------------------------
  SUBROUTINE amoeba(p,y,ftol,func,iter) 
!-----------------------------------------------------------------------
  ! http://www.physics.louisville.edu/help/nr/bookf90pdf/chap10f9.pdf
  !
  ! Minimization of the function func in N dimensions by the downhill simplex method of Nelder and Mead. 
  ! The(N+1)xN matrix p is input. Its N+1 rows are N-dimensional 
  ! vectors that are the vertices of the starting simplex. Also input is the vector y of length 
  ! N+1, whose components must be preinitialized to the values of func evaluated at the 
  ! N+1 vertices(rows) of p; and ftol the fractional convergence tolerance to be achieved 
  ! in the function value(n.b.!). On output, p and y will have been reset to N+1 new points 
  ! all within ftol of a minimum function value, and iter gives the number of function 
  ! evaluations taken. 
  !
    USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap 
    IMPLICIT NONE 
    INTEGER(I4B), INTENT(OUT) :: iter 
    REAL(SP), INTENT(IN) :: ftol 
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: y 
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p 
    INTERFACE 
      FUNCTION func(x) 
        USE nrtype 
        IMPLICIT NONE 
        REAL(SP), DIMENSION(:), INTENT(IN) :: x 
        REAL(SP) :: func 
      END FUNCTION func 
    END INTERFACE 

    ! Parameters: The maximum allowed number of function evaluations, and a small number. 
    INTEGER(I4B), PARAMETER :: ITMAX=500 ! org: 5000 
    REAL(SP), PARAMETER :: TINY=1.0e-10 

    ! Global variables. 
    INTEGER(I4B) :: ihi,ndim 
    REAL(SP), DIMENSION(size(p,2)) :: psum 

    call amoeba_private 

  CONTAINS 

    SUBROUTINE amoeba_private 
      IMPLICIT NONE 
      INTEGER(I4B) :: i,ilo,inhi 
      REAL(SP) ::rtol,ysave,ytry,ytmp 

      ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba') 
      iter=0 
      psum(:)=sum(p(:,:),dim=1) 

      do 
        !Iterationloop. 
        ilo=iminloc(y(:)) 
        !Determine which point is the highest(worst), 
        !next-highest, and lowest(best).

        ihi=imaxloc(y(:)) 
        ytmp=y(ihi) 
        y(ihi)=y(ilo) 
        inhi=imaxloc(y(:)) 
        y(ihi)=ytmp 

        !Compute the fractional range from highest to lowest and return if satisfactory.         
        rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY) 
        
        !! maybe a correction: http://lib.stat.cmu.edu/S/recipes
        !!rtol=2.0_sp*abs(y(ihi))*(abs(y(ihi))-abs(y(ilo)))/(abs(y(ihi))+abs(y(ilo)))
        
        if (rtol < ftol) then 
          !If returning, put best point and value in slot 1
          call swap(y(1),y(ilo)) 
          call swap(p(1,:),p(ilo,:)) 
          RETURN 
        end if 
        if (iter >= ITMAX) call nrerror('ITMAX exceeded inamoeba') 

        !Begin a new iteration. First extrapolate by a factor 
        !_1 through the face of the simplex 
        !across from the highpoint, i.e., reflect the simplex from the highpoint. 

        ytry=amotry(-1.0_sp) 
        iter=iter+1 
        if (ytry <= y(ilo)) then 

          !Gives a result better than the bestpoint, so 
          !try an additional extrapolation by a factor of 2. 

          ytry=amotry(2.0_sp) 
          iter=iter+1 
        else if (ytry >= y(inhi)) then 

          !The reflected point is worse than the sec- 
          !ond highest, so look for an intermediate 
          !lower point, i.e., doa one-dimensional 
          !contraction. 

          ysave=y(ihi) 
          ytry=amotry(0.5_sp) 
          iter=iter+1

          if (ytry >= ysave) then 

            !Can't seemtoget ridof that highpoint. Better contract aroundthe lowest 
            !(best) point. 

            p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1))) 
            do i=1,ndim+1 
              if (i /= ilo) y(i)=func(p(i,:)) 
            end do 
            iter=iter+ndim 

            !Keep track of function evaluations. 

            psum(:)=sum(p(:,:),dim=1) 
          end if 
        end if 
      end do 

      !Go back for the test of doneness and the next iteration.

    END SUBROUTINE amoeba_private 


!-----------------------------------------------------------------------
    FUNCTION amotry(fac) 
!-----------------------------------------------------------------------
      !Extrapolates by a factor fac through the face of the simplex across from the highpoint, 
      !tries it, and replaces the highpoint if the newpoint is better. 
      IMPLICIT NONE 
      REAL(SP), INTENT(IN) :: fac 
      REAL(SP) ::amotry 
      REAL(SP) ::fac1,fac2,ytry 
      REAL(SP), DIMENSION(size(p,2)) :: ptry 

      fac1=(1.0_sp-fac)/ndim 
      fac2=fac1-fac 
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2 
      !Evaluate the function at the trial point. 
      ytry=func(ptry) 


      ! If it's better than the highest, then replace the highest.
      if (ytry <y(ihi)) then 
        y(ihi)=ytry 
        psum(:)=psum(:)-p(ihi,:)+ptry(:) 
        p(ihi,:)=ptry(:) 
      end if 
      amotry=ytry 
    END FUNCTION amotry 
  END SUBROUTINE amoeba 

!      
!!-------------------------------------------------------------------------------------------------------
!! cubic spline routines for respresentation of a discrete function
!!
!! handling: first call the... 
!! subroutine DRSPLN: calculates the coefficients of the representation
!!
!! then find the function values by...
!! function DRSPLE: calculates the function value at a specific location
!
!!      SUBROUTINE DRSPLN
!!$PROG DRSPLN
!      SUBROUTINE DRSPLN(I1,I2,X,Y,Q,F,ilength)
!!
!! !$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
!!
!!   SUBROUTINE RSPLN COMPUTES CUBIC SPLINE INTERPOLATION COEFFICIENTS
!!C   FOR Y(X) BETWEEN GRID POINTS I1 AND I2 SAVING THEM IN Q.  THE
!!C   INTERPOLATION IS CONTINUOUS WITH CONTINUOUS FIRST AND SECOND
!!C   DERIVITIVES.  IT AGREES EXACTLY WITH Y AT GRID POINTS AND WITH THE
!!C   THREE POINT FIRST DERIVITIVES AT BOTH END POINTS (I1 AND I2).
!!C   X MUST BE MONOTONIC BUT IF TWO SUCCESSIVE VALUES OF X ARE EQUAL
!!C   A DISCONTINUITY IS ASSUMED AND SEPERATE INTERPOLATION IS DONE ON
!!C   EACH STRICTLY MONOTONIC SEGMENT.  THE ARRAYS MUST BE DIMENSIONED AT
!!C   LEAST - X(I2), Y(I2), Q(3,I2), AND F(3,I2).  F IS WORKING STORAGE
!!C   FOR RSPLN.
!!C                                                     -RPB
!      implicit double precision(A-H,O-Z)
!      integer ilength
!      DOUBLE PRECISION X,Y,Q,F
!      DIMENSION X(ilength),Y(ilength),Q(3,ilength),F(3,ilength),YY(3)
!      EQUIVALENCE (YY(1),Y0)
!      DATA SMALL/1.D-10/,YY/0.D0,0.D0,0.D0/
!      J1=I1+1
!      Y0=0.D0
!!C   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL.
!      IF(I2-I1)13,17,8
! 8    A0=X(J1-1)
!!C   SEARCH FOR DISCONTINUITIES.
!      DO 3 I=J1,I2
!      B0=A0
!      A0=X(I)
!      IF(DABS((A0-B0)/DMAX1(A0,B0)).LT.SMALL) GO TO 4
! 3    CONTINUE
! 17   J1=J1-1
!      J2=I2-2
!      GO TO 5
! 4    J1=J1-1
!      J2=I-3
!!C   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
! 5    IF(J2+1-J1)9,10,11
!!C   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
! 10   J2=J2+2
!      Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
!      DO 15 J=1,3
!      Q(J,J1)=YY(J)
! 15   Q(J,J2)=YY(J)
!      GO TO 12
!!C   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
! 11   A0=0.
!      H=X(J1+1)-X(J1)
!      H2=X(J1+2)-X(J1)
!      Y0=H*H2*(H2-H)
!      H=H*H
!      H2=H2*H2
!!C   CALCULATE DERIVITIVE AT NEAR END.
!      B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
!      B1=B0
!!C   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
!      DO 1 I=J1,J2
!      H=X(I+1)-X(I)
!      Y0=Y(I+1)-Y(I)
!      H2=H*H
!      HA=H-A0
!      H2A=H-2.D0*A0
!      H3A=2.D0*H-3.D0*A0
!      H2B=H2*B0
!      Q(1,I)=H2/HA
!      Q(2,I)=-HA/(H2A*H2)
!      Q(3,I)=-H*H2A/H3A
!      F(1,I)=(Y0-H*B0)/(H*HA)
!      F(2,I)=(H2B-Y0*(2.D0*H-A0))/(H*H2*H2A)
!      F(3,I)=-(H2B-3.D0*Y0*HA)/(H*H3A)
!      A0=Q(3,I)
! 1    B0=F(3,I)
!!C   TAKE CARE OF LAST TWO ROWS.
!      I=J2+1
!      H=X(I+1)-X(I)
!      Y0=Y(I+1)-Y(I)
!      H2=H*H
!      HA=H-A0
!      H2A=H*HA
!      H2B=H2*B0-Y0*(2.D0*H-A0)
!      Q(1,I)=H2/HA
!      F(1,I)=(Y0-H*B0)/H2A
!      HA=X(J2)-X(I+1)
!      Y0=-H*HA*(HA+H)
!      HA=HA*HA
!!C   CALCULATE DERIVITIVE AT FAR END.
!      Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
!      Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.D0*A0))
!      Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)
!!C   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
!      DO 2 J=J1,J2
!      K=I-1
!      Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
!      Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
!      Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
! 2    I=K
!      Q(1,I)=B1
!!C   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
! 9    J2=J2+2
!      DO 14 J=1,3
! 14   Q(J,J2)=YY(J)
!!C   SEE IF THIS DISCONTINUITY IS THE LAST.
! 12   IF(J2-I2)6,13,13
!!C   NO.  GO BACK FOR MORE.
! 6    J1=J2+2
!      IF(J1-I2)8,8,7
!!C   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
! 7    DO 16 J=1,3
! 16   Q(J,I2)=YY(J)
!!C   FINI.
! 13   RETURN
!      END
!
!
!!C       FUNCTION DRSPLE
!!C$PROG DRSPLE
!      DOUBLE PRECISION FUNCTION DRSPLE(I1,I2,X,Y,Q,ilength,S)
!!C
!!C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
!!C
!!C   RSPLE RETURNS THE VALUE OF THE FUNCTION Y(X) EVALUATED AT POINT S
!!C   USING THE CUBIC SPLINE COEFFICIENTS COMPUTED BY RSPLN AND SAVED IN
!!C   Q.  IF S IS OUTSIDE THE INTERVAL (X(I1),X(I2)) RSPLE EXTRAPOLATES
!!C   USING THE FIRST OR LAST INTERPOLATION POLYNOMIAL.  THE ARRAYS MUST
!!C   BE DIMENSIONED AT LEAST - X(I2), Y(I2), AND Q(3,I2).
!!C
!!C                                                     -RPB
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      integer ilength      
!      DOUBLE PRECISION X,Y,Q,S
!      DIMENSION X(ilength),Y(ilength),Q(3,ilength)
!      DATA I/1/
!      II=I2-1
!!C   GUARANTEE I WITHIN BOUNDS.
!      I=MAX0(I,I1)
!      I=MIN0(I,II)
!!C   SEE IF X IS INCREASING OR DECREASING.
!      IF(X(I2)-X(I1))1,2,2
!!C   X IS DECREASING.  CHANGE I AS NECESSARY.
! 1    IF(S-X(I))3,3,4
! 4    I=I-1
!      IF(I-I1)11,6,1
! 3    IF(S-X(I+1))5,6,6
! 5    I=I+1
!      IF(I-II)3,6,7
!!C   X IS INCREASING.  CHANGE I AS NECESSARY.
! 2    IF(S-X(I+1))8,8,9
! 9    I=I+1
!      IF(I-II)2,6,7
! 8    IF(S-X(I))10,6,6
! 10   I=I-1
!      IF(I-I1)11,6,8
! 7    I=II
!      GO TO 6
! 11   I=I1
!!C   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
! 6    H=S-X(I)
!      DRSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))
!      RETURN
!      END
!      
!
