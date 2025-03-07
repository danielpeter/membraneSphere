
!-----------------------------------------------------------------------
      subroutine correl77(data1,data2,n,ans)
!-----------------------------------------------------------------------
      implicit none
      integer n,NMAX
      double precision data1(n),data2(n),dn2
      double complex ans(n)
      PARAMETER (NMAX=4096) !Maximum anticipated FFTsize.
! USES realft,twofft
! Computes the correlation of two real datasets data1(1:n) and data2(1:n)
! (including any user-supplied zeropadding).
! n MUST be an integer power of two. The answer is returned as the first n points
! in ans stored in wrap-around order, i.e., correlations at increasingly negative lags
! are in ans(n) on down to ans(n/2+1), while correlations at increasingly positive lags
! are in ans(1)(zerolag) on up to ans(n/2).
! Note that ans must be supplied in the calling program with length
! at least 2*n, since it is also used as workingspace. Sign convention of this routine:
! if data1 lags data2, i.e., is shifted to the right of it, then ans will show a peak at positive lags.
      INTEGER i,no2
      double complex fft(NMAX)

      !print *,'data1:',data1(:)

      call twofft77(data1,data2,fft,ans,n) !Transformboth datavectorsatonce.
      dn2 = 1.0d0/n !Normalizationfor inverseFFT.
      no2 = n/2
      do i = 1,no2+1
        ans(i)=fft(i)*conjg(ans(i))*dn2     ! Multiply to find FFT of their correlation.
      enddo

      ans(1)=dcmplx(dreal(ans(1)),dreal(ans(no2+1))) !Pack first and last into one element.

      !do i=1,n
      !  print *,'complex:',i,ans(i)
      !enddo

      !print *,'calling inverse transform ', n
      call drealft77(ans,n,-1) !Inverse transform gives correlation.

      !do i=1,n
      !  print *,'real:',i,ans(i)
      !enddo

      return
      END


!-----------------------------------------------------------------------
      subroutine twofft77(data1,data2,fft1,fft2,n)
!-----------------------------------------------------------------------
      implicit none
      INTEGER n
      double precision data1(n),data2(n)
      double complex fft1(n),fft2(n)
! USESfour1
! Given two real input arrays data1(1:n)and data2(1:n), this routine calls four1and
! returns twocomplexoutputarrays, fft1(1:n)andfft2(1:n), eachof complexlength n
! (i.e., real length 2*n), which contain the discrete Fourier transforms of the
! respective data arrays. n MUST be an integer power of 2.
      INTEGER j,n2
      double complex h1,h2,c1,c2
      c1=dcmplx(0.5,0.0)
      c2=dcmplx(0.0,-0.5)
      do j = 1,n
        fft1(j)=dcmplx(data1(j),data2(j)) !Packthe tworeal arrays intoone complex array.
      enddo
      call dfour177(fft1,n,1) !Transformthecomplexarray.
      fft2(1)=dcmplx(dimag(fft1(1)),0.0)
      fft1(1)=dcmplx(dreal(fft1(1)),0.0)
      n2 = n+2
      do j = 2,n/2+1
        h1 = c1*(fft1(j)+conjg(fft1(n2-j))) !Usesymmetriestoseparatethetwotransforms.
        h2 = c2*(fft1(j)-conjg(fft1(n2-j)))
        fft1(j)=h1 !Ship them out in two complex arrays.
        fft1(n2-j)=conjg(h1)
        fft2(j)=h2
        fft2(n2-j)=conjg(h2)
      enddo
      return
      END



!-----------------------------------------------------------------------
      subroutine drealft77(data,n,isign)
!-----------------------------------------------------------------------
      implicit none
      INTEGER isign,n
      double precision data(n)
! USESfour1
! Calculates the Fourier transform of a set of n real-valued datapoints.
! Replaces this data (which is stored in array data(1:n)) by the positive frequency
! half of its complex Fourier transform. The real-valued first and last components
! of the complex transform are returned as elements data(1) and data(2), respectively.
! n must be a power of 2. This routine also calculates the inverse transform of
! a complex data array if it is the transform of real data (in this case, complex(1:n/2) becomes data(1:n) ).
! (Result in this case must be multiplied by 2/n.)
      INTEGER i,i1,i2,i3,i4,n2p3
      double precision c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      double precision theta,wi,wpi,wpr,wr,wtemp !Doubleprecisionforthetrigonometricrecurrences.
      theta=3.141592653589793d0/dble(n/2)! Initializetherecurrence.
      c1 = 0.5
      if (isign == 1) then
        c2=-0.5
        call dfour177(data,n/2,+1) !Theforwardtransformis here.
      else
        c2 = 0.5 !Otherwisesetupfor aninversetransform.
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr = 1.0d0+wpr
      wi = wpi
      n2p3 = n+3
      do i = 2,n/4 !Casei = 1doneseparatelybelow.
        i1 = 2*i-1
        i2 = i1+1
        i3 = n2p3-i2
        i4 = i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r = c1*(data(i1)+data(i3)) !Thetwoseparatetransformsareseparatedoutof data.
        h1i = c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i = c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i !Heretheyarerecombinedtoformthetruetransformof theoriginal real data.
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp = wr !Therecurrence.
        wr = wr*wpr-wi*wpi+wr
        wi = wi*wpr+wtemp*wpi+wi
      enddo

      if (isign == 1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2) !Squeeze thefirst andlast datatogether toget themall withintheoriginal array.
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call dfour177(data,n/2,-1) !Thisistheinversetransformforthecaseisign=-1.
      endif
      return
      END


!-----------------------------------------------------------------------
      subroutine dfour177(data,nn,isign)
!-----------------------------------------------------------------------
      implicit none
      INTEGER isign,nn
      double precision data(2*nn)
!Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1;
!or replaces data(1:2*nn)by nn times its inverse discrete Fouriertransform, if isign is input as"1.
!datais acomplex array of lengthnnor, equivalently, areal arrayof length2*nn.
!nn MUSTbeaninteger power of 2(this is not checkedfor!).
      INTEGER i,istep,j,m,mmax,n
      double precision tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp !Doubleprecisionforthetrigonometricrecurrences.
      n = 2*nn
      j = 1
      do i = 1,n,2 !Thisisthebit-reversal sectionof theroutine.
        if (j > i) then
          tempr=data(j) !Exchangethetwocomplexnumbers.
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m = nn
 1     if ((m >= 2).and.(j > m)) then
          j = j-m
          m = m/2
          goto 1
        endif
        j = j+m
      enddo

      mmax = 2 !HerebeginstheDanielson-Lanczossectionof theroutine.
 2    if (n > mmax) then! Outer loopexecutedlog2nntimes.
        istep = 2*mmax
        theta=6.28318530717959d0/(isign*mmax) !Initializeforthetrigonometricrecurrence.
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr = 1.d0
        wi = 0.d0
        do m = 1,mmax,2 !Herearethetwonestedinner loops.
          do i = m,n,istep
            j = i+mmax !ThisistheDanielson-Lanczos formula:
            tempr = wr*data(j)-wi*data(j+1)
            tempi = wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
          enddo
          wtemp = wr !Trigonometricrecurrence.
          wr = wr*wpr-wi*wpi+wr
          wi = wi*wpr+wtemp*wpi+wi
        enddo
        mmax = istep
        goto 2 !Not yet done.
      endif !All done.
      return
      end

