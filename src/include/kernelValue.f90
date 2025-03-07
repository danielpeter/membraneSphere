!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!
!      Daniel Peter
!      (c) 2025
!
!      Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      
!=====================================================================

!-----------------------------------------------------------------------
      function getKernelValue(deltaScatterer,heterogeneous,t_lag,phaseVelocityRef,arrivalTime,vperturbation)
!-----------------------------------------------------------------------
! calculates the sensitivity kernel based upon travel time anomalies 
!
! input:
!   deltaScatterer        - vertex index of delta scatterer
!   heterogeneous      - flag to use heterogeneous phase map correction
!   t_lag                      - time lag
!   phaseVelocityRef  - phase reference 
!   arrivalTime            - arrival time
!   vperturbation         - velocity perturbation
!
! returns: kernel value [1/rad^2]
      use cells;use phaseVelocityMap
      implicit none
      real(WP):: getKernelValue,t_lag,phaseVelocityRef,phasevelocity,vertexCellArea,arrivalTime,vperturbation
      integer:: deltaScatterer
      logical:: heterogeneous
      
      ! convert cell area into radians square
      vertexCellArea=cellAreas(deltaScatterer)/EARTHRADIUS_SQUARED
      
      ! get velocity correction
      if( heterogeneous ) then
          phasevelocity = sqrt(phaseVelocitySquare(deltaScatterer))
      else
          phasevelocity=phaseVelocityRef
      endif
      !print*,'phase velocity:',phasevelocity, deltaScatterer
      
      if( phasevelocity .lt. 0.001) then
        print*,'phase velocity correction too small:',phasevelocity
        print*,'vertex:',deltaScatterer
        call stopProgram( 'abort - getKernelValue   ')
      endif
      
      ! kernel value term (becomes unit-less)
      getKernelValue=(t_lag*phasevelocity)/(arrivalTime*vperturbation*vertexCellArea)

      return
      end
      
      
      
!-----------------------------------------------------------------------
      subroutine fitSphericalHarmonics()
!-----------------------------------------------------------------------
! determine spherical harmonics expansion for kernel values
!
! original routine taken from lapo boschi's 'xyz2gsh.f'-file
!
!
! returns: spherical harmnics expansion       
    !---find expansion of xyz function over generalized (N=2 or N=4) or scalar spherical
    !---harmonics (N=0) by solution of least squares problem. needed to express
    !---maps of azimuthal anisotropy in terms of generalized s.h. coefficients.
      use propagationStartup; use verbosity
      implicit none
      integer:: lmx,nunk,ilata,nx,ny,idsize,ncoef,icount,immax,ndat
      integer:: i,k,il,il0,im,ix,iy,ialpha,irow,icl,ierror,nono
      real:: y,ypp,ytp,ypppp,ytppp,realchisq,truerms,denom,varred,syndat,error
      character:: namein*80,nameout*80,namein2*80,line*128
      parameter(lmx=40)
      parameter(nunk=(lmx-1)*(2*lmx+6))
      parameter(ilata=(nunk*(nunk+1))/2)
      parameter(nx=180,ny=90,idsize=nx*ny)  !increase idsize if spacing < 1 deg
      dimension y((lmx+1)**2),ypp((lmx-1)*(2*lmx+6)),ytp((lmx-1)*(2*lmx+6))  &
                ,ypppp((lmx-3)*(2*lmx+10)),ytppp((lmx-3)*(2*lmx+10))
      double precision:: wkspc(9*(2*lmx+1)),fileLat,fileLon,fileKernel,fileTlagrelative,fileTlag
      real:: dat,rowa,ata,atd,x,gstore,ytemp,a,xlon,xlat,damp
      dimension dat(2*idsize),rowa(nunk),ata(ilata),atd(nunk),x(nunk) &
                  ,gstore(ilata),ytemp(nunk),a(2*idsize,nunk)
      real:: ataout(ilata),atdout(nunk),xmasked(nunk)
      integer:: maskv(nunk),maskata(ilata)
      integer*4:: dampcoef(nunk)

      !print*,"term 0, 2, or 4?"
      !read*,ialpha
      !print*,"input model?"
      !read*,namein
      namein=datadirectory(1:len_trim(datadirectory))//'ttkernel.rot.dat'
      open(10,file=namein,status='old')
      ! skip first two lines
      read(10,*)line
      read(10,*)line
      
      ! initialize
      ialpha=0
      ncoef=(lmx+1)**2
      damp=1.
            
      do i=1,ilata
        ata(i)=0.
      enddo

      !----------------read in all the gridpoints and increment matrices
      i=1
  1   continue
      !-------------------------------scalar SH
      !read(1,*,end=20)xlon,xlat,dat(i)      
      read(10,'(2f8.2,e16.6,2e12.4)',end=20) xlon,xlat,fileKernel,fileTlagrelative,fileTlag    
      !xlon=real(fileLon)
      !xlat=real(fileLat)
      dat(i)=real(fileKernel)
      !print*,'readin:',xlon,xlat,dat(i)
      
      call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
      do k=1,ncoef
        rowa(k)=y(k)
        a(i,k)=rowa(k)
      enddo
      call contribution_ata_alt(rowa,dat(i),ncoef,nunk,ilata,ata,atd)

      if( VERBOSE .and. mod(i,1000).eq.0)print*,i," gridpoints read"
      !--go back to read the next datum
      i=i+1
      goto 1

  20	if( VERBOSE ) print*,i-1,' gridpoints read total'  
      ndat=i-1
      
      close(10)
      

      !---------------------------------------regularization
      k=0
      il0=ialpha
      !----------------------MUST BE CHANGED FOR GENERALIZED SH PARAMETERIZATION
      do il=il0,lmx
         immax=2*il+1
         do im=1,immax
            k=k+1
            dampcoef(k)=float(il*(il+1))
         enddo
      enddo
      k=0
      icount=0
      do ix=1,ncoef
         do iy=1,ix
            k=k+1
            if(iy.eq.ix)then
               icount=icount+1
                !ata(k)=ata(k)+damp*dampcoef(icount)
      !------------NON-L-DEPENDENT DAMPING FOR THE TIME BEING
               ata(k)=ata(k)+damp
            endif
         enddo
      enddo

      !--mask coefficients associated with basis functions that cannot be constrained
      !SHOULD NOT BE NEEDED NOW
        !print*,"mask ata matrix and atd vector..."
        !call atamask(ata,atd,ataout,atdout,ilata,nunk,maskata,maskv,nonzero)

      !--find least squares solution x via single-processor cholesky factorization of ATA
      print*,'ATA',ncoef
      print*,ata(1:10)
      print*
      if( VERBOSE ) print*,"cholesky factorization..."
      call choles(ata,gstore,atd,ytemp,x,ncoef,nono)
        !call choles(ataout,gstore,atdout,ytemp,xmasked,nonzero,nono)

      if(nono.ne.0)then
        if( VERBOSE )print*,"cholesky factorization failed ",nono
        !do j=0,nonzero-1
        !istart=j*(j+1)/2+1
        !j1=j+1
        !istop=j1*(j1+1)/2
        !c	write(*,*)j+1,(ata(i),i=istart,istop)
        !c	pause
        !enddo
        call stopProgram('fitSphericalHarmonics() - error with nono   ')
      endif

      !--unmask solution coefficients that were eliminated by atamask
        !print*,"unmask solution..."
        !k=0
        !do i=1,nunk
          !if(maskv(i).eq.1)then
              !k=k+1
              !x(i)=xmasked(k)
          !elseif(maskv(i).eq.0)then
              !x(i)=0.
          !else
              !print*,"illegal maskv value for coeff. ",i
              !stop
          !endif
        !enddo

      !--see how well the least squares solution fits the gridpoints
      if( VERBOSE ) print*,"computing variance reduction..."
      truerms=0.
      realchisq=0.
      denom=0.
      do irow=1,ndat
        error=0.
        syndat=0.
        do icl=1,ncoef
          syndat=syndat+a(irow,icl)*x(icl)
        enddo
        error=(syndat-dat(irow))**2
        realchisq=error+realchisq
        truerms=truerms+abs(syndat-dat(irow))
        denom=denom+(dat(irow)*dat(irow))
      enddo
      truerms=truerms/float(ndat)
      varred=realchisq/denom
      varred=1.-varred
      if( VERBOSE ) then
        print*,'NUMBER OF DATA:',ndat
        print*,'CHI-SQUARE:',realchisq
        print*,'RMS MISFIT OBTAINED:',truerms
        print*,'VARIANCE REDUCTION:',varred
      endif

      ! check fit
      if( varred .lt. 0.5 ) then
        call stopProgram('fitSphericalHarmonics() - variance reduction is not good enough!    ')
      endif

      !open(12,file=nameout)
      !write(12,*)lmx
      !write(12,*)(x(i),i=1,ncoef)
      !close(12)
      
      ! store spherical harmonics coefficients
      SH_lmx=lmx
      SH_ncoef=ncoef
      allocate(SH_coef(SH_ncoef),stat=ierror)
      if( ierror .ne. 0 ) call stopProgram('fitSphericalHarmonics() - error allocating SH_coef    ')
      SH_coef(1:SH_ncoef)=x(1:SH_ncoef)

      end
      
!-----------------------------------------------------------------------
      subroutine choles(a,g,b,y,x,n,nono)
!-----------------------------------------------------------------------
! single-processor Cholesky factorization
      implicit real*4 (a-h, o-z)
      real*4 a(1),g(1)
      real*4 b(1),y(1),x(1)
!
!        a= row-wise p.d. symm. system  n*(n+1)/2
!        g= cholesky storage
!        b= r.h.s. vector               n
!        y= temp. vector
!        x= answer vector
!        n= system dimension
!        nono .gt. 0 is the level at which p.d. failed
!
!        (a,g) and (b,y,x) may be equivalenced.
!
!----------------------------------------------------------
!-----first compute cholesky decomposition

      nono=0
      
      if(a(1).le.0.) then
        nono=1
        return
      endif
      
      g(1)=sqrt(a(1))
      y(1)=b(1)/g(1)

      do 400 i=2,n
      
        kz=(i*(i-1))/2
        g(kz+1)=a(kz+1)/g(1)
        sg=g(kz+1)**2
        y(i)=b(i)-g(kz+1)*y(1)
      
        if(i.gt.2) then
      
          jmax=i-1
      
          do 200 j=2,jmax
      
            gkz=a(kz+j)
            kj=(j*(j-1))/2
            kmax=j-1
      
            do 100 k=1,kmax
              gkz=gkz-g(kz+k)*g(kj+k)
 100        continue

            g(kz+j)=gkz/g(kj+j)
            y(i)=y(i)-g(kz+j)*y(j)
            sg=sg+g(kz+j)**2
      
 200      continue

        endif
      
        gkz=a(kz+i)-sg
      
        if(gkz.le.0.) then
          nono=i
          return
        endif
      
        g(kz+i)=sqrt(gkz)
        y(i)=y(i)/g(kz+i)
      
 400  continue

      kz=(n*(n-1))/2
      x(n)=y(n)/g(kz+n)
      if(n.le.1) return

!-----
!     compute solution for particular rhs
      
      do 600 k=2,n
      
        i=n+1-k
        x(i)=y(i)
        jmin=i+1
      
        do 500 j=jmin,n
          kj=(j*(j-1))/2
          x(i)=x(i)-g(kj+i)*x(j)
 500    continue

        kz=(i*(i+1))/2
        x(i)=x(i)/g(kz)
      
 600  continue

      return
      end
      
!-----------------------------------------------------------------------
      subroutine contribution_ata_alt(rowa,dat,ncoef,nunk,ilata,ata,atd)
!-----------------------------------------------------------------------
! given a row of A (rowa) and the corresponding datum (dat), add their
! contribution to the matrix A^tA and the vector A^td
      dimension rowa(nunk),ata(ilata),atd(nunk)
      k1=0
      do i1=1,ncoef
        atd(i1)=atd(i1)+rowa(i1)*dat
        do j1=1,i1
          k1=k1+1 
          ata(k1)=ata(k1)+rowa(i1)*rowa(j1)
        enddo
      enddo
      end

      
!-----------------------------------------------------------------------
      subroutine writeToRegularGrid()
!-----------------------------------------------------------------------
! writes values of spherical harmonics on a regular latitude/longitude-grid
!
! original routine taken from lapo boschi's 'sh2xyz.f'-file
!
!
! returns: output to file 'ttkernel.regular.xyz'  
      use propagationStartup; use verbosity
      implicit none
      integer:: lmx,lmax,lmaxuser,ialpha,ncoef,idsize,nx,ny
      double precision:: phasepercent,scale,vphase,lat,lon,Vsource(3),Vreceiver(3),vtmp(3),rot(3,3)
      double precision:: blockdata
      character:: nameout*80,chlmax*3,chialpha*1
      real:: xlon,xlat,x,z,xincr
      integer:: k,igrid
      parameter(lmx=80)  !maximum harmonic expansion 
      real y((lmx+1)**2)    
      double precision wkspc(9*(2*lmx+1))
      real cmod((LMX+1)**2)
      parameter(nx=180,ny=90,idsize=4*nx*ny)

      !print*,"what input model?"
      !read*,fileName
      !print*,"what alpha (0, 1, 2, 3, 4)?"
      !read*,ialpha
      ialpha=0    ! scalar spherical harmonics
      !if(ialpha.lt.0.or.ialpha.gt.4)stop "non-existing term"
      !print*,"what increment for grid?"
      !read*,xincr
      xincr=1.0
      !c--------------------following line useful e.g. if model given as
      !c--------------------relative perturbation and we need percent
      !print*,"scale model by? (1 if no scale)"
      !read*,scale      
  
      !read(1,*)lmax
      lmax=SH_lmx
      print*,"  maximum l=",lmax
      
      lmaxuser=80
      print*,"  interrupt expansion at l=",lmaxuser

      !if(ialpha.eq.0)then
      ncoef=(lmax+1)**2
      !elseif(ialpha.le.2)then
      !  ncoef=(lmax-1)*(2*lmax+6)
      !else
      !  ncoef=(lmax-3)*(2*lmax+10)
      !endif

      !read(1,*)(cmod(i),i=1,ncoef)      
      !close(1)
      cmod(1:ncoef)=SH_coef(1:ncoef)
      
      lmax=lmaxuser		! user can filter-lapo 2004-04-03
      
      ! file output
      write(chlmax,'(i3.3)')lmax
      write(chialpha,"(i1.1)")ialpha
      nameout=datadirectory(1:len_trim(datadirectory))//'ttkernel.regular.xyz'
      open(12,file=nameout)

      igrid=0
      if(lmax.gt.lmx)stop "l too big"
      do xlat=-89.,89,xincr
        do xlon=0.,359.,xincr
          ! console status display
          igrid=igrid+1          
          if(mod(igrid,1000).eq.0)print*,igrid," grid points done"
          
          ! determine value
          call ylmv0(xlat,xlon,lmx,y,wkspc)
          
          ! scalar spherical harmonic function
          z=0.
          ncoef=(lmax+1)**2
          do k=1,ncoef
            z=z+y(k)*cmod(k)
          enddo
          
          ! file output
          write(12,*)xlon,xlat,z
        enddo
      enddo
            
      ! close output-file
      close(12)
      
      ! console output
      if( VERBOSE ) then
        print*,'  total points read in: ',igrid
        print*,'  regular grid stored as: '//nameout
      endif
      end
      

      
!-----------------------------------------------------------------------
      subroutine interpolateGrid()
!-----------------------------------------------------------------------
      use propagationStartup; use cells
      implicit none
      character*128::filename
      integer:: ierror,i,triangleIndex,lat,lon,corners(3),lengths(3)
      real(WP):: kernelvalue,distances(3),xlat,xlon
      
      do lat=-89,89
        do lon=0,359
          ! determine triangle
          xlat=lat
          xlon=lon
          !call getClosestTriangle(xlat,xlon,triangleIndex,distances,corners,lengths)
          
          ! interpolate
          ! .. to do, need kernel values at each grid point
          !call interpolateLinear(distances,corners,value)
        enddo
      enddo
      
      end
      
      
!-----------------------------------------------------------------------
      subroutine getClosestTriangle(lat,lon,triangleIndex,distances,corners,lengths)
!-----------------------------------------------------------------------
! finds the triangle which contains the desired location
!
! input:
!     lat/lon            - desired location (in degrees)
!     triangleIndex - triangle which contains location
!     distances       - great circle distances between location and triangles corners (in radian)
!     corners          - corner indices of triangle
!     lengths           - triangle side lengths (in radian)
! returns: tirangleIndex, distances and corners

      use cells
      implicit none
      real(WP):: lat,lon,totaldistance,distanceMin,vectorA(3),vectorB(3),distances(3)
      real(WP):: distance,lengths(3),tmpdistances(3)
      integer:: triangleIndex,triangleMin,i,j,corners(3),tmpcorners(3)
      
      ! get vector on sphere of location
      call getVector(lat,lon,vectorA(1),vectorA(2),vectorA(3))      
      !print*,'vector:',lat,lon
      !print*,'    ',vectorA(:)
      
      ! determine triangle where the corners are closest
      distanceMin=1.0e+4
      triangleMin=0
      do i=1,numTriangleFaces
        ! distance to triangle's corners
        totaldistance=0.0
        do j=1,3
          vectorB(:) = vertices(cellTriangleFace(i,j),:)
          call greatCircleDistance(vectorA,vectorB,distance)
          totaldistance = totaldistance+distance
          ! store values
          tmpdistances(j) = distance          
          tmpcorners(j) = cellTriangleFace(i,j)
        enddo
        
        ! check for minimum
        if( distanceMin .gt. totaldistance ) then
          distanceMin = totaldistance
          triangleMin = i
          distances(:) = tmpdistances(:)
          corners(:) = tmpcorners(:)
          !print*,'  new distance :',totaldistance,i
          !print*,'      distances:',(distances(j)*180.0/PI,j=1,3)
          !print*,'      corners  :',(corners(j),j=1,3)          
        endif
      enddo
      
      ! set as result
      triangleIndex=triangleMin
      
      ! check result
      if( triangleIndex .lt. 1 .or. triangleIndex .gt. numTriangleFaces ) then
        print*,'triangle not found',lat,lon,triangleIndex
        call stopProgram('no fitting triangle   ')
      endif      

      ! side lengths of triangle
      vectorA(:) = vertices(corners(1),:)
      vectorB(:) = vertices(corners(2),:)
      call greatCircleDistance(vectorA,vectorB,lengths(1))

      vectorA(:) = vertices(corners(1),:)
      vectorB(:) = vertices(corners(3),:)
      call greatCircleDistance(vectorA,vectorB,lengths(2))

      vectorA(:) = vertices(corners(2),:)
      vectorB(:) = vertices(corners(3),:)
      call greatCircleDistance(vectorA,vectorB,lengths(3))
      
      end
      
!-----------------------------------------------------------------------
      subroutine interpolateLinear(distances,lengths,corners,value)
!-----------------------------------------------------------------------
! linear interpolation of displacement at a location inside a given triangle
!
! input:
!     distances   - great circle distances to corners
!     lengths       - triangle side lengths by order (corner indices 1 to 2, 1 to 3 and 2 to 3)
!     corners      - corner indices
!     value         - interpolated displacement
!
! returns: value
      use displacements
      implicit none
      real(WP):: f_1,f_2,f_3,distance_12,distance_13,distance_23
      real(WP):: distances(3),lengths(3),value
      integer:: corners(3)
      
      ! displacements at triangle corners
      f_1 = newdisplacement(corners(1))
      f_2 = newdisplacement(corners(2))      
      f_3 = newdisplacement(corners(3))
      
      ! triangle side lengths
      distance_12 = lengths(1)
      distance_13 = lengths(2)
      distance_23 = lengths(3)      
      
      ! linear interpolation depending on distances to corners
      value = distances(2)*distances(3)/(distance_12*distance_13)*f_1 &
              + distances(1)*distances(3)/(distance_12*distance_23)*f_2 &
              + distances(1)*distances(2)/(distance_13*distance_23)*f_3

      end
      
!----------------------------------------------------------------
      subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!        subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
!
!  evaluates the vectors ypp,ytp, of length lenv=(lmax-1)*(2*lmax+6)
!  and vectors ypppp,ytppp, of length lenv4=(lmax-3)*(2*lmax+10).
!
!  ypp,ytp represent contributions to the symmetric trace-free
!  second rank tensor c.
!
!  ypppp,ytppp represent contributions to the completely symmetric trace-free
!  fourth rank tensor e.
!
!    In spherical coordinates ('t'=theta=colatitiude, 'p'=phi=longitude)
!
!    c  =  ( c_tt  c_tp )
!          ( c_pt  c_pp )
!
!    and  c_pp   = -c_tt   = ypp . coefs
!         c_tp    = c_pt   = ytp . coefs
!
!     Similarly for e:
!     e_pppp =  e_tttt = -e_pptt   = ypppp . coefs
!     e_tppp = -e_tttp             = ytppp . coefs
!       
! 
!     Scalar spherical harmonics are also calculated and
!     placed in y(1) -- y( (lmax+1)**2 )
!------------------------------------------------------------------
! I believe scalar harmonics here are orthogonal but NOT
! orthonormal: a factor sqrt(2) is missing from harmonics with
! nonzero m. see Dahlen and Tromp eq. B.72
!------------------------------------------------------------------
!
!     The companion routine ylmavv4() (q.v.) calculates the
!     contribution to the average value of  c_ll (= l.c.l) and e_llll
!     along the minor arc and the complete great circle,
!     where 'l' represents the tangent to the path:
!         
!           l_t = -cos(az)      l_p = sin(az)
! 
!     and az is the local azimuth of the path.
!
!     Thus:
!            c_ll   =  - c_pp * cos(2*az)  - c_tp * sin(2*az)
!     and  
!            e_llll =  e_pppp * cos(4*az) + e_tppp * sin(4*az)
!----------------------------------------------------------------
!     'd' is d.p. workspace. (notice that it's bigger than in ylmv)
      dimension y((lmax+1)**2)              &
     &         ,ypp((lmax-1)*(2*lmax+6))    &
     &         ,ytp((lmax-1)*(2*lmax+6))      &
     &         ,ypppp((lmax-3)*(2*lmax+10))   &
     &         ,ytppp((lmax-3)*(2*lmax+10)) 
      double precision theta,d(9*(2*lmax+1))  
      complex cfac,dfac
      data radian/57.2957795/,rr4pi/0.28209479/
      theta=(90.-xlat)/radian
      dfac=cexp(cmplx(0.,xlon/radian))
      k=0
      ka=0
      ka4=0
      do l=0,lmax
        if(l.lt.2) then
          call rotmx2(0,l,theta,d,1,2*l+1)
          ind=l
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            ind=ind+1
            y(k)=d(ind)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
            endif
            cfac=cfac*dfac
          enddo
        else if(l.lt.4) then
          call rotmx2(2,l,theta,d,5,2*l+1)
          ind=5*l+3
          indp=ind+2
          indm=indp
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
            endif
            ind=ind+5
            indp=indp+5
            indm=indm-5
            cfac=cfac*dfac
          enddo
        else
          call rotmx2(4,l,theta,d,9,2*l+1)
          ind=9*l+5
          indp=ind+2
          indm=indp
          indp4=ind+4
          indm4=indp4
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            ka4=ka4+1
            ypppp(ka4)=-d(indp4)*real(cfac)
            ytppp(ka4)=-d(indp4)*aimag(cfac)
            ka4=ka4+1
            ypppp(ka4)=+d(indp4)*aimag(cfac)
            ytppp(ka4)=-d(indp4)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*real(cfac)
              ytppp(ka4)=+d(indm4)*aimag(cfac)
              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*aimag(cfac)
              ytppp(ka4)=-d(indm4)*real(cfac)
            endif
            ind=ind+9
            indp=indp+9
            indm=indm-9
            indp4=indp4+9
            indm4=indm4-9
            cfac=cfac*dfac
          enddo
        endif
      enddo
      return
      end

