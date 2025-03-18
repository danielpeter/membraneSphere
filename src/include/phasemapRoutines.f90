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
      subroutine readGSHPhasemap(fileName,given_as_percent)
!-----------------------------------------------------------------------
! reads in phase map file with format of spherical harmonics and fill the array phaseMap(:)
! for the non-adjoint simulation such that the phase map is rotated to the source/receiver on equator setup
!
! original routine taken from lapo boschi's 'sh2xyz.f'-file
!
! format of the file must be:
!     id1            blockdata1
!     id2            blockdata2
!
! input:
!     fileName      - pixel map file name
!
! returns: phaseMap() filled with file data
      use phaseVelocityMap; use propagationStartup; use cells; use adjointVariables
      use phaseBlockData; use verbosity
      implicit none
      character(len=128),intent(in):: fileName
      logical,intent(in):: given_as_percent
      ! local parameters
      integer:: i,ioerror
      integer:: lmax,lmaxuser,ialpha,ncoef
      real(WP):: scale,vphase,lat,lon
      real(WP):: Vsource(3),Vreceiver(3),vtmp(3),rot(3,3)
      character(len=3):: chlmax
      character(len=1):: chialpha
      real:: xlon,xlat,z,zmin,zmax
      integer:: k,igrid
      integer, parameter:: lmx = 55  !maximum harmonic expansion
      real:: y((lmx+1)**2)
      double precision:: wkspc(9*(2*lmx+1))
      real:: cmod((LMX+1)**2)

      !parameter(nx=180,ny=90,idsize=4*nx*ny)

      !print *,"what input model?"
      !read*,fileName
      !print *,"what alpha (0, 1, 2, 3, 4)?"
      !read*,ialpha
      ialpha = 0    ! scalar spherical harmonics
      !if (ialpha < 0.or.ialpha>4)stop "non-existing term"
      !print *,"what increment for grid?"
      !read*,xincr
      !c--------------------following line useful e.g. if model given as
      !c--------------------relative perturbation and we need percent
      !print *,"scale model by? (1 if no scale)"
      !read*,scale

      ! read phase map
      if (VERBOSE) then
        print *,'  reading gsh model: '
        print *,'    file: ',trim(fileName)
        if (given_as_percent) then
          print *,'    values given in percentage %'
        else
          print *,'    values given in absolute range'
        endif
      endif

      !open file
      open(IIN,file=trim(fileName),status='old',iostat=ioerror)
      if (ioerror /= 0) then
        print *,'Error: could not open file ',trim(fileName)
        print *,'       Please check if file exists...'
        call stopProgram( 'abort - could not open GSHPhasemap    ')
      endif
      read(IIN,*) lmax

      cmod(:) = 0.0
      lmaxuser = gsh_maximum_expansion
      if (lmaxuser > lmax) call stopProgram('error lmax expansion')

      ncoef = (lmaxuser+1)**2
      read(IIN,*) (cmod(i),i=1,ncoef)
      close(IIN)

      if (VERBOSE) then
        print *,"    maximum expansion     : l = ",lmax
        print *,"    interrupt expansion at: l = ",lmaxuser
        if (rotate_frame) print *,"  rotate frame to equator"
      endif

      ! debug check values
      !open(20,file=trim(fileName)//'-readin.check.dat')
      !write(20,*) lmaxuser, lmax, ncoef
      !write(20,*) (cmod(i),i=1,ncoef)
      !close(20)
      !print *,'  check file:',trim(fileName)//'-readin.check.dat'

      ! sets maxiumum expansion
      lmax = lmaxuser

      ! file output
      write(chlmax,'(i3.3)') lmax
      write(chialpha,"(i1.1)") ialpha

      if (lmax > lmx) then
        print *,'Error: l too big:',lmax,lmx
        call stopProgram( "l too big    ")
      endif

      ! debug phase speed values on grid
      !nameout=trim(fileName)//"-"//chlmax//'.spheregrid.xyz'
      !open(22,file=nameout)
      !open(33,file=trim(fileName)//"-"//chlmax//'.spheregrid.relative.xyz')

      ! phase speeds
      open(IOUT,file=trim(datadirectory)//'PhaseMap.percent.dat')
      ! header comment
      write(IOUT,*) '# Phase map - perturbations'
      write(IOUT,*) '# reference velocity = ',cphaseRef,'(km/s)'
      write(IOUT,*) '# format:'
      write(IOUT,*) '#lon #lat #phase-velocity-perturbation (in percent)'

      ! determine rotation matrix from (1,0,0)/... to source/receiver frame
      if (rotate_frame) then
        Vsource(:) = vertices(originSourceVertex,:)
        Vreceiver(:) = vertices(originReceiverVertex,:)
        call getRotationMatrix(Vsource,Vreceiver,rot)
      endif

      igrid = 0
      zmin = 0.0
      zmax = 0.0
      do i = 1,numVertices
        ! vector to vertex
        vtmp(:) = vertices(i,:)

        ! get rotated vector
        if (rotate_frame) then
          call rotateVector(rot,vtmp,vtmp)
        endif

        ! get latitude/longitude of its location on the phase map
        call getSphericalCoordinates(vtmp,lat,lon)

        xlat = real(lat,kind=4)
        xlon = real(lon,kind=4)

        ! status display
        igrid = igrid+1
        !if (mod(igrid,10000)==0) print *,"  grid points ",igrid

        ! determine value
        call ylmv0(xlat,xlon,lmx,y,wkspc)

        z = 0.0
        ! scalar SH
        ncoef = (lmax+1)**2
        do k = 1,ncoef
          z = z + y(k) * cmod(k)
        enddo

        ! set corresponding phase velocity
        scale = cphaseRef
        if (given_as_percent) then
          vphase = z/100.0 * scale + cphaseRef
        else
          vphase = z * scale + cphaseRef
        endif
        phaseMap(i) = vphase

        ! statistics
        if (z < zmin) zmin = z
        if (z > zmax) zmax = z

        ! debug file output
        !write(20,*)xlon,xlat,vphase
        !write(30,*)xlon,xlat,z

        ! file output
        ! format: #lon #lat #phase-velocity-perturbation (in percent)
        if (given_as_percent) then
          write(IOUT,*) xlon,xlat,z
        else
          write(IOUT,*) xlon,xlat,z*100.0
        endif
      enddo
      close(IOUT)

      ! debug close output-file
      !close(22)
      !close(33)

      if (VERBOSE) then
        print *,'    total points readin   : ',igrid
        print *,'    values minimum/maximum: ',zmin,zmax
        print *,'    phase map stored as   : '//trim(datadirectory)//'PhaseMap.percent.dat'
        ! debug
        !print *,'  debug phase maps stored as: '
        !print *,'    ',trim(nameout)
        print *
      endif

      end subroutine


!-----------------------------------------------------------------------
      subroutine readPixelPhasemap(fileName)
!-----------------------------------------------------------------------
! reads in pixel map file and fill the array phaseBlock(:)
! format of the file must be:
!     id1            blockdata1
!     id2            blockdata2
!
! input:
!     fileName      - pixel map file name
!
! returns: phaseBlock filled with file data blockdata1,..
      use phaseBlockData
      implicit none
      character(len=128),intent(in):: fileName
      ! local parameters
      integer:: i,ioerror,id
      real(WP):: blockdata

      ! read phase map
      print *,'  reading pixel phase map:'
      print *,'  file: ',trim(fileName)

      !open file
      open(IIN, file=trim(fileName),status='old',iostat=ioerror)
      if (ioerror /= 0) call stopProgram( 'abort - readPixelPhasemap()  ')

      !get file length
      numBlocks = 0
      ioerror = 0
      do while (ioerror == 0 )
        read(IIN,*,iostat=ioerror) id, blockdata
        numBlocks = numBlocks + 1
      enddo
      numBlocks = numBlocks - 1
      !print *,'file length', numBlocks

      !resize phaseBlock array
      allocate(phaseBlock(numBlocks), stat=ioerror)
      if (ioerror /= 0) then
        print *,'Error: allocating phaseBlock array'
        call stopProgram( 'abort - readPixelPhasemap    ')
      endif

      !fill in values
      rewind(IIN)
      do i = 1, numBlocks
        read(IIN, *, iostat=ioerror) id, phaseBlock(i)
        if (ioerror /= 0) then
          exit
        endif
        !print *,'block',id,phaseBlock(i)
      enddo

      print *,'number of phase blocks: ',numBlocks
      print *
      close(IIN)

      end

!-----------------------------------------------------------------------
      subroutine determineBlock(reflat,reflon,blockindex)
!-----------------------------------------------------------------------
! this routine is taken from file blk2xyz_reg.f
! further information can be asked from Lapo Boschi
!
! input:
!     relat,reflon  - reference latitude/longitude
!     blockindex    - index of blockdata which covers referenced location
!
! returns: blockindex
      use phaseBlockData, only: heterogeneousPixelsize,compatible_refgridXrefgrid,numBlocks
      implicit none
      real,intent(in):: reflat, reflon
      integer,intent(out):: blockindex
      ! local parameters
      !--------this version creates an xyz file with regularly spaced grid
      !--the following parameters need to be changed depending on the model
      integer,parameter:: ifa      = 1
      real,parameter:: westbo     = -10.
      real,parameter:: eastbo     = 45.
      real,parameter:: southbo    = 30.
      real,parameter:: rthnobo    = 20.
      real,parameter:: PI         = 3.1415926535897931
      !--set iswit=1 to make the grid compatible with a refgridXrefgrid grid
      !   0 for L0035.jwkb.lsqr.mod.pcn-.2-.2, or L150.crust.2degree.gsh.modelVector.*.pcn.blk
      !   1 for L150.crust
      integer,parameter:: iswit     = compatible_refgridXrefgrid
      real:: EQ_INCR,refgrid,eq_hi
      integer:: nlatzones,nlatzohi
      integer,allocatable,dimension(:):: nsqrs,nsqtot,nsqrsh,nsqtoth
      integer,dimension(:),allocatable:: inew,iold,inewh,ioldh,iresolu,ire2
      integer:: numto,numhi,kfine,ifine,kireso,kire2,ilat,ilon,icoarse,ifila,ico
      integer:: isqrl,isqre,superISQRE,icoar,icoblo,iche,nprime,iof,jsqre,isqrh,ifilo
      integer:: i,j,k
      real:: colat,theta,deltalon,parall,rmerid
      real:: xlafi,xlofi,xlat,xlon,xlamin,xlamax,xlomin,xlomax

      ! initialize
      EQ_INCR = heterogeneousPixelsize
      refgrid = eq_incr*1.
      eq_hi = EQ_INCR/ifa
      nlatzones = 180./EQ_INCR
      nlatzohi = 180./eq_hi

      allocate(nsqrs(nlatzones),nsqtot(nlatzones+1))
      allocate(nsqrsh(nlatzohi),nsqtoth(nlatzohi+1))

      allocate(inew(50000),iold(50000),inewh(50000),ioldh(50000),iresolu(10000),ire2(50000))

      !--determine nsqrs,nsqrsh
      NUMTO = 0
      numhi = 0
      colat = -eq_incr/2.

      do k = 1,nlatzones
         !print *,'k',k,nlatzones
         colat = colat+eq_incr
         theta = (colat/180.)*PI
         deltalon = eq_incr/(sin(theta))
         NSQRS(k) = (360./deltalon)+1
         if (MOD(NSQRS(K),2) /= 0) nsqrs(k) = nsqrs(k)-1

         !--if requested, correct nsqrs(k) so the grid is compatible to reference grid
         if (iswit == 1) then
          if (360./nsqrs(k) >= refgrid) then
100         if (mod(360./nsqrs(k),refgrid) /= 0) then
              nsqrs(k) = nsqrs(k)+1
              goto 100
            endif
          else if (360./nsqrs(k) < refgrid) then
101         if (mod(refgrid,360./nsqrs(k)) /= 0) then
              nsqrs(k) = nsqrs(k)+1
              goto 101
            endif
          endif
         endif

         do j = 1,ifa
            kfine = ((k-1)*ifa)+j
            nsqrsh(kfine) = NSQRS(k)*ifa
            nsqtoth(kfine) = numhi
            numhi = numhi+nsqrsh(kfine)
         enddo
         NSQTOT(K) = NUMTO
         NUMTO = NUMTO+NSQRS(K)
      enddo
      nsqtot(nlatzones+1) = numto
      nsqtoth(nlatzohi+1) = numhi
      !print *,'numto=',numto,'  numhi=',numhi

      !--determine iresolu,ire2,kireso,kire2
      kireso = 0
      kire2 = 0
      ! old style
      !do parall = rthnobo,southbo,-eq_incr
      ! new style with loop over integers
      do i = int(rthnobo/eq_incr),int(southbo/eq_incr),-1
         parall = i * eq_incr
         ! old style
         !do rmerid = westbo,eastbo,eq_incr
         ! new style with loop over integers
         do j = int(westbo/eq_incr),int(eastbo/eq_incr)
            rmerid = j * eq_incr
            kireso = kireso+1
            ilat = parall*100.+0.5
            ilon = rmerid*100.+0.5

            iresolu(kireso) = superISQRE(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)
            icoarse = iresolu(kireso)

            call coord_range(icoarse,XLAMIN,XLAMAX,XLOMIN,XLOMAX,nsqrs,nsqtot,nlatzones,eq_incr)

            do ifila = 1,ifa
              do ifilo = 1,ifa
                 kire2 = kire2+1
                 xlafi = xlamin+((xlamax-xlamin)/ifa)*(ifila-0.5)
                 xlofi = xlomin+((xlomax-xlomin)/ifa)*(ifilo-0.5)
                 ilat = xlafi*100.+0.5
                 ilon = xlofi*100.+0.5
                 ire2(kire2) = superISQRE(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
              enddo
            enddo
         enddo
      enddo

      !--count blocks used in parameterization
      icoar = 0
      do icoblo = 1,numto
        do iche = 1,kireso
          if (icoblo == iresolu(iche)) then
            icoar = icoar+1
            goto 37
          endif
        enddo
37      continue
      enddo
      ifine = 0
      do icoblo = 1,numhi
        do iche = 1,kire2
          if (icoblo == ire2(iche)) then
            ifine = ifine+1
            goto 38
          endif
        enddo
38      continue
      enddo
      nprime = numto-icoar+ifine

      ! check number of pixels
      if (nprime /= numBlocks) then
        print *,'Error: pixel grid error'
        print *,'       file blocks',numBlocks
        print *,'       parameterization: total',nprime,' blocks'
        call stopProgram('error - determineBlock')
      endif

      !--determine inew,inewh,iold,ioldh
      call correspc(iresolu,10000,kireso,numto,0,inew,50000,iold,50000,ico)
      iof = numto-icoar
      call corresph(ire2,10000,kire2,numhi,iof,inewh,50000,ioldh,50000,ico)

      !---------find coarse block to which location belongs
      xlat = reflat
      xlon = reflon

      ilat = int(xlat*100.)
      ilon = int(xlon*100.)
      !print *,'looking for:',ilat,ilon

      isqrl = isqre(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)
      jsqre = inew(isqrl)

      !-----check whether such block belongs to the high resolution region
      do iche = 1,kireso
        if (isqrl == iresolu(iche)) then
          isqrh = isqre(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
          jsqre = inewh(isqrh)
        endif
      enddo

      ! return this block id
      blockindex = jsqre

      ! free temporary arrays
      deallocate(inew,iold,inewh,ioldh,iresolu,ire2)
      deallocate(nsqrs,nsqtot)
      deallocate(nsqrsh,nsqtoth)

      end

!-----------------------------------------------------------------------
      subroutine coord_range(NSQ,XLAMIN,XLAMAX,XLOMIN,XLOMAX,nsqrs,nsqtot,nlatzones,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE COORDINATE RANGE OF SQUARE NUMBER 'NSQ'
      implicit none
      integer,intent(in):: nlatzones,nsq
      integer,intent(in):: NSQRS(NLATZONES),NSQTOT(NLATZONES+1)
      real,intent(in):: eq_incr
      real,intent(out):: xlamin,xlamax,xlomin,xlomax
      ! local parameters
      integer:: lazone,nnsq
      real:: grsize

      !print *,'rang', nlatzones
      LAZONE = 2
      DO WHILE (NSQ > NSQTOT(LAZONE))
         LAZONE = LAZONE+1
      enddo
      LAZONE = LAZONE-1

      NNSQ = NSQ-NSQTOT(LAZONE)
      XLAMIN = 90.-LAZONE*eq_incr
      XLAMAX = XLAMIN+eq_incr
      GRSIZE = 360./NSQRS(LAZONE)
      XLOMAX = NNSQ*GRSIZE
      XLOMIN = XLOMAX-GRSIZE
      return
      END

!-----------------------------------------------------------------------
      integer function superisqre(LAT,LON,nsqrs,nsqtot,nlatzones,n,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE NUMBER OF THE SQUARE WHERE (LAT,LON) IS
      implicit none
      integer,intent(in):: nlatzones,n,lat,lon
      integer,intent(in):: NSQRS(nlatzones),NSQTOT(nlatzones+1)
      real,intent(in):: eq_incr
      ! local parameters
      real:: incr
      integer:: llon,lazone
      !print *,'superisqre',nlatzones
      incr = eq_incr*100.
      LAZONE = (9000-LAT)/incr+1
      if (LAZONE > nlatzones) LAZONE = nlatzones

      LLON = LON
      if (LLON < 0) LLON = 36000+LLON
      superISQRE = (LLON*NSQRS(LAZONE))/36000+1
      superISQRE = superISQRE+NSQTOT(LAZONE)

      if (superISQRE > n) superISQRE = n
      return
      END

!-----------------------------------------------------------------------
      integer function isqre(LAT,LON,nsqrs,nsqtot,nlatzones,n,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE NUMBER OF THE SQUARE WHERE (LAT,LON) IS
      implicit none
      integer,intent(in):: nlatzones,n,lat,lon
      integer,intent(in):: NSQRS(nlatzones),NSQTOT(nlatzones+1)
      real,intent(in):: eq_incr
      ! local parameters
      real:: incr
      integer:: llon,lazone

      incr = eq_incr*100.
      LAZONE = (9000-LAT)/incr+1
      if (LAZONE > nlatzones) LAZONE = nlatzones

      LLON = LON
      if (LLON < 0) LLON = 36000+LLON

      ISQRE = (LLON*NSQRS(LAZONE))/36000+1
      ISQRE = ISQRE+NSQTOT(LAZONE)

      if (ISQRE > n) ISQRE = n
      return
      END

!-----------------------------------------------------------------------
      subroutine correspc(iresolu,niresolu,kireso,n,ioffset,inew,ninew,iold,niold,ico)
!-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: niresolu,kireso,n,ioffset,ninew,niold
      integer,intent(in):: iresolu(niresolu)
      integer,intent(inout):: inew(ninew),iold(niold)
      integer,intent(out):: ico
      ! local parameters
      integer:: i,k
      !print *,'correspc',n
      ico = 0
      do i = 1,n
         do k = 1,kireso
         if (i == iresolu(k)) then
            ico = ico+1
            inew(i) = -1
            goto 42
         endif
         enddo
         inew(i) = (i+ioffset)-ico

42       iold(inew(i)) = i
      enddo
      return
      end

!-----------------------------------------------------------------------
      subroutine corresph(iresolu,niresolu,kireso,n,ioffset,inew,ninew,iold,niold,ico)
!-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: niresolu,kireso,n,ioffset,ninew,niold
      integer,intent(in):: iresolu(niresolu)
      integer,intent(inout):: inew(ninew),iold(niold)
      integer,intent(out):: ico
      ! local parameters
      integer:: i,k

      !print *,'corresph',n
      ico = 0
      do i = 1,n
         do k = 1,kireso
         if (i == iresolu(k)) then
            inew(i) = (i+ioffset)-ico
            goto 42
         endif
         enddo
         ico = ico+1
         inew(i) = -1

42       iold(inew(i)) = i
      enddo
      return
      end

!---------------------------------------------------------------------------
      subroutine ylmv0(xlat,xlon,lmax,y,d)
!---------------------------------------------------------------------------
      implicit none
      integer,intent(in):: lmax
      real,intent(in):: xlat,xlon
      real,intent(out):: y((lmax+1)**2)
      double precision,intent(out):: d(9*(2*lmax+1))
      ! local parameters
      integer:: l,k,ind,m
      double precision:: theta
      complex:: cfac,dfac
      double precision,parameter:: radian = 57.2957795
      double precision,parameter:: rr4pi  = 0.28209479

      theta = dble((90.d0- xlat)/radian)
      dfac = cexp(cmplx(0.,xlon/radian))
      k = 0
      do l = 0,lmax
         call rotmx2(0,l,theta,d,1,2*l+1)
         ind = l
         cfac=rr4pi*sqrt(2.*l+1)
         do m = 0,l
            k = k+1
            ind = ind+1
            y(k) = d(ind)*real(cfac)
            if (m /= 0) then
               k = k+1
               y(k) = d(ind)*aimag(cfac)
            endif
            cfac = cfac*dfac
         enddo
      enddo
      return
      end

!---------------------------------------------------------------------------
!prog rotmx2
!xref
  subroutine rotmx2(nmax, l, theta, d, id1, id2)
!---------------------------------------------------------------------------
! subroutine to calculate rotation matrix elements
  implicit none
  ! Input/output parameters
  integer, intent(in) :: nmax, l, id1, id2
  real(kind=8), intent(in) :: theta
  real(kind=8), intent(inout) :: d(id1, id2)
  ! local parameters
  integer :: im1ct, im2, im1, in1, in1ct, in2, in2ct, im2ct
  integer :: m1, nm2, nm2p1, nit, m2, mult, i
  integer :: isup, nm, nmp1, lp1, lm1, lp2, nrow, nmaxp1, lmn
  real(kind=8) :: shth, chth, sth, cth, dlogf, dlogs, t1, temp, rmod
  real(kind=8) :: csum, fac, sfac, rm1, sgn
  real(kind=8) :: th
  ! Constants
  real(kind=8), parameter :: big = 1.0d25
  real(kind=8), parameter :: small = 1.0d-25
  real(kind=8), parameter :: dlbig = 25.0d0
  real(kind=8), parameter :: dlsml = -25.0d0
  real(kind=8), parameter :: pi = 3.14159265358979d0
  real(kind=8), parameter :: EPS = 1.0e-8

  ! Check for l = 0 case
  if (l == 0) then
    d(1+nmax, l+1) = 1.0d0
    return
  endif

  ! Set up parameters based on angle
  isup = 1
  th = theta

  if (abs(th-pi) < EPS) th = pi

  if ((th > pi) .or. (th < 0.0d0)) then
    print *, 'Error: illegal arg in rotmx2:', th, pi
    call stopProgram('illegal rotmx2      ')
  endif

  if (th > pi/2.0d0) then
    th = pi - th
    isup = -1
  endif

  nm = 2*l + 1
  nmp1 = nm + 1
  lp1 = l + 1
  lm1 = l - 1
  lp2 = l + 2
  nrow = 2*nmax + 1
  nmaxp1 = nmax + 1
  lmn = l - nmax

  ! Special case for th = 0
  if (th == 0.0d0) then
    do im1ct = 1, nrow
      im1 = im1ct + lmn
      do im2 = lp1, nm
        d(im1ct, im2) = 0.0d0
        if (im1 == im2) d(im1ct, im2) = 1.0d0
      enddo
    enddo
    goto 400
  endif

  ! Zero left hand side of matrix
  d(1:nrow, 1:lp1) = 0.0d0

  ! Set up parameters
  shth = sin(0.5d0 * th)
  chth = cos(0.5d0 * th)
  sth = 2.0d0 * shth * chth
  cth = 2.0d0 * chth * chth - 1.0d0
  dlogf = log10(chth / shth)
  dlogs = log10(shth)

  ! Iterate from last column using 1.0 as starting value
  do im1ct = 1, nrow
    im1 = im1ct + lmn
    m1 = im1 - lp1
    rm1 = real(m1, kind=8)
    nm2 = min(im1-1, nm-im1)
    d(im1ct, nm) = 1.0d0

    if (nm2 /= 0) then
      do nit = 1, nm2
        m2 = l - nit
        im2 = m2 + lp1

        if (m2 /= lm1) then
          t1 = -sqrt(real((im2+1)*(l-m2-1), kind=8)) * d(im1ct, im2+2)
        else
          t1 = 0.0d0
        endif

        d(im1ct, im2) = t1 - (2.0d0/sth) * (cth*real(m2+1, kind=8)-rm1) * d(im1ct, im2+1)
        d(im1ct, im2) = d(im1ct, im2) / sqrt(real(im2*(l-m2), kind=8))

        temp = d(im1ct, im2)
        rmod = abs(temp)

        if (rmod >= big .and. nit /= nm2) then
          d(im1ct, nit+1) = dlbig
          d(im1ct, im2) = d(im1ct, im2) / big
          d(im1ct, im2+1) = d(im1ct, im2+1) / big
        endif
      enddo
    endif
  enddo

  ! Set up normalization for rightmost column
  t1 = real(2*l, kind=8) * dlogs

  if (lmn /= 0) then
    do i = 1, lmn
      m1 = i - l
      t1 = dlogf + 0.5d0 * log10(real(lp1-m1, kind=8)/real(l+m1, kind=8)) + t1
    enddo
  endif

  d(1, 1) = t1

  if (nrow /= 1) then
    do im1ct = 2, nrow
      m1 = im1ct - nmaxp1
      d(im1ct, 1) = dlogf + 0.5d0 * log10(real(l-m1+1, kind=8)/real(l+m1, kind=8)) + d(im1ct-1, 1)
    enddo
  endif

  sgn = -1.0d0
  if (mod(lmn, 2) /= 0) sgn = 1.0d0

  ! Renormalize rows
  do im1ct = 1, nrow
    im1 = im1ct + lmn
    sgn = -sgn
    csum = d(im1ct, 1)
    mult = 1

    do while (abs(csum) >= dlbig)
      mult = mult * 2
      csum = 0.5 * csum
    enddo

    fac = 10.0d0**csum
    sfac = small / fac
    nm2 = min(im1-1, nm-im1)
    nm2p1 = nm2 + 1

    do im2 = 1, nm2p1
      if ((d(im1ct, im2+1) /= 0.0d0) .and. (im2 < nm2)) then
        csum = csum * real(mult, kind=8) + d(im1ct, im2+1)
        mult = 1

        do while (abs(csum) >= dlbig)
          mult = mult * 2
          csum = 0.5d0 * csum
        enddo

        fac = 10.0d0**csum
        sfac = small / fac
      endif

      in2 = nmp1 - im2

      do i = 1, mult
        temp = d(im1ct, in2)
        rmod = abs(temp)

        if (rmod <= sfac) then
          d(im1ct, in2) = 0.0d0
          exit
        endif

        d(im1ct, in2) = d(im1ct, in2) * fac
      enddo

      d(im1ct, in2) = sgn * d(im1ct, in2)
    enddo
  enddo

  ! Fill rest of matrix
400 continue
  if (isup <= 0) then
    sgn = -1.0d0
    if (mod(lmn, 2) /= 0) sgn = 1.0d0

    do im1ct = 1, nrow
      sgn = -sgn
      im1 = im1ct + lmn
      nm2 = min(im1, nmp1-im1)

      do in2 = 1, nm2
        im2 = nmp1 - in2
        d(im1ct, in2) = sgn * d(im1ct, im2)
      enddo
    enddo

    do im1ct = 1, nrow
      im1 = im1ct + lmn
      in1 = nmp1 - im1
      in1ct = in1 - lmn
      sgn = -1.0d0
      nm2 = min(im1, in1)

      do nit = 1, nm2
        sgn = -sgn
        im2 = 1 + nm2 - nit
        in2 = nmp1 - im2
        im2ct = im2 - lmn
        in2ct = in2 - lmn

        d(in1ct, in2) = sgn * d(im1ct, im2)

        if (in2ct <= nrow) then
          d(im2ct, im1) = d(in1ct, in2)
          d(in2ct, in1) = d(im1ct, im2)
        endif
      enddo
    enddo
  else
    do im1ct = 1, nrow
      im1 = im1ct + lmn
      in1 = nmp1 - im1
      in1ct = in1 - lmn
      sgn = -1.0d0
      nm2 = min(im1, in1)

      do nit = 1, nm2
        sgn = -sgn
        im2 = nm - nm2 + nit
        in2 = nmp1 - im2
        im2ct = im2 - lmn
        in2ct = in2 - lmn

        d(in1ct, in2) = sgn * d(im1ct, im2)

        if (im2ct <= nrow) then
          d(im2ct, im1) = d(in1ct, in2)
          d(in2ct, in1) = d(im1ct, im2)
        endif
      enddo
    enddo
  endif

end subroutine rotmx2

!-----------------------------------------------------------------------
      subroutine makeCheckerboard()
!-----------------------------------------------------------------------
! constructs a phaseMap with the spherical harmonics L=20, M=10
!
! returns: phaseMap() filled with file data
      use phaseVelocityMap; use propagationStartup; use cells
      use parallel; use verbosity
      implicit none
      ! local parameters
      integer:: i,l,m,igrid
      character(len=2):: chL,chM
      real(WP):: vphase,xcolat,xlon,z,zmin,zmax,lat_degree,lon_degree
      double precision:: colatitude,longitude,harmonical,Nfactor
      integer,external:: factorial
      !double precision, external:: generalizedLegendre
      double precision, external:: associatedLegendre

      ! checks that only main process is executing this
      if (.not. MAIN_PROCESS) return

      ! file output
      write(chL,'(i2.2)') MAP_DEGREE_L
      write(chM,'(i2.2)') MAP_DEGREE_M
      open(20,file=trim(datadirectory)//'PhaseMap.L'//chL//'-M'//chM//'.dat')
      open(40,file=trim(datadirectory)//'PhaseMap.L'//chL//'-M'//chM//'.percent.dat')

      igrid = 0
      l = MAP_DEGREE_L
      m = MAP_DEGREE_M

      ! gets maximum
      zmax = 0.0
      do i = 1,numVertices
        ! get colat/lon of vertex
        call getSphericalCoord(i,xcolat,xlon)
        colatitude = xcolat
        harmonical = associatedLegendre(l,m,dcos(colatitude))
        if (dabs(harmonical) > zmax) zmax = dabs(harmonical)
      enddo
      Nfactor = 1.0d0/zmax
      if (MAIN_PROCESS .and. VERBOSE) print *,'normalization: ',Nfactor

      ! creates phase map
      zmax = 0.0
      zmin = 0.0
      do i = 1,numVertices
        igrid = igrid+1
        ! get colat/lon of vertex
        call getSphericalCoord(i,xcolat,xlon)
        lat_degree = 90.0 - xcolat*180./PI
        lon_degree = xlon*180./PI

        !calculate the real part of spherical harmonic function
        z = 0.0
        colatitude = xcolat
        longitude = xlon

        !harmonical = generalizedLegendre(l,m,colatitude)
        harmonical = associatedLegendre(l,m,dcos(colatitude))

        ! scale to [-1,1]
        harmonical = harmonical*Nfactor

        !print *,'colat/lat',colatitude,longitude
        !print *,'  harmonical = ',harmonical

        if (m /= 0) then
          harmonical = harmonical*dsin( m*longitude)
        endif
        z = harmonical

        ! scale amplitude
        z = MAP_PERCENT_AMPLITUDE * z

        !print *,'  perturbation = ',z

        ! set corresponding phase velocity
        ! with perturbation given_as_percent
        vphase = z/100.0_WP*cphaseRef + cphaseRef
        phaseMap(i)=vphase

        ! statistics
        if (z < zmin) zmin = z
        if (z > zmax) zmax = z

        ! file output
        write(20,*)lon_degree,lat_degree,phaseMap(i)
        ! given_as_percent
        write(40,*)lon_degree,lat_degree,z
      enddo

      ! close output-file
      close(20)
      close(40)
      if (MAIN_PROCESS .and. VERBOSE) then
        print *,'  total points readin: ',igrid
        print *,'    percent values minimum/maximum: ',zmin,zmax
        print *,'    phase velocity min/max:',zmin/100.0*cphaseRef+cphaseRef, &
                                      zmax/100.0*cphaseRef+cphaseRef
        print *,'  phase maps stored as: '
        print *,'    ',trim(datadirectory)//'PhaseMap.L'//chL//'-M'//chM//'.dat'
      endif

      end subroutine

!-----------------------------------------------------------------------
!      subroutine rotatePhaseMap()
!-----------------------------------------------------------------------
! rotates phase map from frame source/receiver to (1,0,0)/(0,1,0)/(0,0,1)
!
! returns: rotated cell index
!      use propagationStartup;use phaseVelocityMap
!      implicit none
!      ! local parameters
!      double precision,allocatable,dimension(:):: rotPhaseMap
!      integer i,index,ierror
!      double precision Vsource(3),Vreceiver(3),vtmp(3)
!      double precision rot(3,3),lat,lon
!
!      allocate(rotPhaseMap(numVertices),stat=ierror)
!      if (ierror > 0) then
!        print *,'error in allocating rotPhaseMap array'
!        stop 'abort - rotatePhaseMap'
!      endif
!
!      ! determine rotation matrix from source/receiver to (1,0,0)/... frame
!      Vsource(:)=vertices(originSourceVertex,:)
!      Vreceiver(:)=vertices(originReceiverVertex,:)
!      call getInverseRotationMatrix(Vsource,Vreceiver,rot)
!
!      ! create rotated phase map
!      do i=1,numVertices
!        index=i
!        vtmp(:)=vertices(i,:)
!        call rotateVector(rot,vtmp,vtmp)
!        call getSphericalCoordinates(vtmp,lat,lon)
!        call findVertex(lat,lon,index)
!
!        !call getRotatedIndex(index)
!        !print *,'   index:',i,index
!        !print *,'   phase value:',phaseMap(i)
!        !print *,'   rotated value:',rotPhaseMap(index)
!        if (index < 0 .or. index > numVertices) then
!          print *,'   rotated index violates boundaries:',index, i
!          stop 'abort - rotatePhaseMap'
!        endif
!
!        if (rotPhaseMap(index) < 0.000000001) then
!          rotPhaseMap(index)=phaseMap(i)
!        else
!          print *,'phase map rotation: ',i,' already mapped to', index, rotPhaseMap(index)
!          rotPhaseMap(index)=phaseMap(i)
!          !stop
!        endif
!      enddo
!
!      ! copy it over
!      phaseMap=rotPhaseMap
!
!      ! output
!      open(10,file='tmpRotatedPhaseMap.dat')
!      do i=1, numVertices
!        write(10,*) i, phaseMap(i)
!      enddo
!      close(10)
!      end subroutine
