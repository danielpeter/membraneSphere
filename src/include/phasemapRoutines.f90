!=====================================================================
!
!       m e m b r a n e S p h e r e  1 . 0
!       --------------------------------------------------
!
!      Daniel Peter
!      ETH Zurich - Institute of Geophysics
!      (c) ETH July 2006
!
!      Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      
!=====================================================================

!-----------------------------------------------------------------------
      subroutine readGSHPhasemap(fileName)
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
      use phaseBlockData
      implicit none
      integer:: i,ioerror
      character*128:: fileName
      integer:: lmx,lmax,lmaxuser,ialpha,ncoef
      real(WP):: phasepercent,scale,vphase,lat,lon,Vsource(3),Vreceiver(3),vtmp(3),rot(3,3)
      real(WP):: blockdata
      character:: nameout*80,chlmax*3,chialpha*1
      real:: xlon,xlat,x,z
      integer:: k,igrid
      parameter(lmx=55)  !maximum harmonic expansion 
      real:: y((lmx+1)**2)    
      double precision:: wkspc(9*(2*lmx+1))
      real:: cmod((LMX+1)**2)
      
      !parameter(nx=180,ny=90,idsize=4*nx*ny)

      !print*,"what input model?"
      !read*,fileName
      !print*,"what alpha (0, 1, 2, 3, 4)?"
      !read*,ialpha
      ialpha=0    ! scalar spherical harmonics
      !if(ialpha.lt.0.or.ialpha.gt.4)stop "non-existing term"
      !print*,"what increment for grid?"
      !read*,xincr
      !c--------------------following line useful e.g. if model given as
      !c--------------------relative perturbation and we need percent
      !print*,"scale model by? (1 if no scale)"
      !read*,scale
      
      ! read phase map
      fileName = trim(fileName)
      print*,'reading gsh file: '//fileName   
      
      !open file      
      open(10, file= fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) then
        print*,'  could not open file: '//fileName
        call stopProgram( 'abort - readGSHPhasemap    ')
      endif
  
      read(10,*)lmax
      print*,"  maximum l=",lmax
      lmaxuser = heterogeneousPixelsize
      print*,"  interrupt expansion at l=",lmaxuser

      !if(ialpha.eq.0)then
      ncoef=(lmax+1)**2
      !elseif(ialpha.le.2)then
      !  ncoef=(lmax-1)*(2*lmax+6)
      !else
      !  ncoef=(lmax-3)*(2*lmax+10)
      !endif

      read(10,*)(cmod(i),i=1,ncoef)
      close(10)
      
      
      !do k=1,80
      !  if(fileName(k:k).eq.' ')goto2
      !enddo
      !2	continue

      lmax=lmaxuser		! user can filter-lapo 2004-04-03
      
      ! file output
      write(chlmax,'(i3.3)')lmax
      write(chialpha,"(i1.1)")ialpha
      nameout=fileName(1:len_trim(fileName) )//"-"//chlmax//'.spheregrid.xyz'
      open(20,file=nameout)

      igrid=0
      if(lmax.gt.lmx) then
        print*,'l too big:',lmax,lmx
        call stopProgram( "l too big    ")      
      endif
      !do xlat=-89.,89,xincr
      !do xlon=1.,359.,xincr
      
      ! determine rotation matrix from (1,0,0)/... to source/receiver frame
      if( rotate_frame ) then
        Vsource(:)= vertices(originSourceVertex,:)
        Vreceiver(:)= vertices(originReceiverVertex,:)
        call getRotationMatrix(Vsource,Vreceiver,rot)
      endif
      
      do i=1,numVertices
        ! vector to vertex
        vtmp(:)=vertices(i,:)
        
        ! get rotated vector
        if( rotate_frame ) then
          call rotateVector(rot,vtmp,vtmp)
        endif
        
        ! get latitude/longitude of its location on the phase map
        call getSphericalCoordinates(vtmp,lat,lon)
        if( WP .eq. 4 ) then
          xlat=lat
          xlon=lon
        else
          xlat=real(lat)
          xlon=real(lon)          
        endif
        
        ! status display 
        igrid=igrid+1
        !if(mod(igrid,10000).eq.0)print*,"  grid points ",igrid
        
        ! determine value
        call ylmv0(xlat,xlon,lmx,y,wkspc)
        
        z=0.
        !if(ialpha.eq.0)then
        !c-------------------------------scalar SH
        ncoef=(lmax+1)**2
        do k=1,ncoef
          z=z+y(k)*cmod(k)
        enddo
        !elseif(ialpha.eq.1)then
        !c-------------------------------generalized SH for second order tensor
        !  ncoef=(lmax-1)*(2*lmax+6)
        ! do k=1,ncoef
        !  z=z-ypp(k)*cmod(k)
        ! enddo
        !elseif(ialpha.eq.2)then
        ! ncoef=(lmax-1)*(2*lmax+6)
        ! do k=1,ncoef
	      !   z=z-ytp(k)*cmod(k)
        ! enddo
        !elseif(ialpha.eq.3)then
        !c-------------------------------generalized SH for fourth order tensor
        !c	   ncoef=(lmax-3)*(2*lmax+10)
        !   do k=1,ncoef
        !    z=z+ypppp(k)*cmod(k)
        !   enddo
        !elseif(ialpha.eq.4)then
        ! c	   ncoef=(lmax-3)*(2*lmax+10)
        ! do k=1,ncoef
	      !   z=z+ytppp(k)*cmod(k)
        ! enddo
        !else
        !  stop "error in variable ialpha"
        !endif
        
        
        ! set corresponding phase velocity
        scale=cphaseRef
        vphase=z*scale/100.0_WP+cphaseRef
        phaseMap(i)=vphase        
        
        write(20,*)xlon,xlat,vphase
      enddo
      
      ! close output-file
      close(20)
      
      print*,'  total points read in: ',igrid
      print*,'  phase map stored as: '//nameout
      end

      
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
      integer:: i,ioerror,id
      character*128:: fileName
      real(WP):: phasepercent,blockdata
  
      ! read phase map
      fileName = trim(fileName)
      print*,fileName   
      
      !open file      
      open(10, file= fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readPixelPhasemap()  ')

      !get file length 
      numBlocks = 0     
      ioerror = 0
      do while (ioerror .eq. 0 )
        read(10,*,iostat=ioerror) id, blockdata
        numBlocks = numBlocks + 1
      enddo
      numBlocks = numBlocks - 1
      !print*,'file length', numBlocks
      
      !resize phaseBlock array
      allocate(phaseBlock(numBlocks), stat=ioerror)
      if( ioerror .ne. 0) then
        print*,'error allocating phaseBlock array'
        call stopProgram( 'abort - readPixelPhasemap    ')
      endif
      
      !fill in values
      rewind(10)
      do i=1, numBlocks
        read(10, *, iostat=ioerror) id, phaseBlock(i)
        if( ioerror .ne. 0) then
          exit
        endif
        !print*,'block',id,phaseBlock(i)
      enddo    
      
      print*,'number of phase blocks: ',numBlocks
      print*
      close(10)
      
      end
  
!-----------------------------------------------------------------------
      subroutine determineBlock(reflat,reflon,blockindex)
!-----------------------------------------------------------------------
! this routine is taken from file blk2xyz_reg.f
! further information can be asked from Lapo Boschi
!
! input:
!     relat,reflon - reference latitude/longitude
!     blockindex - index of blockdata which covers referenced location
!
! returns: blockindex
      use phaseBlockData
      implicit none
      real:: pixelsize
      !--------this version creates an xyz file with regularly spaced grid
      !--the following parameters need to be changed depending on the model
      integer,parameter:: ifa      = 1
      real,parameter:: westbo     = -10.
      real,parameter:: eastbo     = 45.
      real,parameter:: southbo    = 30.
      real,parameter:: rthnobo    = 20.
      
      !--set iswit=1 to make the grid compatible with a refgridXrefgrid grid
      !   0 for L0035.jwkb.lsqr.mod.pcn-.2-.2, or L150.crust.2degree.gsh.modelVector.*.pcn.blk
      !   1 for L150.crust
      integer,parameter::iswit     = compatible_refgridXrefgrid

      real:: EQ_INCR,refgrid,eq_hi
      integer:: nlatzones,nlatzohi
      integer,allocatable,dimension(:):: nsqrs,nsqtot,nsqrsh,nsqtoth
      integer:: blockindex
      real:: reflat, reflon      
      integer:: iresolu(10000),ire2(50000)
      integer:: inew(50000),iold(50000),inewh(50000),ioldh(50000)
      character:: filename*80,dummy*80,outpo*80,colorname*30
      integer:: irgb(100,3),rgb(3)
      integer:: numto,numhi,k,j,kfine,ifine,kireso,kire2,ilat,ilon,icoarse,n,ifila,ico
      integer:: isqrl,isqre,superISQRE,icoar,icoblo,iche,nprime,iof,jsqre,isqrh,ifilo
      real:: colat,theta,deltalon,parall,rmerid,v(100)
      real:: xlafi,xlofi,xlat,xlon,xlamin,xlamax,xlomin,xlomax

      ! initialize
      EQ_INCR = heterogeneousPixelsize
      refgrid = eq_incr*1.  
      eq_hi = EQ_INCR/ifa
      nlatzones = 180./EQ_INCR
      nlatzohi = 180./eq_hi
      allocate(nsqrs(nlatzones),nsqtot(nlatzones+1))
      allocate(nsqrsh(nlatzohi),nsqtoth(nlatzohi+1))

      !--determine nsqrs,nsqrsh
      NUMTO=0
      numhi=0
      colat=-eq_incr/2.
      
      do k=1,nlatzones
         !print*,'k',k,nlatzones
         colat=colat+eq_incr
         theta=(colat/180.)*PI
         deltalon=eq_incr/(sin(theta))
         NSQRS(k)=(360./deltalon)+1
         if(MOD(NSQRS(K),2).ne.0)nsqrs(k)=nsqrs(k)-1

         !--if requested, correct nsqrs(k) so the grid is compatible to reference grid
         if(iswit.eq.1)then
          if(360./nsqrs(k).ge.refgrid)then
 100	     if(mod(360./nsqrs(k),refgrid).ne.0)then
            nsqrs(k)=nsqrs(k)+1
            goto 100
           else
           endif
          elseif(360./nsqrs(k).lt.refgrid)then
 101	     if(mod(refgrid,360./nsqrs(k)).ne.0)then
            nsqrs(k)=nsqrs(k)+1
            goto 101
           else
           endif
          endif
         endif

         do j=1,ifa
            kfine=((k-1)*ifa)+j
            nsqrsh(kfine)=NSQRS(k)*ifa
            nsqtoth(kfine)=numhi
            numhi=numhi+nsqrsh(kfine)
         enddo
         NSQTOT(K)=NUMTO
         NUMTO=NUMTO+NSQRS(K)
      ENDDO
      nsqtot(nlatzones+1)=numto
      nsqtoth(nlatzohi+1)=numhi
      !print*,'numto=',numto,'  numhi=',numhi

      !--determine iresolu,ire2,kireso,kire2
      kireso=0
      kire2=0
      do parall=rthnobo,southbo,-eq_incr
         do rmerid=westbo,eastbo,eq_incr
            kireso=kireso+1
            ilat=parall*100.+0.5
            ilon=rmerid*100.+0.5
            iresolu(kireso)=superISQRE(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)
            icoarse=iresolu(kireso)
            call rang(icoarse,XLAMIN,XLAMAX,XLOMIN,XLOMAX,nsqrs,nsqtot,nlatzones,n,eq_incr)
            do ifila=1,ifa
            do ifilo=1,ifa
               kire2=kire2+1
               xlafi=xlamin+((xlamax-xlamin)/ifa)*(ifila-0.5)
               xlofi=xlomin+((xlomax-xlomin)/ifa)*(ifilo-0.5)
               ilat=xlafi*100.+0.5
               ilon=xlofi*100.+0.5
               ire2(kire2)=superISQRE(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
            enddo
            enddo
         enddo
      enddo

      !--count blocks used in parameterization
      icoar=0
      do icoblo=1,numto
         do iche=1,kireso
            IF(icoblo.eq.iresolu(iche))THEN
               icoar=icoar+1
               goto 37
            endif
         enddo
 37   continue
      enddo
      ifine=0
      do icoblo=1,numhi
         do iche=1,kire2
            IF(icoblo.eq.ire2(iche))THEN
               ifine=ifine+1
               goto 38
            endif
         enddo
 38   continue
      enddo
      nprime=numto-icoar+ifine
      
      ! check number of pixels
      if( nprime .ne. numBlocks ) then
        print*,'pixel grid error'
        print*,'    file blocks',numBlocks
        print*,'    parameterization: total',nprime,' blocks'
        call stopProgram('error - determineBlock')
      endif

      !--determine inew,inewh,iold,ioldh
      call correspc(iresolu,10000,kireso,numto,0,inew,50000,iold,50000,ico)
      iof=numto-icoar
      call corresph(ire2,10000,kire2,numhi,iof,inewh,50000,ioldh,50000,ico)

      !---------find coarse block to which location belongs
      xlat = reflat
      xlon = reflon
      
      ilat=int(xlat*100.)
      ilon=int(xlon*100.)
      !print*,'looking for:',ilat,ilon
            
      isqrl=isqre(ilat,ilon,nsqrs,nsqtot,nlatzones,numto,eq_incr)                
      jsqre=inew(isqrl)               
               
      !-----check whether such block belongs to the high resolution region
      do iche=1,kireso
        if(isqrl.eq.iresolu(iche))then
          isqrh=isqre(ilat,ilon,nsqrsh,nsqtoth,nlatzohi,numhi,eq_hi)
          jsqre=inewh(isqrh)
        endif
      enddo       
      
      ! return this block id
      blockindex = jsqre                       
      end    
      
!-----------------------------------------------------------------------
      SUBROUTINE RANG(NSQ,XLAMIN,XLAMAX,XLOMIN,XLOMAX,nsqrs,nsqtot,nlatzones,n,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE COORDINATE RANGE OF SQUARE NUMBER 'NSQ'
      DIMENSION NSQRS(NLATZONES),NSQTOT(NLATZONES+1)
      !print*,'rang', nlatzones
      LAZONE=2
      DO WHILE (NSQ.GT.NSQTOT(LAZONE))
         LAZONE=LAZONE+1
      ENDDO
      LAZONE=LAZONE-1
      NNSQ=NSQ-NSQTOT(LAZONE)
      XLAMIN=90.-LAZONE*eq_incr
      XLAMAX=XLAMIN+eq_incr
      GRSIZE=360./NSQRS(LAZONE)
      XLOMAX=NNSQ*GRSIZE
      XLOMIN=XLOMAX-GRSIZE
      RETURN
      END
!-----------------------------------------------------------------------
      FUNCTION superisqre(LAT,LON,nsqrs,nsqtot,nlatzones,n,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE NUMBER OF THE SQUARE WHERE (LAT,LON) IS
      dimension NSQRS(nlatzones),NSQTOT(nlatzones+1)
      !print*,'superisqre',nlatzones
      incr=eq_incr*100.
      LAZONE=(9000-LAT)/incr+1
      IF(LAZONE.GT.nlatzones)LAZONE=nlatzones
      LLON=LON
      IF(LLON.LT.0)LLON=36000+LLON
      superISQRE=(LLON*NSQRS(LAZONE))/36000+1
      superISQRE=superISQRE+NSQTOT(LAZONE)
      IF(superISQRE.GT.n)ISQRE=n
      RETURN
      END
      
!-----------------------------------------------------------------------
      function isqre(LAT,LON,nsqrs,nsqtot,nlatzones,n,eq_incr)
!-----------------------------------------------------------------------
! FINDS THE NUMBER OF THE SQUARE WHERE (LAT,LON) IS
      DIMENSION NSQRS(nlatzones),NSQTOT(nlatzones+1)
      !print*,'isqre',LAT,LON,nlatzonesn,eq_incr
      incr=eq_incr*100.
      LAZONE=(9000-LAT)/incr+1
      IF(LAZONE.GT.nlatzones)LAZONE=nlatzones
      LLON=LON
      IF(LLON.LT.0)LLON=36000+LLON
      ISQRE=(LLON*NSQRS(LAZONE))/36000+1
      ISQRE=ISQRE+NSQTOT(LAZONE)
      IF(ISQRE.GT.n) ISQRE=n
      !print*,'isqre returns:',isqre
      RETURN
      END
      
!-----------------------------------------------------------------------
      subroutine correspc(iresolu,niresolu,kireso,n,ioffset,inew,ninew,iold,niold,ico)
!-----------------------------------------------------------------------
      dimension iresolu(niresolu),inew(ninew),iold(niold)
      !print*,'correspc',n
      ico=0
      do i=1,n
         do k=1,kireso
         if(i.eq.iresolu(k))then
            ico=ico+1
            inew(i)=-1
            goto 42
         endif
         enddo
         inew(i)=(i+ioffset)-ico
 42	   iold(inew(i))=i
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine corresph(iresolu,niresolu,kireso,n,ioffset,inew,ninew,iold,niold,ico)
!-----------------------------------------------------------------------
      dimension iresolu(niresolu),inew(ninew),iold(niold)
      !print*,'corresph',n
      ico=0
      do i=1,n
         do k=1,kireso
         if(i.eq.iresolu(k))then
            inew(i)=(i+ioffset)-ico
            goto 42
         endif
         enddo
         ico=ico+1
         inew(i)=-1
 42	   iold(inew(i))=i
      enddo
      return
      end



!---------------------------------------------------------------------------
      subroutine ylmv0(xlat,xlon,lmax,y,d)
!---------------------------------------------------------------------------
      integer:: lmax,l
      real:: xlat,xlon,y((lmax+1)**2)
      double precision:: d(9*(2*lmax+1)),theta,radian,rr4pi
      complex:: cfac,dfac
      data radian/57.2957795/,rr4pi/0.28209479/
      
      theta=(90.d0-dble(xlat))/radian
      dfac=cexp(cmplx(0.,xlon/radian))
      k=0
      do l=0,lmax
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
      enddo
      return
      end

!---------------------------------------------------------------------------
!prog rotmx2
!xref
      subroutine rotmx2(nmax,l,theta,d,id1,id2)
!---------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer:: nmax,l,id1,id2
      double precision:: theta, d(id1,id2),big,small,dlbig,dlsml,pi,th
      
!     data big,small,dlbig,dlsml/1.d35,1.d-35,35.d0,-35.d0/
      data big,small,dlbig,dlsml/1.d25,1.d-25,25.d0,-25.d0/
      data pi/3.14159265358979d0/
      
      dfloat(n)=n
      th=theta
      if((th .gt. pi).or.(th .lt. 0.0d0)) then
        print*,'illegal arg in rotmx2:',th
        call stopProgram( 'illegal arg in rotmx2      ')
      endif
      
      if(l.ne.0) goto 350
      d(1+nmax,l+1)=1.d0
      return
 350   isup=1
      if(th.le.pi/2.d0) goto 310
      th=pi-th
      isup=-1
 310   nm=2*l+1
      nmp1=nm+1
      lp1=l+1
      lm1=l-1
      lp2=l+2
      nrow=2*nmax+1
      nmaxp1=nmax+1
      lmn=l-nmax
      if(th.ne.0.d0) goto 320
      do 330 im1ct=1,nrow
      im1=im1ct+lmn
      do 330 im2=lp1,nm
      d(im1ct,im2)=0.d0
      if(im1.eq.im2) d(im1ct,im2)=1.d0
 330   continue
      goto 400
 320   continue
!
!     zero l.h.s. of matrix
!
      do 340 im1=1,nrow
      do 340 im2=1,lp1
 340   d(im1,im2)=0.d0
!
!        set up parameters
!
        shth=dsin(0.5d0*th)
        chth=dcos(0.5d0*th)
        sth=2.d0*shth*chth
        cth=2.d0*chth*chth-1.d0
        dlogf=dlog10(chth/shth)
        dlogs=dlog10(shth)
!
!       iterate from last column using 1. as starting value
!
        do 10 im1ct=1,nrow
          im1=im1ct+lmn
          m1=im1-lp1
          rm1=m1
          nm2=min0(im1-1,nm-im1)
          d(im1ct,nm)=1.d0
          if(nm2.eq.0) goto 10
          do 20 nit=1,nm2
            m2=l-nit
            im2=m2+lp1
            if(m2.ne.lm1) goto 70
            t1=0.d0
            goto 30
 70         t1=-dsqrt(dfloat((im2+1)*(l-m2-1)))*d(im1ct,im2+2)
 30         d(im1ct,im2)=t1-(2.d0/sth)*(cth*dfloat(m2+1)-rm1)*d(im1ct,im2+1)
            d(im1ct,im2)=d(im1ct,im2)/dsqrt(dfloat(im2*(l-m2)))
            temp=d(im1ct,im2)
            rmod=dabs(temp)
            if(rmod.lt.big) goto 20
            if(nit.eq.nm2) goto 20
            d(im1ct,nit+1)=dlbig
            d(im1ct,im2)=d(im1ct,im2)/big
            d(im1ct,im2+1)=d(im1ct,im2+1)/big
 20       continue
 10     continue
!
!        set up normalization for rightmost column
!
        t1=dfloat(2*l)*dlogs
        if(lmn.eq.0) goto 720
        do 710 i=1,lmn
          m1=i-l
          t1=dlogf+0.5d0*dlog10(dfloat(lp1-m1)/dfloat(l+m1))+t1
 710    continue
 720    d(1,1)=t1
        if(nrow.eq.1) goto 730
        do 110 im1ct=2,nrow
          m1=im1ct-nmaxp1
 110      d(im1ct,1)=dlogf+0.5d0*dlog10(dfloat(l-m1+1)/dfloat(l+m1))+d(im1ct-1,1)
 730      sgn=-1.d0
          if((lmn/2)*2.ne.lmn) sgn=1.d0
!
!       renormalize rows
!
          do 120 im1ct=1,nrow
            im1=im1ct+lmn
            sgn=-sgn
            csum=d(im1ct,1)
            mult=1
 520        if(dabs(csum).lt.dlbig) goto 510
            mult=mult*2
            csum=0.5*csum
            goto 520
 510        fac=10.d0**csum
            sfac=small/fac
            nm2=min0(im1-1,nm-im1)
            nm2p1=nm2+1
            do 130 im2=1,nm2p1
              if((d(im1ct,im2+1).eq.0.d0).or.(im2.ge.nm2)) goto 250
              csum=csum*dfloat(mult)+d(im1ct,im2+1)
              mult=1
 220          if(dabs(csum).lt.dlbig) goto 210
              mult=mult*2
              csum=0.5d0*csum
              goto 220
 210          fac=10.d0**csum
              sfac=small/fac
 250          in2=nmp1-im2
              do 270 i=1,mult
                temp=d(im1ct,in2)
                rmod=dabs(temp)
                if(rmod.gt.sfac) goto 260
                d(im1ct,in2)=0.d0
                goto 130
 260            d(im1ct,in2)=d(im1ct,in2)*fac
 270          continue
              d(im1ct,in2)=sgn*d(im1ct,in2)
 130        continue
 120      continue
!
!       fill rest of matrix
!
 400      if(isup.gt.0) goto 410
          sgn=-1.d0
          if((lmn/2)*2.ne.lmn) sgn=1.d0
          do 420 im1ct=1,nrow
            sgn=-sgn
            im1=im1ct+lmn
            nm2=min0(im1,nmp1-im1)
            do 420 in2=1,nm2
              im2=nmp1-in2
 420          d(im1ct,in2)=sgn*d(im1ct,im2)

      do 430 im1ct=1,nrow
        im1=im1ct+lmn
        in1=nmp1-im1
        in1ct=in1-lmn
        sgn=-1.d0
        nm2=min0(im1,in1)
        do 440 nit=1,nm2
          sgn=-sgn
          im2=1+nm2-nit
          in2=nmp1-im2
          im2ct=im2-lmn
          in2ct=in2-lmn
          d(in1ct,in2)=sgn*d(im1ct,im2)
          if(in2ct.gt.nrow) goto 440
          d(im2ct,im1)=d(in1ct,in2)
          d(in2ct,in1)=d(im1ct,im2)
 440    continue
 430  continue
      return
 410  do 450 im1ct=1,nrow
        im1=im1ct+lmn
        in1=nmp1-im1
        in1ct=in1-lmn
        sgn=-1.d0
        nm2=min0(im1,in1)
        do 460 nit=1,nm2
          sgn=-sgn
          im2=nm-nm2+nit
          in2=nmp1-im2
          im2ct=im2-lmn
          in2ct=in2-lmn
          d(in1ct,in2)=sgn*d(im1ct,im2)
          if(im2ct.gt.nrow) goto 460
          d(im2ct,im1)=d(in1ct,in2)
          d(in2ct,in1)=d(im1ct,im2)
 460    continue
 450  continue
      return
      end

