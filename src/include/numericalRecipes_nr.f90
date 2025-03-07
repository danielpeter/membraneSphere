!-----------------------------------------------------------------------
! nrtype.f90
      MODULE nrtype
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
        INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
        INTEGER, PARAMETER :: SP = KIND(1.0)
        INTEGER, PARAMETER :: DP = KIND(1.0D0)
        INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
        INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
        INTEGER, PARAMETER :: LGT = KIND(.true.)
        REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
        REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
        REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
        REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
        REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
        REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
        REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
        REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
        TYPE sprs2_sp
          INTEGER(I4B) :: n,len
          REAL(SP), DIMENSION(:), POINTER :: val
          INTEGER(I4B), DIMENSION(:), POINTER :: irow
          INTEGER(I4B), DIMENSION(:), POINTER :: jcol
        END TYPE sprs2_sp
        TYPE sprs2_dp
          INTEGER(I4B) :: n,len
          REAL(DP), DIMENSION(:), POINTER :: val
          INTEGER(I4B), DIMENSION(:), POINTER :: irow
          INTEGER(I4B), DIMENSION(:), POINTER :: jcol
        END TYPE sprs2_dp
      END MODULE nrtype


!-----------------------------------------------------------------------
! nr.f90
      MODULE nr
        INTERFACE
          subroutine airy(x,ai,bi,aip,bip)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP), INTENT(OUT) :: ai,bi,aip,bip
          end subroutine airy
        END INTERFACE
        INTERFACE
          subroutine amebsa(p,y,pb,yb,ftol,func,iter,temptr)
          USE nrtype
          INTEGER(I4B), INTENT(INOUT) :: iter
          REAL(SP), INTENT(INOUT) :: yb
          REAL(SP), INTENT(IN) :: ftol,temptr
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine amebsa
        END INTERFACE
        INTERFACE
          subroutine amoeba(p,y,ftol,func,iter)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: iter
          REAL(SP), INTENT(IN) :: ftol
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine amoeba
        END INTERFACE
        INTERFACE
          subroutine anneal(x,y,iorder)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          end subroutine anneal
        END INTERFACE
        INTERFACE
          subroutine asolve(b,x,itrnsp)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: b
          REAL(DP), DIMENSION(:), INTENT(OUT) :: x
          INTEGER(I4B), INTENT(IN) :: itrnsp
          end subroutine asolve
        END INTERFACE
        INTERFACE
          subroutine atimes(x,r,itrnsp)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(OUT) :: r
          INTEGER(I4B), INTENT(IN) :: itrnsp
          end subroutine atimes
        END INTERFACE
        INTERFACE
          subroutine avevar(data,ave,var)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data
          REAL(SP), INTENT(OUT) :: ave,var
          end subroutine avevar
        END INTERFACE
        INTERFACE
          subroutine balanc(a)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          end subroutine balanc
        END INTERFACE
        INTERFACE
          subroutine banbks(a,m1,m2,al,indx,b)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: m1,m2
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
          end subroutine banbks
        END INTERFACE
        INTERFACE
          subroutine bandec(a,m1,m2,al,indx,d)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: m1,m2
          INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
          REAL(SP), INTENT(OUT) :: d
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
          end subroutine bandec
        END INTERFACE
        INTERFACE
          subroutine banmul(a,m1,m2,x,b)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: m1,m2
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(OUT) :: b
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          end subroutine banmul
        END INTERFACE
        INTERFACE
          subroutine bcucof(y,y1,y2,y12,d1,d2,c)
          USE nrtype
          REAL(SP), INTENT(IN) :: d1,d2
          REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
          REAL(SP), DIMENSION(4,4), INTENT(OUT) :: c
          end subroutine bcucof
        END INTERFACE
        INTERFACE
          subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy, &
            ansy1,ansy2)
          USE nrtype
          REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
          REAL(SP), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
          REAL(SP), INTENT(OUT) :: ansy,ansy1,ansy2
          end subroutine bcuint
        END INTERFACE
        INTERFACE beschb
          subroutine beschb_s(x,gam1,gam2,gampl,gammi)
          USE nrtype
          REAL(DP), INTENT(IN) :: x
          REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
          end subroutine beschb_s
      !BL
          subroutine beschb_v(x,gam1,gam2,gampl,gammi)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
          end subroutine beschb_v
        END INTERFACE
        INTERFACE bessi
          function bessi_s(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessi_s
          end function bessi_s
      !BL
          function bessi_v(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessi_v
          end function bessi_v
        END INTERFACE
        INTERFACE bessi0
          function bessi0_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessi0_s
          end function bessi0_s
      !BL
          function bessi0_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessi0_v
          end function bessi0_v
        END INTERFACE
        INTERFACE bessi1
          function bessi1_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessi1_s
          end function bessi1_s
      !BL
          function bessi1_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessi1_v
          end function bessi1_v
        END INTERFACE
        INTERFACE
          subroutine bessik(x,xnu,ri,rk,rip,rkp)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,xnu
          REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
          end subroutine bessik
        END INTERFACE
        INTERFACE bessj
          function bessj_s(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessj_s
          end function bessj_s
      !BL
          function bessj_v(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessj_v
          end function bessj_v
        END INTERFACE
        INTERFACE bessj0
          function bessj0_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessj0_s
          end function bessj0_s
      !BL
          function bessj0_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessj0_v
          end function bessj0_v
        END INTERFACE
        INTERFACE bessj1
          function bessj1_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessj1_s
          end function bessj1_s
      !BL
          function bessj1_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessj1_v
          end function bessj1_v
        END INTERFACE
        INTERFACE bessjy
          subroutine bessjy_s(x,xnu,rj,ry,rjp,ryp)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,xnu
          REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
          end subroutine bessjy_s
      !BL
          subroutine bessjy_v(x,xnu,rj,ry,rjp,ryp)
          USE nrtype
          REAL(SP), INTENT(IN) :: xnu
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
          end subroutine bessjy_v
        END INTERFACE
        INTERFACE bessk
          function bessk_s(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessk_s
          end function bessk_s
      !BL
          function bessk_v(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessk_v
          end function bessk_v
        END INTERFACE
        INTERFACE bessk0
          function bessk0_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessk0_s
          end function bessk0_s
      !BL
          function bessk0_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessk0_v
          end function bessk0_v
        END INTERFACE
        INTERFACE bessk1
          function bessk1_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessk1_s
          end function bessk1_s
      !BL
          function bessk1_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessk1_v
          end function bessk1_v
        END INTERFACE
        INTERFACE bessy
          function bessy_s(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessy_s
          end function bessy_s
      !BL
          function bessy_v(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessy_v
          end function bessy_v
        END INTERFACE
        INTERFACE bessy0
          function bessy0_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessy0_s
          end function bessy0_s
      !BL
          function bessy0_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessy0_v
          end function bessy0_v
        END INTERFACE
        INTERFACE bessy1
          function bessy1_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: bessy1_s
          end function bessy1_s
      !BL
          function bessy1_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: bessy1_v
          end function bessy1_v
        END INTERFACE
        INTERFACE beta
          function beta_s(z,w)
          USE nrtype
          REAL(SP), INTENT(IN) :: z,w
          REAL(SP) :: beta_s
          end function beta_s
      !BL
          function beta_v(z,w)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
          REAL(SP), DIMENSION(size(z)) :: beta_v
          end function beta_v
        END INTERFACE
        INTERFACE betacf
          function betacf_s(a,b,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b,x
          REAL(SP) :: betacf_s
          end function betacf_s
      !BL
          function betacf_v(a,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
          REAL(SP), DIMENSION(size(x)) :: betacf_v
          end function betacf_v
        END INTERFACE
        INTERFACE betai
          function betai_s(a,b,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b,x
          REAL(SP) :: betai_s
          end function betai_s
      !BL
          function betai_v(a,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
          REAL(SP), DIMENSION(size(a)) :: betai_v
          end function betai_v
        END INTERFACE
        INTERFACE bico
          function bico_s(n,k)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n,k
          REAL(SP) :: bico_s
          end function bico_s
      !BL
          function bico_v(n,k)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n,k
          REAL(SP), DIMENSION(size(n)) :: bico_v
          end function bico_v
        END INTERFACE
        INTERFACE
          function bnldev(pp,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: pp
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP) :: bnldev
          end function bnldev
        END INTERFACE
        INTERFACE
          function brent(ax,bx,cx,func,tol,xmin)
          USE nrtype
          REAL(SP), INTENT(IN) :: ax,bx,cx,tol
          REAL(SP), INTENT(OUT) :: xmin
          REAL(SP) :: brent
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function brent
        END INTERFACE
        INTERFACE
          subroutine broydn(x,check)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          LOGICAL(LGT), INTENT(OUT) :: check
          end subroutine broydn
        END INTERFACE
        INTERFACE
          subroutine bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
          REAL(SP), INTENT(INOUT) :: x
          REAL(SP), INTENT(IN) :: htry,eps
          REAL(SP), INTENT(OUT) :: hdid,hnext
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine bsstep
        END INTERFACE
        INTERFACE
          subroutine caldat(julian,mm,id,iyyy)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: julian
          INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
          end subroutine caldat
        END INTERFACE
        INTERFACE
          function chder(a,b,c)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(IN) :: c
          REAL(SP), DIMENSION(size(c)) :: chder
          end function chder
        END INTERFACE
        INTERFACE chebev
          function chebev_s(a,b,c,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b,x
          REAL(SP), DIMENSION(:), INTENT(IN) :: c
          REAL(SP) :: chebev_s
          end function chebev_s
      !BL
          function chebev_v(a,b,c,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
          REAL(SP), DIMENSION(size(x)) :: chebev_v
          end function chebev_v
        END INTERFACE
        INTERFACE
          function chebft(a,b,n,func)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(n) :: chebft
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end function chebft
        END INTERFACE
        INTERFACE
          function chebpc(c)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: c
          REAL(SP), DIMENSION(size(c)) :: chebpc
          end function chebpc
        END INTERFACE
        INTERFACE
          function chint(a,b,c)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(IN) :: c
          REAL(SP), DIMENSION(size(c)) :: chint
          end function chint
        END INTERFACE
        INTERFACE
          subroutine choldc(a,p)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: p
          end subroutine choldc
        END INTERFACE
        INTERFACE
          subroutine cholsl(a,p,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          end subroutine cholsl
        END INTERFACE
        INTERFACE
          subroutine chsone(bins,ebins,knstrn,df,chsq,prob)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: knstrn
          REAL(SP), INTENT(OUT) :: df,chsq,prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
          end subroutine chsone
        END INTERFACE
        INTERFACE
          subroutine chstwo(bins1,bins2,knstrn,df,chsq,prob)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: knstrn
          REAL(SP), INTENT(OUT) :: df,chsq,prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: bins1,bins2
          end subroutine chstwo
        END INTERFACE
        INTERFACE
          subroutine cisi(x,ci,si)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP), INTENT(OUT) :: ci,si
          end subroutine cisi
        END INTERFACE
        INTERFACE
          subroutine cntab1(nn,chisq,df,prob,cramrv,ccc)
          USE nrtype
          INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
          REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
          end subroutine cntab1
        END INTERFACE
        INTERFACE
          subroutine cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
          USE nrtype
          INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
          REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
          end subroutine cntab2
        END INTERFACE
        INTERFACE
          function convlv(data,respns,isign)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data
          REAL(SP), DIMENSION(:), INTENT(IN) :: respns
          INTEGER(I4B), INTENT(IN) :: isign
          REAL(SP), DIMENSION(size(data)) :: convlv
          end function convlv
        END INTERFACE
        INTERFACE correl
          function correl_dp(data1,data2)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(DP), DIMENSION(size(data1)) :: correl_dp
          end function correl_dp
        !BL
          function correl_sp(data1,data2)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(SP), DIMENSION(size(data1)) :: correl_sp
          end function correl_sp
        END INTERFACE correl
        INTERFACE
          subroutine cosft1(y)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          end subroutine cosft1
        END INTERFACE
        INTERFACE
          subroutine cosft2(y,isign)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine cosft2
        END INTERFACE
        INTERFACE
          subroutine covsrt(covar,maska)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
          LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
          end subroutine covsrt
        END INTERFACE
        INTERFACE
          subroutine cyclic(a,b,c,alpha,beta,r,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
          REAL(SP), INTENT(IN) :: alpha,beta
          REAL(SP), DIMENSION(:), INTENT(OUT):: x
          end subroutine cyclic
        END INTERFACE
        INTERFACE
          subroutine daub4(a,isign)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine daub4
        END INTERFACE
        INTERFACE dawson
          function dawson_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: dawson_s
          end function dawson_s
      !BL
          function dawson_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: dawson_v
          end function dawson_v
        END INTERFACE
        INTERFACE
          function dbrent(ax,bx,cx,func,dbrent_dfunc,tol,xmin)
          USE nrtype
          REAL(SP), INTENT(IN) :: ax,bx,cx,tol
          REAL(SP), INTENT(OUT) :: xmin
          REAL(SP) :: dbrent
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
      !BL
            function dbrent_dfunc(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: dbrent_dfunc
            end function dbrent_dfunc
          END INTERFACE
          end function dbrent
        END INTERFACE
        INTERFACE
          subroutine ddpoly(c,x,pd)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(IN) :: c
          REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
          end subroutine ddpoly
        END INTERFACE
        INTERFACE
          function decchk(string,ch)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
          CHARACTER(1), INTENT(OUT) :: ch
          LOGICAL(LGT) :: decchk
          end function decchk
        END INTERFACE
        INTERFACE
          subroutine dfpmin(p,gtol,iter,fret,func,dfunc)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: iter
          REAL(SP), INTENT(IN) :: gtol
          REAL(SP), INTENT(OUT) :: fret
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
          INTERFACE
            function func(p)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP) :: func
            end function func
      !BL
            function dfunc(p)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP), DIMENSION(size(p)) :: dfunc
            end function dfunc
          END INTERFACE
          end subroutine dfpmin
        END INTERFACE
        INTERFACE
          function dfridr(func,x,h,err)
          USE nrtype
          REAL(DP), INTENT(IN) :: x,h
          REAL(DP), INTENT(OUT) :: err
          REAL(DP) :: dfridr
          INTERFACE
            function func(x)
            USE nrtype
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
            end function func
          END INTERFACE
          end function dfridr
        END INTERFACE
        INTERFACE
          subroutine dftcor(w,delta,a,b,endpts,corre,corim,corfac)
          USE nrtype
          REAL(SP), INTENT(IN) :: w,delta,a,b
          REAL(SP), INTENT(OUT) :: corre,corim,corfac
          REAL(SP), DIMENSION(:), INTENT(IN) :: endpts
          end subroutine dftcor
        END INTERFACE
        INTERFACE
          subroutine dftint(func,a,b,w,cosint,sinint)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b,w
          REAL(SP), INTENT(OUT) :: cosint,sinint
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end subroutine dftint
        END INTERFACE
        INTERFACE
          subroutine difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: is1,isf,jsf,k,k1,k2
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: s
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: y
          end subroutine difeq
        END INTERFACE
        INTERFACE
          function eclass(lista,listb,n)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: lista,listb
          INTEGER(I4B), INTENT(IN) :: n
          INTEGER(I4B), DIMENSION(n) :: eclass
          end function eclass
        END INTERFACE
        INTERFACE
          function eclazz(equiv,n)
          USE nrtype
          INTERFACE
            function equiv(i,j)
            USE nrtype
            LOGICAL(LGT) :: equiv
            INTEGER(I4B), INTENT(IN) :: i,j
            end function equiv
          END INTERFACE
          INTEGER(I4B), INTENT(IN) :: n
          INTEGER(I4B), DIMENSION(n) :: eclazz
          end function eclazz
        END INTERFACE
        INTERFACE
          function ei(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: ei
          end function ei
        END INTERFACE
        INTERFACE
          subroutine eigsrt(d,v)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
          end subroutine eigsrt
        END INTERFACE
        INTERFACE elle
          function elle_s(phi,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: phi,ak
          REAL(SP) :: elle_s
          end function elle_s
      !BL
          function elle_v(phi,ak)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
          REAL(SP), DIMENSION(size(phi)) :: elle_v
          end function elle_v
        END INTERFACE
        INTERFACE ellf
          function ellf_s(phi,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: phi,ak
          REAL(SP) :: ellf_s
          end function ellf_s
      !BL
          function ellf_v(phi,ak)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
          REAL(SP), DIMENSION(size(phi)) :: ellf_v
          end function ellf_v
        END INTERFACE
        INTERFACE ellpi
          function ellpi_s(phi,en,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: phi,en,ak
          REAL(SP) :: ellpi_s
          end function ellpi_s
      !BL
          function ellpi_v(phi,en,ak)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
          REAL(SP), DIMENSION(size(phi)) :: ellpi_v
          end function ellpi_v
        END INTERFACE
        INTERFACE
          subroutine elmhes(a)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          end subroutine elmhes
        END INTERFACE
        INTERFACE erf
          function erf_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: erf_s
          end function erf_s
      !BL
          function erf_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: erf_v
          end function erf_v
        END INTERFACE
        INTERFACE erfc
          function erfc_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: erfc_s
          end function erfc_s
      !BL
          function erfc_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: erfc_v
          end function erfc_v
        END INTERFACE
        INTERFACE erfcc
          function erfcc_s(x)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: erfcc_s
          end function erfcc_s
      !BL
          function erfcc_v(x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: erfcc_v
          end function erfcc_v
        END INTERFACE
        INTERFACE
          subroutine eulsum(sum,term,jterm)
          USE nrtype
          REAL(SP), INTENT(INOUT) :: sum
          REAL(SP), INTENT(IN) :: term
          INTEGER(I4B), INTENT(IN) :: jterm
          end subroutine eulsum
        END INTERFACE
        INTERFACE
          function evlmem(fdt,d,xms)
          USE nrtype
          REAL(SP), INTENT(IN) :: fdt,xms
          REAL(SP), DIMENSION(:), INTENT(IN) :: d
          REAL(SP) :: evlmem
          end function evlmem
        END INTERFACE
        INTERFACE expdev
          subroutine expdev_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine expdev_s
      !BL
          subroutine expdev_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine expdev_v
        END INTERFACE
        INTERFACE
          function expint(n,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: expint
          end function expint
        END INTERFACE
        INTERFACE factln
          function factln_s(n)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP) :: factln_s
          end function factln_s
      !BL
          function factln_v(n)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
          REAL(SP), DIMENSION(size(n)) :: factln_v
          end function factln_v
        END INTERFACE
        INTERFACE factrl
          function factrl_s(n)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP) :: factrl_s
          end function factrl_s
      !BL
          function factrl_v(n)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
          REAL(SP), DIMENSION(size(n)) :: factrl_v
          end function factrl_v
        END INTERFACE
        INTERFACE
          subroutine fasper(x,y,ofac,hifac,px,py,jmax,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), INTENT(IN) :: ofac,hifac
          INTEGER(I4B), INTENT(OUT) :: jmax
          REAL(SP), INTENT(OUT) :: prob
          REAL(SP), DIMENSION(:), POINTER :: px,py
          end subroutine fasper
        END INTERFACE
        INTERFACE
          subroutine fdjac(x,fvec,df)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
          end subroutine fdjac
        END INTERFACE
        INTERFACE
          subroutine fgauss(x,a,y,dyda)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: y
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
          end subroutine fgauss
        END INTERFACE
        INTERFACE
          subroutine fit(x,y,a,b,siga,sigb,chi2,q,sig)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
          REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
          end subroutine fit
        END INTERFACE
        INTERFACE
          subroutine fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
          REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
          end subroutine fitexy
        END INTERFACE
        INTERFACE
          subroutine fixrts(d)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
          end subroutine fixrts
        END INTERFACE
        INTERFACE
          function fleg(x,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(n) :: fleg
          end function fleg
        END INTERFACE
        INTERFACE
          subroutine flmoon(n,nph,jd,frac)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n,nph
          INTEGER(I4B), INTENT(OUT) :: jd
          REAL(SP), INTENT(OUT) :: frac
          end subroutine flmoon
        END INTERFACE
        INTERFACE four1
          subroutine four1_dp(data,isign)
          USE nrtype
          COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four1_dp
      !BL
          subroutine four1_sp(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four1_sp
        END INTERFACE
        INTERFACE
          subroutine four1_alt(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four1_alt
        END INTERFACE
        INTERFACE
          subroutine four1_gather(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four1_gather
        END INTERFACE
        INTERFACE
          subroutine four2(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
          INTEGER(I4B),INTENT(IN) :: isign
          end subroutine four2
        END INTERFACE
        INTERFACE
          subroutine four2_alt(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four2_alt
        END INTERFACE
        INTERFACE
          subroutine four3(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
          INTEGER(I4B),INTENT(IN) :: isign
          end subroutine four3
        END INTERFACE
        INTERFACE
          subroutine four3_alt(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine four3_alt
        END INTERFACE
        INTERFACE
          subroutine fourcol(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourcol
        END INTERFACE
        INTERFACE
          subroutine fourcol_3d(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourcol_3d
        END INTERFACE
        INTERFACE
          subroutine fourn_gather(data,nn,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourn_gather
        END INTERFACE
        INTERFACE fourrow
          subroutine fourrow_dp(data,isign)
          USE nrtype
          COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourrow_dp
      !BL
          subroutine fourrow_sp(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourrow_sp
        END INTERFACE
        INTERFACE
          subroutine fourrow_3d(data,isign)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine fourrow_3d
        END INTERFACE
        INTERFACE
          function fpoly(x,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(n) :: fpoly
          end function fpoly
        END INTERFACE
        INTERFACE
          subroutine fred2(a,b,t,f,w,g,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
          INTERFACE
            function g(t)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: t
            REAL(SP), DIMENSION(size(t)) :: g
            end function g
      !BL
            function ak(t,s)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(size(t),size(s)) :: ak
            end function ak
          END INTERFACE
          end subroutine fred2
        END INTERFACE
        INTERFACE
          function fredin(x,a,b,t,f,w,g,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
          REAL(SP), DIMENSION(size(x)) :: fredin
          INTERFACE
            function g(t)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: t
            REAL(SP), DIMENSION(size(t)) :: g
            end function g
      !BL
            function ak(t,s)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(size(t),size(s)) :: ak
            end function ak
          END INTERFACE
          end function fredin
        END INTERFACE
        INTERFACE
          subroutine frenel(x,s,c)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP), INTENT(OUT) :: s,c
          end subroutine frenel
        END INTERFACE
        INTERFACE
          subroutine frprmn(p,ftol,iter,fret)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: iter
          REAL(SP), INTENT(IN) :: ftol
          REAL(SP), INTENT(OUT) :: fret
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
          end subroutine frprmn
        END INTERFACE
        INTERFACE
          subroutine ftest(data1,data2,f,prob)
          USE nrtype
          REAL(SP), INTENT(OUT) :: f,prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          end subroutine ftest
        END INTERFACE
        INTERFACE
          function gamdev(ia)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: ia
          REAL(SP) :: gamdev
          end function gamdev
        END INTERFACE
        INTERFACE gammln
          function gammln_s(xx)
          USE nrtype
          REAL(SP), INTENT(IN) :: xx
          REAL(SP) :: gammln_s
          end function gammln_s
      !BL
          function gammln_v(xx)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xx
          REAL(SP), DIMENSION(size(xx)) :: gammln_v
          end function gammln_v
        END INTERFACE
        INTERFACE gammp
          function gammp_s(a,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,x
          REAL(SP) :: gammp_s
          end function gammp_s
      !BL
          function gammp_v(a,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
          REAL(SP), DIMENSION(size(a)) :: gammp_v
          end function gammp_v
        END INTERFACE
        INTERFACE gammq
          function gammq_s(a,x)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,x
          REAL(SP) :: gammq_s
          end function gammq_s
      !BL
          function gammq_v(a,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
          REAL(SP), DIMENSION(size(a)) :: gammq_v
          end function gammq_v
        END INTERFACE
        INTERFACE gasdev
          subroutine gasdev_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine gasdev_s
      !BL
          subroutine gasdev_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine gasdev_v
        END INTERFACE
        INTERFACE
          subroutine gaucof(a,b,amu0,x,w)
          USE nrtype
          REAL(SP), INTENT(IN) :: amu0
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
          end subroutine gaucof
        END INTERFACE
        INTERFACE
          subroutine gauher(x,w)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
          end subroutine gauher
        END INTERFACE
        INTERFACE
          subroutine gaujac(x,w,alf,bet)
          USE nrtype
          REAL(SP), INTENT(IN) :: alf,bet
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
          end subroutine gaujac
        END INTERFACE
        INTERFACE
          subroutine gaulag(x,w,alf)
          USE nrtype
          REAL(SP), INTENT(IN) :: alf
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
          end subroutine gaulag
        END INTERFACE
        INTERFACE
          subroutine gauleg(x1,x2,x,w)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
          end subroutine gauleg
        END INTERFACE
        INTERFACE
          subroutine gaussj(a,b)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
          end subroutine gaussj
        END INTERFACE
        INTERFACE gcf
          function gcf_s(a,x,gln)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,x
          REAL(SP), OPTIONAL, INTENT(OUT) :: gln
          REAL(SP) :: gcf_s
          end function gcf_s
      !BL
          function gcf_v(a,x,gln)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
          REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
          REAL(SP), DIMENSION(size(a)) :: gcf_v
          end function gcf_v
        END INTERFACE
        INTERFACE
          function golden(ax,bx,cx,func,tol,xmin)
          USE nrtype
          REAL(SP), INTENT(IN) :: ax,bx,cx,tol
          REAL(SP), INTENT(OUT) :: xmin
          REAL(SP) :: golden
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function golden
        END INTERFACE
        INTERFACE gser
          function gser_s(a,x,gln)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,x
          REAL(SP), OPTIONAL, INTENT(OUT) :: gln
          REAL(SP) :: gser_s
          end function gser_s
      !BL
          function gser_v(a,x,gln)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
          REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
          REAL(SP), DIMENSION(size(a)) :: gser_v
          end function gser_v
        END INTERFACE
        INTERFACE
          subroutine hqr(a,wr,wi)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          end subroutine hqr
        END INTERFACE
        INTERFACE
          subroutine hunt(xx,x,jlo)
          USE nrtype
          INTEGER(I4B), INTENT(INOUT) :: jlo
          REAL(SP), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(IN) :: xx
          end subroutine hunt
        END INTERFACE
        INTERFACE
          subroutine hypdrv(s,ry,rdyds)
          USE nrtype
          REAL(SP), INTENT(IN) :: s
          REAL(SP), DIMENSION(:), INTENT(IN) :: ry
          REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
          end subroutine hypdrv
        END INTERFACE
        INTERFACE
          function hypgeo(a,b,c,z)
          USE nrtype
          COMPLEX(SPC), INTENT(IN) :: a,b,c,z
          COMPLEX(SPC) :: hypgeo
          end function hypgeo
        END INTERFACE
        INTERFACE
          subroutine hypser(a,b,c,z,series,deriv)
          USE nrtype
          COMPLEX(SPC), INTENT(IN) :: a,b,c,z
          COMPLEX(SPC), INTENT(OUT) :: series,deriv
          end subroutine hypser
        END INTERFACE
        INTERFACE
          function icrc(crc,buf,jinit,jrev)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
          INTEGER(I2B), INTENT(IN) :: crc,jinit
          INTEGER(I4B), INTENT(IN) :: jrev
          INTEGER(I2B) :: icrc
          end function icrc
        END INTERFACE
        INTERFACE
          function igray(n,is)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n,is
          INTEGER(I4B) :: igray
          end function igray
        END INTERFACE
        INTERFACE
          RECURSIVE subroutine index_bypack(arr,index,partial)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: arr
          INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: index
          INTEGER, OPTIONAL, INTENT(IN) :: partial
          end subroutine index_bypack
        END INTERFACE
        INTERFACE indexx
          subroutine indexx_sp(arr,index)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: arr
          INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
          end subroutine indexx_sp
          subroutine indexx_i4b(iarr,index)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
          INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
          end subroutine indexx_i4b
        END INTERFACE
        INTERFACE
          function interp(uc)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
          REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
          end function interp
        END INTERFACE
        INTERFACE
          function rank(indx)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
          INTEGER(I4B), DIMENSION(size(indx)) :: rank
          end function rank
        END INTERFACE
        INTERFACE
          function irbit1(iseed)
          USE nrtype
          INTEGER(I4B), INTENT(INOUT) :: iseed
          INTEGER(I4B) :: irbit1
          end function irbit1
        END INTERFACE
        INTERFACE
          function irbit2(iseed)
          USE nrtype
          INTEGER(I4B), INTENT(INOUT) :: iseed
          INTEGER(I4B) :: irbit2
          end function irbit2
        END INTERFACE
        INTERFACE
          subroutine jacobi(a,d,v,nrot)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: nrot
          REAL(SP), DIMENSION(:), INTENT(OUT) :: d
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
          end subroutine jacobi
        END INTERFACE
        INTERFACE
          subroutine jacobn(x,y,dfdx,dfdy)
          USE nrtype
          REAL(SP), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(IN) :: y
          REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
          end subroutine jacobn
        END INTERFACE
        INTERFACE
          function julday(mm,id,iyyy)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: mm,id,iyyy
          INTEGER(I4B) :: julday
          end function julday
        END INTERFACE
        INTERFACE
          subroutine kendl1(data1,data2,tau,z,prob)
          USE nrtype
          REAL(SP), INTENT(OUT) :: tau,z,prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          end subroutine kendl1
        END INTERFACE
        INTERFACE
          subroutine kendl2(tab,tau,z,prob)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
          REAL(SP), INTENT(OUT) :: tau,z,prob
          end subroutine kendl2
        END INTERFACE
        INTERFACE
          function kermom(y,m)
          USE nrtype
          REAL(DP), INTENT(IN) :: y
          INTEGER(I4B), INTENT(IN) :: m
          REAL(DP), DIMENSION(m) :: kermom
          end function kermom
        END INTERFACE
        INTERFACE
          subroutine ks2d1s(x1,y1,quadvl,d1,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
          REAL(SP), INTENT(OUT) :: d1,prob
          INTERFACE
            subroutine quadvl(x,y,fa,fb,fc,fd)
            USE nrtype
            REAL(SP), INTENT(IN) :: x,y
            REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
            end subroutine quadvl
          END INTERFACE
          end subroutine ks2d1s
        END INTERFACE
        INTERFACE
          subroutine ks2d2s(x1,y1,x2,y2,d,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
          REAL(SP), INTENT(OUT) :: d,prob
          end subroutine ks2d2s
        END INTERFACE
        INTERFACE
          subroutine ksone(data,func,d,prob)
          USE nrtype
          REAL(SP), INTENT(OUT) :: d,prob
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end subroutine ksone
        END INTERFACE
        INTERFACE
          subroutine kstwo(data1,data2,d,prob)
          USE nrtype
          REAL(SP), INTENT(OUT) :: d,prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          end subroutine kstwo
        END INTERFACE
        INTERFACE
          subroutine laguer(a,x,its)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: its
          COMPLEX(SPC), INTENT(INOUT) :: x
          COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
          end subroutine laguer
        END INTERFACE
        INTERFACE
          subroutine lfit(x,y,sig,a,maska,covar,chisq,funcs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
          REAL(SP), INTENT(OUT) :: chisq
          INTERFACE
            subroutine funcs(x,arr)
            USE nrtype
            REAL(SP),INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
            end subroutine funcs
          END INTERFACE
          end subroutine lfit
        END INTERFACE
        INTERFACE
          subroutine linbcg(b,x,itol,tol,itmax,iter,err)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: b
          REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
          INTEGER(I4B), INTENT(IN) :: itol,itmax
          REAL(DP), INTENT(IN) :: tol
          INTEGER(I4B), INTENT(OUT) :: iter
          REAL(DP), INTENT(OUT) :: err
          end subroutine linbcg
        END INTERFACE
        INTERFACE
          subroutine linmin(p,xi,fret)
          USE nrtype
          REAL(SP), INTENT(OUT) :: fret
          REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
          end subroutine linmin
        END INTERFACE
        INTERFACE
          subroutine lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
          REAL(SP), INTENT(IN) :: fold,stpmax
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x
          REAL(SP), INTENT(OUT) :: f
          LOGICAL(LGT), INTENT(OUT) :: check
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP) :: func
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            end function func
          END INTERFACE
          end subroutine lnsrch
        END INTERFACE
        INTERFACE
          function locate(xx,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xx
          REAL(SP), INTENT(IN) :: x
          INTEGER(I4B) :: locate
          end function locate
        END INTERFACE
        INTERFACE
          function lop(u)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
          REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
          end function lop
        END INTERFACE
        INTERFACE
          subroutine lubksb(a,indx,b)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
          end subroutine lubksb
        END INTERFACE
        INTERFACE
          subroutine ludcmp(a,indx,d)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
          REAL(SP), INTENT(OUT) :: d
          end subroutine ludcmp
        END INTERFACE
        INTERFACE
          subroutine machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp, &
            maxexp,eps,epsneg,xmin,xmax)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp, &
            minexp,negep,ngrd
          REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
          end subroutine machar
        END INTERFACE
        INTERFACE
          subroutine medfit(x,y,a,b,abdev)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), INTENT(OUT) :: a,b,abdev
          end subroutine medfit
        END INTERFACE
        INTERFACE
          subroutine memcof(data,xms,d)
          USE nrtype
          REAL(SP), INTENT(OUT) :: xms
          REAL(SP), DIMENSION(:), INTENT(IN) :: data
          REAL(SP), DIMENSION(:), INTENT(OUT) :: d
          end subroutine memcof
        END INTERFACE
        INTERFACE
          subroutine mgfas(u,maxcyc)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
          INTEGER(I4B), INTENT(IN) :: maxcyc
          end subroutine mgfas
        END INTERFACE
        INTERFACE
          subroutine mglin(u,ncycle)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
          INTEGER(I4B), INTENT(IN) :: ncycle
          end subroutine mglin
        END INTERFACE
        INTERFACE
          subroutine midexp(funk,aa,bb,s,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: aa,bb
          REAL(SP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function funk(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
            end function funk
          END INTERFACE
          end subroutine midexp
        END INTERFACE
        INTERFACE
          subroutine midinf(funk,aa,bb,s,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: aa,bb
          REAL(SP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function funk(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
            end function funk
          END INTERFACE
          end subroutine midinf
        END INTERFACE
        INTERFACE
          subroutine midpnt(func,a,b,s,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end subroutine midpnt
        END INTERFACE
        INTERFACE
          subroutine midsql(funk,aa,bb,s,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: aa,bb
          REAL(SP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function funk(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
            end function funk
          END INTERFACE
          end subroutine midsql
        END INTERFACE
        INTERFACE
          subroutine midsqu(funk,aa,bb,s,n)
          USE nrtype
          REAL(SP), INTENT(IN) :: aa,bb
          REAL(SP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function funk(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
            end function funk
          END INTERFACE
          end subroutine midsqu
        END INTERFACE
        INTERFACE
          RECURSIVE subroutine miser(func,regn,ndim,npts,dith,ave,var)
          USE nrtype
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP) :: func
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            end function func
          END INTERFACE
          REAL(SP), DIMENSION(:), INTENT(IN) :: regn
          INTEGER(I4B), INTENT(IN) :: ndim,npts
          REAL(SP), INTENT(IN) :: dith
          REAL(SP), INTENT(OUT) :: ave,var
          end subroutine miser
        END INTERFACE
        INTERFACE
          subroutine mmid(y,dydx,xs,htot,nstep,yout,derivs)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: nstep
          REAL(SP), INTENT(IN) :: xs,htot
          REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine mmid
        END INTERFACE
        INTERFACE
          subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
          USE nrtype
          REAL(SP), INTENT(INOUT) :: ax,bx
          REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine mnbrak
        END INTERFACE
        INTERFACE
          subroutine mnewt(ntrial,x,tolx,tolf,usrfun)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: ntrial
          REAL(SP), INTENT(IN) :: tolx,tolf
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          INTERFACE
            subroutine usrfun(x,fvec,fjac)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
            end subroutine usrfun
          END INTERFACE
          end subroutine mnewt
        END INTERFACE
        INTERFACE
          subroutine moment(data,ave,adev,sdev,var,skew,curt)
          USE nrtype
          REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
          REAL(SP), DIMENSION(:), INTENT(IN) :: data
          end subroutine moment
        END INTERFACE
        INTERFACE
          subroutine mp2dfr(a,s,n,m)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          INTEGER(I4B), INTENT(OUT) :: m
          CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
          CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
          end subroutine mp2dfr
        END INTERFACE
        INTERFACE
          subroutine mpdiv(q,r,u,v,n,m)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
          INTEGER(I4B), INTENT(IN) :: n,m
          end subroutine mpdiv
        END INTERFACE
        INTERFACE
          subroutine mpinv(u,v,n,m)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
          INTEGER(I4B), INTENT(IN) :: n,m
          end subroutine mpinv
        END INTERFACE
        INTERFACE
          subroutine mpmul(w,u,v,n,m)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
          CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
          INTEGER(I4B), INTENT(IN) :: n,m
          end subroutine mpmul
        END INTERFACE
        INTERFACE
          subroutine mppi(n)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          end subroutine mppi
        END INTERFACE
        INTERFACE
          subroutine mprove(a,alud,indx,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
          REAL(SP), DIMENSION(:), INTENT(IN) :: b
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          end subroutine mprove
        END INTERFACE
        INTERFACE
          subroutine mpsqrt(w,u,v,n,m)
          USE nrtype
          CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
          CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
          INTEGER(I4B), INTENT(IN) :: n,m
          end subroutine mpsqrt
        END INTERFACE
        INTERFACE
          subroutine mrqcof(x,y,sig,a,maska,alpha,beta,chisq,funcs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,a,sig
          REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
          REAL(SP), INTENT(OUT) :: chisq
          LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
          INTERFACE
            subroutine funcs(x,a,yfit,dyda)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
            REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
            end subroutine funcs
          END INTERFACE
          end subroutine mrqcof
        END INTERFACE
        INTERFACE
          subroutine mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
          REAL(SP), INTENT(OUT) :: chisq
          REAL(SP), INTENT(INOUT) :: alamda
          LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
          INTERFACE
            subroutine funcs(x,a,yfit,dyda)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
            REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
            end subroutine funcs
          END INTERFACE
          end subroutine mrqmin
        END INTERFACE
        INTERFACE
          subroutine newt(x,check)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
          LOGICAL(LGT), INTENT(OUT) :: check
          end subroutine newt
        END INTERFACE
        INTERFACE
          subroutine odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
          REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
      !BL
            subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
            REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
            REAL(SP), INTENT(INOUT) :: x
            REAL(SP), INTENT(IN) :: htry,eps
            REAL(SP), INTENT(OUT) :: hdid,hnext
              INTERFACE
              subroutine derivs(x,y,dydx)
                USE nrtype
                REAL(SP), INTENT(IN) :: x
                REAL(SP), DIMENSION(:), INTENT(IN) :: y
                REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
                end subroutine derivs
              END INTERFACE
            end subroutine rkqs
          END INTERFACE
          end subroutine odeint
        END INTERFACE
        INTERFACE
          subroutine orthog(anu,alpha,beta,a,b)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
          REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
          end subroutine orthog
        END INTERFACE
        INTERFACE
          subroutine pade(cof,resid)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
          REAL(SP), INTENT(OUT) :: resid
          end subroutine pade
        END INTERFACE
        INTERFACE
          function pccheb(d)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: d
          REAL(SP), DIMENSION(size(d)) :: pccheb
          end function pccheb
        END INTERFACE
        INTERFACE
          subroutine pcshft(a,b,d)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
          end subroutine pcshft
        END INTERFACE
        INTERFACE
          subroutine pearsn(x,y,r,prob,z)
          USE nrtype
          REAL(SP), INTENT(OUT) :: r,prob,z
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          end subroutine pearsn
        END INTERFACE
        INTERFACE
          subroutine period(x,y,ofac,hifac,px,py,jmax,prob)
          USE nrtype
          INTEGER(I4B), INTENT(OUT) :: jmax
          REAL(SP), INTENT(IN) :: ofac,hifac
          REAL(SP), INTENT(OUT) :: prob
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), DIMENSION(:), POINTER :: px,py
          end subroutine period
        END INTERFACE
        INTERFACE plgndr
          function plgndr_s(l,m,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: l,m
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: plgndr_s
          end function plgndr_s
      !BL
          function plgndr_v(l,m,x)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: l,m
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(size(x)) :: plgndr_v
          end function plgndr_v
        END INTERFACE
        INTERFACE
          function poidev(xm)
          USE nrtype
          REAL(SP), INTENT(IN) :: xm
          REAL(SP) :: poidev
          end function poidev
        END INTERFACE
        INTERFACE
          function polcoe(x,y)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), DIMENSION(size(x)) :: polcoe
          end function polcoe
        END INTERFACE
        INTERFACE
          function polcof(xa,ya)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
          REAL(SP), DIMENSION(size(xa)) :: polcof
          end function polcof
        END INTERFACE
        INTERFACE
          subroutine poldiv(u,v,q,r)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
          REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
          end subroutine poldiv
        END INTERFACE
        INTERFACE
          subroutine polin2(x1a,x2a,ya,x1,x2,y,dy)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
          REAL(SP), INTENT(IN) :: x1,x2
          REAL(SP), INTENT(OUT) :: y,dy
          end subroutine polin2
        END INTERFACE
        INTERFACE
          subroutine polint(xa,ya,x,y,dy)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
          REAL(DP), INTENT(IN) :: x
          REAL(DP), INTENT(OUT) :: y,dy
          end subroutine polint
        END INTERFACE
        INTERFACE
          subroutine powell(p,xi,ftol,iter,fret)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
          INTEGER(I4B), INTENT(OUT) :: iter
          REAL(SP), INTENT(IN) :: ftol
          REAL(SP), INTENT(OUT) :: fret
          end subroutine powell
        END INTERFACE
        INTERFACE
          function predic(data,d,nfut)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data,d
          INTEGER(I4B), INTENT(IN) :: nfut
          REAL(SP), DIMENSION(nfut) :: predic
          end function predic
        END INTERFACE
        INTERFACE
          function probks(alam)
          USE nrtype
          REAL(SP), INTENT(IN) :: alam
          REAL(SP) :: probks
          end function probks
        END INTERFACE
        INTERFACE psdes
          subroutine psdes_s(lword,rword)
          USE nrtype
          INTEGER(I4B), INTENT(INOUT) :: lword,rword
          end subroutine psdes_s
      !BL
          subroutine psdes_v(lword,rword)
          USE nrtype
          INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
          end subroutine psdes_v
        END INTERFACE
        INTERFACE
          subroutine pwt(a,isign)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine pwt
        END INTERFACE
        INTERFACE
          subroutine pwtset(n)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          end subroutine pwtset
        END INTERFACE
        INTERFACE pythag
          function pythag_dp(a,b)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP) :: pythag_dp
          end function pythag_dp
      !BL
          function pythag_sp(a,b)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP) :: pythag_sp
          end function pythag_sp
        END INTERFACE
        INTERFACE
          subroutine pzextr(iest,xest,yest,yz,dy)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: iest
          REAL(SP), INTENT(IN) :: xest
          REAL(SP), DIMENSION(:), INTENT(IN) :: yest
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
          end subroutine pzextr
        END INTERFACE
        INTERFACE
          subroutine qrdcmp(a,c,d,sing)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
          LOGICAL(LGT), INTENT(OUT) :: sing
          end subroutine qrdcmp
        END INTERFACE
        INTERFACE
          function qromb(func,a,b)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP) :: qromb
          INTERFACE
            function func(x)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end function qromb
        END INTERFACE
        INTERFACE
          function qromo(func,a,b,choose)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP) :: qromo
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          INTERFACE
            subroutine choose(funk,aa,bb,s,n)
            USE nrtype
            REAL(SP), INTENT(IN) :: aa,bb
            REAL(SP), INTENT(INOUT) :: s
            INTEGER(I4B), INTENT(IN) :: n
            INTERFACE
              function funk(x)
              USE nrtype
              REAL(SP), DIMENSION(:), INTENT(IN) :: x
              REAL(SP), DIMENSION(size(x)) :: funk
              end function funk
            END INTERFACE
            end subroutine choose
          END INTERFACE
          end function qromo
        END INTERFACE
        INTERFACE
          subroutine qroot(p,b,c,eps)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: p
          REAL(SP), INTENT(INOUT) :: b,c
          REAL(SP), INTENT(IN) :: eps
          end subroutine qroot
        END INTERFACE
        INTERFACE
          subroutine qrsolv(a,c,d,b)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
          end subroutine qrsolv
        END INTERFACE
        INTERFACE
          subroutine qrupdt(r,qt,u,v)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
          REAL(SP), DIMENSION(:), INTENT(IN) :: v
          end subroutine qrupdt
        END INTERFACE
        INTERFACE
          function qsimp(func,a,b)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP) :: qsimp
          INTERFACE
            function func(x)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end function qsimp
        END INTERFACE
        INTERFACE
          function qtrap(func,a,b)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP) :: qtrap
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end function qtrap
        END INTERFACE
        INTERFACE
          subroutine quadct(x,y,xx,yy,fa,fb,fc,fd)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y
          REAL(SP), DIMENSION(:), INTENT(IN) :: xx,yy
          REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
          end subroutine quadct
        END INTERFACE
        INTERFACE
          subroutine quadmx(a)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
          end subroutine quadmx
        END INTERFACE
        INTERFACE
          subroutine quadvl(x,y,fa,fb,fc,fd)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y
          REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
          end subroutine quadvl
        END INTERFACE
        INTERFACE
          function ran(idum)
          INTEGER(selected_int_kind(9)), INTENT(INOUT) :: idum
          REAL :: ran
          end function ran
        END INTERFACE
        INTERFACE ran0
          subroutine ran0_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine ran0_s
      !BL
          subroutine ran0_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine ran0_v
        END INTERFACE
        INTERFACE ran1
          subroutine ran1_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine ran1_s
      !BL
          subroutine ran1_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine ran1_v
        END INTERFACE
        INTERFACE ran2
          subroutine ran2_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine ran2_s
      !BL
          subroutine ran2_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine ran2_v
        END INTERFACE
        INTERFACE ran3
          subroutine ran3_s(harvest)
          USE nrtype
          REAL(SP), INTENT(OUT) :: harvest
          end subroutine ran3_s
      !BL
          subroutine ran3_v(harvest)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
          end subroutine ran3_v
        END INTERFACE
        INTERFACE
          subroutine ratint(xa,ya,x,y,dy)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
          REAL(SP), INTENT(IN) :: x
          REAL(SP), INTENT(OUT) :: y,dy
          end subroutine ratint
        END INTERFACE
        INTERFACE
          subroutine ratlsq(func,a,b,mm,kk,cof,dev)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          INTEGER(I4B), INTENT(IN) :: mm,kk
          REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
          REAL(DP), INTENT(OUT) :: dev
          INTERFACE
            function func(x)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end subroutine ratlsq
        END INTERFACE
        INTERFACE ratval
          function ratval_s(x,cof,mm,kk)
          USE nrtype
          REAL(DP), INTENT(IN) :: x
          INTEGER(I4B), INTENT(IN) :: mm,kk
          REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
          REAL(DP) :: ratval_s
          end function ratval_s
      !BL
          function ratval_v(x,cof,mm,kk)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          INTEGER(I4B), INTENT(IN) :: mm,kk
          REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
          REAL(DP), DIMENSION(size(x)) :: ratval_v
          end function ratval_v
        END INTERFACE
        INTERFACE rc
          function rc_s(x,y)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y
          REAL(SP) :: rc_s
          end function rc_s
      !BL
          function rc_v(x,y)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), DIMENSION(size(x)) :: rc_v
          end function rc_v
        END INTERFACE
        INTERFACE rd
          function rd_s(x,y,z)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y,z
          REAL(SP) :: rd_s
          end function rd_s
      !BL
          function rd_v(x,y,z)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
          REAL(SP), DIMENSION(size(x)) :: rd_v
          end function rd_v
        END INTERFACE
        INTERFACE realft
          subroutine realft_dp(data,isign,zdata)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
          end subroutine realft_dp
      !BL
          subroutine realft_sp(data,isign,zdata)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
          INTEGER(I4B), INTENT(IN) :: isign
          COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
          end subroutine realft_sp
        END INTERFACE
        INTERFACE
          RECURSIVE function recur1(a,b) RESULT(u)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
          REAL(SP), DIMENSION(size(a)) :: u
          end function recur1
        END INTERFACE
        INTERFACE
          function recur2(a,b,c)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
          REAL(SP), DIMENSION(size(a)) :: recur2
          end function recur2
        END INTERFACE
        INTERFACE
          subroutine relax(u,rhs)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
          end subroutine relax
        END INTERFACE
        INTERFACE
          subroutine relax2(u,rhs)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
          end subroutine relax2
        END INTERFACE
        INTERFACE
        function resid(u,rhs)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
          REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
          end function resid
        END INTERFACE
        INTERFACE rf
          function rf_s(x,y,z)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y,z
          REAL(SP) :: rf_s
          end function rf_s
      !BL
          function rf_v(x,y,z)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
          REAL(SP), DIMENSION(size(x)) :: rf_v
          end function rf_v
        END INTERFACE
        INTERFACE rj
          function rj_s(x,y,z,p)
          USE nrtype
          REAL(SP), INTENT(IN) :: x,y,z,p
          REAL(SP) :: rj_s
          end function rj_s
      !BL
          function rj_v(x,y,z,p)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
          REAL(SP), DIMENSION(size(x)) :: rj_v
          end function rj_v
        END INTERFACE
        INTERFACE
          subroutine rk4(y,dydx,x,h,yout,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
          REAL(SP), INTENT(IN) :: x,h
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine rk4
        END INTERFACE
        INTERFACE
          subroutine rkck(y,dydx,x,h,yout,yerr,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
          REAL(SP), INTENT(IN) :: x,h
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine rkck
        END INTERFACE
        INTERFACE
          subroutine rkdumb(vstart,x1,x2,nstep,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: vstart
          REAL(SP), INTENT(IN) :: x1,x2
          INTEGER(I4B), INTENT(IN) :: nstep
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine rkdumb
        END INTERFACE
        INTERFACE
          subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
          REAL(SP), INTENT(INOUT) :: x
          REAL(SP), INTENT(IN) :: htry,eps
          REAL(SP), INTENT(OUT) :: hdid,hnext
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine rkqs
        END INTERFACE
        INTERFACE
          subroutine rlft2(data,spec,speq,isign)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
          COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: spec
          COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: speq
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine rlft2
        END INTERFACE
        INTERFACE
          subroutine rlft3(data,spec,speq,isign)
          USE nrtype
          REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
          COMPLEX(SPC), DIMENSION(:,:,:), INTENT(OUT) :: spec
          COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: speq
          INTEGER(I4B), INTENT(IN) :: isign
          end subroutine rlft3
        END INTERFACE
        INTERFACE
          subroutine rotate(r,qt,i,a,b)
          USE nrtype
          REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
          INTEGER(I4B), INTENT(IN) :: i
          REAL(SP), INTENT(IN) :: a,b
          end subroutine rotate
        END INTERFACE
        INTERFACE
          subroutine rsolv(a,d,b)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          REAL(SP), DIMENSION(:), INTENT(IN) :: d
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
          end subroutine rsolv
        END INTERFACE
        INTERFACE
          function rstrct(uf)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
          REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
          end function rstrct
        END INTERFACE
        INTERFACE
          function rtbis(func,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: rtbis
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function rtbis
        END INTERFACE
        INTERFACE
          function rtflsp(func,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: rtflsp
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function rtflsp
        END INTERFACE
        INTERFACE
          function rtnewt(funcd,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: rtnewt
          INTERFACE
            subroutine funcd(x,fval,fderiv)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), INTENT(OUT) :: fval,fderiv
            end subroutine funcd
          END INTERFACE
          end function rtnewt
        END INTERFACE
        INTERFACE
          function rtsafe(funcd,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: rtsafe
          INTERFACE
            subroutine funcd(x,fval,fderiv)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), INTENT(OUT) :: fval,fderiv
            end subroutine funcd
          END INTERFACE
          end function rtsafe
        END INTERFACE
        INTERFACE
          function rtsec(func,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: rtsec
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function rtsec
        END INTERFACE
        INTERFACE
          subroutine rzextr(iest,xest,yest,yz,dy)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: iest
          REAL(SP), INTENT(IN) :: xest
          REAL(SP), DIMENSION(:), INTENT(IN) :: yest
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
          end subroutine rzextr
        END INTERFACE
        INTERFACE
          function savgol(nl,nrr,ld,m)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
          REAL(SP), DIMENSION(nl+nrr+1) :: savgol
          end function savgol
        END INTERFACE
        INTERFACE
          subroutine scrsho(func)
          USE nrtype
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine scrsho
        END INTERFACE
        INTERFACE
          function select(k,arr)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: k
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          REAL(SP) :: select
          end function select
        END INTERFACE
        INTERFACE
          function select_bypack(k,arr)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: k
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          REAL(SP) :: select_bypack
          end function select_bypack
        END INTERFACE
        INTERFACE
          subroutine select_heap(arr,heap)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: arr
          REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
          end subroutine select_heap
        END INTERFACE
        INTERFACE
          function select_inplace(k,arr)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: k
          REAL(SP), DIMENSION(:), INTENT(IN) :: arr
          REAL(SP) :: select_inplace
          end function select_inplace
        END INTERFACE
        INTERFACE
          subroutine simplx(a,m1,m2,m3,icase,izrov,iposv)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          INTEGER(I4B), INTENT(IN) :: m1,m2,m3
          INTEGER(I4B), INTENT(OUT) :: icase
          INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: izrov,iposv
          end subroutine simplx
        END INTERFACE
        INTERFACE
          subroutine simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
          USE nrtype
          REAL(SP), INTENT(IN) :: xs,htot
          REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
          INTEGER(I4B), INTENT(IN) :: nstep
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine simpr
        END INTERFACE
        INTERFACE
          subroutine sinft(y)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          end subroutine sinft
        END INTERFACE
        INTERFACE
          subroutine slvsm2(u,rhs)
          USE nrtype
          REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
          REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
          end subroutine slvsm2
        END INTERFACE
        INTERFACE
          subroutine slvsml(u,rhs)
          USE nrtype
          REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
          REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
          end subroutine slvsml
        END INTERFACE
        INTERFACE
          subroutine sncndn(uu,emmc,sn,cn,dn)
          USE nrtype
          REAL(SP), INTENT(IN) :: uu,emmc
          REAL(SP), INTENT(OUT) :: sn,cn,dn
          end subroutine sncndn
        END INTERFACE
        INTERFACE
          function snrm(sx,itol)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: sx
          INTEGER(I4B), INTENT(IN) :: itol
          REAL(DP) :: snrm
          end function snrm
        END INTERFACE
        INTERFACE
          subroutine sobseq(x,init)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x
          INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
          end subroutine sobseq
        END INTERFACE
        INTERFACE
          subroutine solvde(itmax,conv,slowc,scalv,indexv,nb,y)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: itmax,nb
          REAL(SP), INTENT(IN) :: conv,slowc
          REAL(SP), DIMENSION(:), INTENT(IN) :: scalv
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: y
          end subroutine solvde
        END INTERFACE
        INTERFACE
          subroutine sor(a,b,c,d,e,f,u,rjac)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
          REAL(DP), INTENT(IN) :: rjac
          end subroutine sor
        END INTERFACE
        INTERFACE
          subroutine sort(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort
        END INTERFACE
        INTERFACE
          subroutine sort2(arr,slave)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
          end subroutine sort2
        END INTERFACE
        INTERFACE
          subroutine sort3(arr,slave1,slave2)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
          end subroutine sort3
        END INTERFACE
        INTERFACE
          subroutine sort_bypack(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_bypack
        END INTERFACE
        INTERFACE
          subroutine sort_byreshape(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_byreshape
        END INTERFACE
        INTERFACE
          subroutine sort_heap(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_heap
        END INTERFACE
        INTERFACE
          subroutine sort_pick(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_pick
        END INTERFACE
        INTERFACE
          subroutine sort_radix(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_radix
        END INTERFACE
        INTERFACE
          subroutine sort_shell(arr)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          end subroutine sort_shell
        END INTERFACE
        INTERFACE
          subroutine spctrm(p,k,ovrlap,unit,n_window)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(OUT) :: p
          INTEGER(I4B), INTENT(IN) :: k
          LOGICAL(LGT), INTENT(IN) :: ovrlap
          INTEGER(I4B), OPTIONAL, INTENT(IN) :: n_window,unit
          end subroutine spctrm
        END INTERFACE
        INTERFACE
          subroutine spear(data1,data2,d,zd,probd,rs,probrs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
          end subroutine spear
        END INTERFACE
        INTERFACE sphbes
          subroutine sphbes_s(n,x,sj,sy,sjp,syp)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: x
          REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
          end subroutine sphbes_s
      !BL
          subroutine sphbes_v(n,x,sj,sy,sjp,syp)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
          end subroutine sphbes_v
        END INTERFACE
        INTERFACE
          subroutine splie2(x1a,x2a,ya,y2a)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
          end subroutine splie2
        END INTERFACE
        INTERFACE
          function splin2(x1a,x2a,ya,y2a,x1,x2)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
          REAL(SP), INTENT(IN) :: x1,x2
          REAL(SP) :: splin2
          end function splin2
        END INTERFACE
        INTERFACE
          subroutine spline(x,y,yp1,ypn,y2)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
          REAL(SP), INTENT(IN) :: yp1,ypn
          REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
          end subroutine spline
        END INTERFACE
        INTERFACE
          function splint(xa,ya,y2a,x)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
          REAL(SP), INTENT(IN) :: x
          REAL(SP) :: splint
          end function splint
        END INTERFACE
        INTERFACE sprsax
          subroutine sprsax_dp(sa,x,b)
          USE nrtype
          TYPE(sprs2_dp), INTENT(IN) :: sa
          REAL(DP), DIMENSION (:), INTENT(IN) :: x
          REAL(DP), DIMENSION (:), INTENT(OUT) :: b
          end subroutine sprsax_dp
      !BL
          subroutine sprsax_sp(sa,x,b)
          USE nrtype
          TYPE(sprs2_sp), INTENT(IN) :: sa
          REAL(SP), DIMENSION (:), INTENT(IN) :: x
          REAL(SP), DIMENSION (:), INTENT(OUT) :: b
          end subroutine sprsax_sp
        END INTERFACE
        INTERFACE sprsdiag
          subroutine sprsdiag_dp(sa,b)
          USE nrtype
          TYPE(sprs2_dp), INTENT(IN) :: sa
          REAL(DP), DIMENSION(:), INTENT(OUT) :: b
          end subroutine sprsdiag_dp
      !BL
          subroutine sprsdiag_sp(sa,b)
          USE nrtype
          TYPE(sprs2_sp), INTENT(IN) :: sa
          REAL(SP), DIMENSION(:), INTENT(OUT) :: b
          end subroutine sprsdiag_sp
        END INTERFACE
        INTERFACE sprsin
          subroutine sprsin_sp(a,thresh,sa)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
          REAL(SP), INTENT(IN) :: thresh
          TYPE(sprs2_sp), INTENT(OUT) :: sa
          end subroutine sprsin_sp
      !BL
          subroutine sprsin_dp(a,thresh,sa)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
          REAL(DP), INTENT(IN) :: thresh
          TYPE(sprs2_dp), INTENT(OUT) :: sa
          end subroutine sprsin_dp
        END INTERFACE
        INTERFACE
          subroutine sprstp(sa)
          USE nrtype
          TYPE(sprs2_sp), INTENT(INOUT) :: sa
          end subroutine sprstp
        END INTERFACE
        INTERFACE sprstx
          subroutine sprstx_dp(sa,x,b)
          USE nrtype
          TYPE(sprs2_dp), INTENT(IN) :: sa
          REAL(DP), DIMENSION (:), INTENT(IN) :: x
          REAL(DP), DIMENSION (:), INTENT(OUT) :: b
          end subroutine sprstx_dp
      !BL
          subroutine sprstx_sp(sa,x,b)
          USE nrtype
          TYPE(sprs2_sp), INTENT(IN) :: sa
          REAL(SP), DIMENSION (:), INTENT(IN) :: x
          REAL(SP), DIMENSION (:), INTENT(OUT) :: b
          end subroutine sprstx_sp
        END INTERFACE
        INTERFACE
          subroutine stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
          REAL(SP), INTENT(IN) :: htry,eps
          REAL(SP), INTENT(INOUT) :: x
          REAL(SP), INTENT(OUT) :: hdid,hnext
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine stifbs
        END INTERFACE
        INTERFACE
          subroutine stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
          REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
          REAL(SP), INTENT(INOUT) :: x
          REAL(SP), INTENT(IN) :: htry,eps
          REAL(SP), INTENT(OUT) :: hdid,hnext
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine stiff
        END INTERFACE
        INTERFACE
          subroutine stoerm(y,d2y,xs,htot,nstep,yout,derivs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
          REAL(SP), INTENT(IN) :: xs,htot
          INTEGER(I4B), INTENT(IN) :: nstep
          REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
          INTERFACE
            subroutine derivs(x,y,dydx)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
            end subroutine derivs
          END INTERFACE
          end subroutine stoerm
        END INTERFACE
        INTERFACE svbksb
          subroutine svbksb_dp(u,w,v,b,x)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
          REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
          REAL(DP), DIMENSION(:), INTENT(OUT) :: x
          end subroutine svbksb_dp
      !BL
          subroutine svbksb_sp(u,w,v,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
          REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x
          end subroutine svbksb_sp
        END INTERFACE
        INTERFACE svdcmp
          subroutine svdcmp_dp(a,w,v)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(DP), DIMENSION(:), INTENT(OUT) :: w
          REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
          end subroutine svdcmp_dp
      !BL
          subroutine svdcmp_sp(a,w,v)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: w
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
          end subroutine svdcmp_sp
        END INTERFACE
        INTERFACE
          subroutine svdfit(x,y,sig,a,v,w,chisq,funcs)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
          REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
          REAL(SP), INTENT(OUT) :: chisq
          INTERFACE
            function funcs(x,n)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            INTEGER(I4B), INTENT(IN) :: n
            REAL(SP), DIMENSION(n) :: funcs
            end function funcs
          END INTERFACE
          end subroutine svdfit
        END INTERFACE
        INTERFACE
          subroutine svdvar(v,w,cvm)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
          REAL(SP), DIMENSION(:), INTENT(IN) :: w
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
          end subroutine svdvar
        END INTERFACE
        INTERFACE
          function toeplz(r,y)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
          REAL(SP), DIMENSION(size(y)) :: toeplz
          end function toeplz
        END INTERFACE
        INTERFACE
          subroutine tptest(data1,data2,t,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(SP), INTENT(OUT) :: t,prob
          end subroutine tptest
        END INTERFACE
        INTERFACE
          subroutine tqli(d,e,z)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
          REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
          end subroutine tqli
        END INTERFACE
        INTERFACE
          subroutine trapzd(func,a,b,s,n)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP), INTENT(INOUT) :: s
          INTEGER(I4B), INTENT(IN) :: n
          INTERFACE
            function func(x)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
            end function func
          END INTERFACE
          end subroutine trapzd
        END INTERFACE
        INTERFACE
          subroutine tred2(a,d,e,novectors)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
          LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
          end subroutine tred2
        END INTERFACE
      !	On a purely serial machine, for greater efficiency, remove
      !	the generic name tridag from the following interface,
      !	and put it on the next one after that.
        INTERFACE tridag
          RECURSIVE subroutine tridag_par(a,b,c,r,u)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
          REAL(SP), DIMENSION(:), INTENT(OUT) :: u
          end subroutine tridag_par
        END INTERFACE
        INTERFACE
          subroutine tridag_ser(a,b,c,r,u)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
          REAL(SP), DIMENSION(:), INTENT(OUT) :: u
          end subroutine tridag_ser
        END INTERFACE
        INTERFACE
          subroutine ttest(data1,data2,t,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(SP), INTENT(OUT) :: t,prob
          end subroutine ttest
        END INTERFACE
        INTERFACE
          subroutine tutest(data1,data2,t,prob)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          REAL(SP), INTENT(OUT) :: t,prob
          end subroutine tutest
        END INTERFACE
        INTERFACE
          subroutine twofft(data1,data2,fft1,fft2)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
          COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
          end subroutine twofft
        END INTERFACE
        INTERFACE
          function vander(x,q)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
          REAL(DP), DIMENSION(size(x)) :: vander
          end function vander
        END INTERFACE
        INTERFACE
          subroutine vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: region
          INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
          REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
          INTERFACE
            function func(pt,wgt)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(IN) :: pt
            REAL(SP), INTENT(IN) :: wgt
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine vegas
        END INTERFACE
        INTERFACE
          subroutine voltra(t0,h,t,f,g,ak)
          USE nrtype
          REAL(SP), INTENT(IN) :: t0,h
          REAL(SP), DIMENSION(:), INTENT(OUT) :: t
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: f
          INTERFACE
            function g(t)
            USE nrtype
            REAL(SP), INTENT(IN) :: t
            REAL(SP), DIMENSION(:), POINTER :: g
            end function g
      !BL
            function ak(t,s)
            USE nrtype
            REAL(SP), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(:,:), POINTER :: ak
            end function ak
          END INTERFACE
          end subroutine voltra
        END INTERFACE
        INTERFACE
          subroutine wt1(a,isign,wtstep)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          INTEGER(I4B), INTENT(IN) :: isign
          INTERFACE
            subroutine wtstep(a,isign)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
            INTEGER(I4B), INTENT(IN) :: isign
            end subroutine wtstep
          END INTERFACE
          end subroutine wt1
        END INTERFACE
        INTERFACE
          subroutine wtn(a,nn,isign,wtstep)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
          INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
          INTEGER(I4B), INTENT(IN) :: isign
          INTERFACE
            subroutine wtstep(a,isign)
            USE nrtype
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
            INTEGER(I4B), INTENT(IN) :: isign
            end subroutine wtstep
          END INTERFACE
          end subroutine wtn
        END INTERFACE
        INTERFACE
          function wwghts(n,h,kermom)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          REAL(SP), INTENT(IN) :: h
          REAL(SP), DIMENSION(n) :: wwghts
          INTERFACE
            function kermom(y,m)
            USE nrtype
            REAL(DP), INTENT(IN) :: y
            INTEGER(I4B), INTENT(IN) :: m
            REAL(DP), DIMENSION(m) :: kermom
            end function kermom
          END INTERFACE
          end function wwghts
        END INTERFACE
        INTERFACE
          subroutine zbrac(func,x1,x2,succes)
          USE nrtype
          REAL(SP), INTENT(INOUT) :: x1,x2
          LOGICAL(LGT), INTENT(OUT) :: succes
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine zbrac
        END INTERFACE
        INTERFACE
          subroutine zbrak(func,x1,x2,n,xb1,xb2,nb)
          USE nrtype
          INTEGER(I4B), INTENT(IN) :: n
          INTEGER(I4B), INTENT(OUT) :: nb
          REAL(SP), INTENT(IN) :: x1,x2
          REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end subroutine zbrak
        END INTERFACE
        INTERFACE
          function zbrent(func,x1,x2,tol)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,tol
          REAL(SP) :: zbrent
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function zbrent
        END INTERFACE
        INTERFACE
          subroutine zrhqr(a,rtr,rti)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
          end subroutine zrhqr
        END INTERFACE
        INTERFACE
          function zriddr(func,x1,x2,xacc)
          USE nrtype
          REAL(SP), INTENT(IN) :: x1,x2,xacc
          REAL(SP) :: zriddr
          INTERFACE
            function func(x)
            USE nrtype
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
            end function func
          END INTERFACE
          end function zriddr
        END INTERFACE
        INTERFACE
          subroutine zroots(a,roots,polish)
          USE nrtype
          COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
          COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
          LOGICAL(LGT), INTENT(IN) :: polish
          end subroutine zroots
        END INTERFACE
      END MODULE nr


!-----------------------------------------------------------------------
! nrutil.f90
      MODULE nrutil
        USE nrtype
        implicit none
        INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
        INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
        INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
        INTEGER(I4B), PARAMETER :: NPAR_POLY=8
        INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
        INTERFACE array_copy
          MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
        END INTERFACE
        INTERFACE swap
          MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
            swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
            masked_swap_rs,masked_swap_rv,masked_swap_rm
        END INTERFACE
        INTERFACE reallocate
          MODULE PROCEDURE reallocate_rv,reallocate_rm, &
            reallocate_iv,reallocate_im,reallocate_hv
        END INTERFACE
        INTERFACE imaxloc
          MODULE PROCEDURE imaxloc_r,imaxloc_i
        END INTERFACE
        INTERFACE assert
          MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
        END INTERFACE
        INTERFACE assert_eq
          MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
        END INTERFACE
        INTERFACE arth
          MODULE PROCEDURE arth_r, arth_d, arth_i
        END INTERFACE
        INTERFACE geop
          MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
        END INTERFACE
        INTERFACE cumsum
          MODULE PROCEDURE cumsum_r,cumsum_i
        END INTERFACE
        INTERFACE poly
          MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv, &
            poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
        END INTERFACE
        INTERFACE poly_term
          MODULE PROCEDURE poly_term_rr,poly_term_cc
        END INTERFACE
        INTERFACE outerprod
          MODULE PROCEDURE outerprod_r,outerprod_d
        END INTERFACE
        INTERFACE outerdiff
          MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
        END INTERFACE
        INTERFACE scatter_add
          MODULE PROCEDURE scatter_add_r,scatter_add_d
        END INTERFACE
        INTERFACE scatter_max
          MODULE PROCEDURE scatter_max_r,scatter_max_d
        END INTERFACE
        INTERFACE diagadd
          MODULE PROCEDURE diagadd_rv,diagadd_r
        END INTERFACE
        INTERFACE diagmult
          MODULE PROCEDURE diagmult_rv,diagmult_r
        END INTERFACE
        INTERFACE get_diag
          MODULE PROCEDURE get_diag_rv, get_diag_dv
        END INTERFACE
        INTERFACE put_diag
          MODULE PROCEDURE put_diag_rv, put_diag_r
        END INTERFACE
      CONTAINS
      !BL
        subroutine array_copy_r(src,dest,n_copied,n_not_copied)
        REAL(SP), DIMENSION(:), INTENT(IN) :: src
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied = min(size(src),size(dest))
        n_not_copied=size(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        end subroutine array_copy_r
      !BL
        subroutine array_copy_d(src,dest,n_copied,n_not_copied)
        REAL(DP), DIMENSION(:), INTENT(IN) :: src
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied = min(size(src),size(dest))
        n_not_copied=size(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        end subroutine array_copy_d
      !BL
        subroutine array_copy_i(src,dest,n_copied,n_not_copied)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied = min(size(src),size(dest))
        n_not_copied=size(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        end subroutine array_copy_i
      !BL
      !BL
        subroutine swap_i(a,b)
        INTEGER(I4B), INTENT(INOUT) :: a,b
        INTEGER(I4B) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_i
      !BL
        subroutine swap_r(a,b)
        REAL(SP), INTENT(INOUT) :: a,b
        REAL(SP) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_r
      !BL
        subroutine swap_rv(a,b)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(SP), DIMENSION(SIZE(a)) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_rv
      !BL
        subroutine swap_c(a,b)
        COMPLEX(SPC), INTENT(INOUT) :: a,b
        COMPLEX(SPC) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_c
      !BL
        subroutine swap_cv(a,b)
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_cv
      !BL
        subroutine swap_cm(a,b)
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_cm
      !BL
        subroutine swap_z(a,b)
        COMPLEX(DPC), INTENT(INOUT) :: a,b
        COMPLEX(DPC) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_z
      !BL
        subroutine swap_zv(a,b)
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_zv
      !BL
        subroutine swap_zm(a,b)
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum = a
        a = b
        b = dum
        end subroutine swap_zm
      !BL
        subroutine masked_swap_rs(a,b,mask)
        REAL(SP), INTENT(INOUT) :: a,b
        LOGICAL(LGT), INTENT(IN) :: mask
        REAL(SP) :: swp
        if (mask) then
          swp = a
          a = b
          b = swp
        endif
        end subroutine masked_swap_rs
      !BL
        subroutine masked_swap_rv(a,b,mask)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a)) :: swp
        where (mask)
          swp = a
          a = b
          b = swp
        end where
        end subroutine masked_swap_rv
      !BL
        subroutine masked_swap_rm(a,b,mask)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
          swp = a
          a = b
          b = swp
        end where
        end subroutine masked_swap_rm
      !BL
      !BL
        function reallocate_rv(p,n)
        REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_rv(n),stat=ierr)
        if (ierr /= 0) call &
          nrerror('reallocate_rv: problem in attempt to allocate memory')
        if (.not. associated(p)) return
        nold=size(p)
        reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        end function reallocate_rv
      !BL
        function reallocate_iv(p,n)
        INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_iv(n),stat=ierr)
        if (ierr /= 0) call &
          nrerror('reallocate_iv: problem in attempt to allocate memory')
        if (.not. associated(p)) return
        nold=size(p)
        reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        end function reallocate_iv
      !BL
        function reallocate_hv(p,n)
        CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_hv(n),stat=ierr)
        if (ierr /= 0) call &
          nrerror('reallocate_hv: problem in attempt to allocate memory')
        if (.not. associated(p)) return
        nold=size(p)
        reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        end function reallocate_hv
      !BL
        function reallocate_rm(p,n,m)
        REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_rm(n,m),stat=ierr)
        if (ierr /= 0) call &
          nrerror('reallocate_rm: problem in attempt to allocate memory')
        if (.not. associated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm(1:min(nold,n),1:min(mold,m))=&
          p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
        end function reallocate_rm
      !BL
        function reallocate_im(p,n,m)
        INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_im(n,m),stat=ierr)
        if (ierr /= 0) call &
          nrerror('reallocate_im: problem in attempt to allocate memory')
        if (.not. associated(p)) return
        nold=size(p,1)
        mold=size(p,2)
        reallocate_im(1:min(nold,n),1:min(mold,m))=&
          p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
        end function reallocate_im
      !BL
        function ifirstloc(mask)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        INTEGER(I4B) :: ifirstloc
        INTEGER(I4B), DIMENSION(1) :: loc
        loc = maxloc(merge(1,0,mask))
        ifirstloc=loc(1)
        if (.not. mask(ifirstloc)) ifirstloc = size(mask)+1
        end function ifirstloc
      !BL
        function imaxloc_r(arr)
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B) :: imaxloc_r
        INTEGER(I4B), DIMENSION(1) :: imax
        imax = maxloc(arr(:))
        imaxloc_r=imax(1)
        end function imaxloc_r
      !BL
        function imaxloc_i(iarr)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
        INTEGER(I4B), DIMENSION(1) :: imax
        INTEGER(I4B) :: imaxloc_i
        imax = maxloc(iarr(:))
        imaxloc_i=imax(1)
        end function imaxloc_i
      !BL
        function iminloc(arr)
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(1) :: imin
        INTEGER(I4B) :: iminloc
        imin=minloc(arr(:))
        iminloc=imin(1)
        end function iminloc
      !BL
        subroutine assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
          write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
          STOP 'program terminated by assert1'
        endif
        end subroutine assert1
      !BL
        subroutine assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
          write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
          STOP 'program terminated by assert2'
        endif
        end subroutine assert2
      !BL
        subroutine assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
          write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
          STOP 'program terminated by assert3'
        endif
        end subroutine assert3
      !BL
        subroutine assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
          write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
          STOP 'program terminated by assert4'
        endif
        end subroutine assert4
      !BL
        subroutine assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
          write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
          STOP 'program terminated by assert_v'
        endif
        end subroutine assert_v
      !BL
        function assert_eq2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2
        INTEGER :: assert_eq2
        if (n1 == n2) then
          assert_eq2 = n1
        else
          write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
          STOP 'program terminated by assert_eq2'
        endif
        end function assert_eq2
      !BL
        function assert_eq3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3
        INTEGER :: assert_eq3
        if (n1 == n2 .and. n2 == n3) then
          assert_eq3 = n1
        else
          write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
          STOP 'program terminated by assert_eq3'
        endif
        end function assert_eq3
      !BL
        function assert_eq4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3,n4
        INTEGER :: assert_eq4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
          assert_eq4 = n1
        else
          write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
          STOP 'program terminated by assert_eq4'
        endif
        end function assert_eq4
      !BL
        function assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER :: assert_eqn
        if (all(nn(2:) == nn(1))) then
          assert_eqn=nn(1)
        else
          write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
          STOP 'program terminated by assert_eqn'
        endif
        end function assert_eqn
      !BL
        subroutine nrerror(string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        write (*,*) 'nrerror: ',string
        STOP 'program terminated by nrerror'
        end subroutine nrerror
      !BL
        function arth_r(first,increment,n)
        REAL(SP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: arth_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        if (n > 0) arth_r(1)=first
        if (n <= NPAR_ARTH) then
          do k = 2,n
            arth_r(k)=arth_r(k-1)+increment
          enddo
        else
          do k = 2,NPAR2_ARTH
            arth_r(k)=arth_r(k-1)+increment
          enddo
          temp = increment*NPAR2_ARTH
          k = NPAR2_ARTH
          do
            if (k >= n) exit
            k2 = k+k
            arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
            temp = temp+temp
            k = k2
          enddo
        endif
        end function arth_r
      !BL
        function arth_d(first,increment,n)
        REAL(DP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: arth_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        if (n > 0) arth_d(1)=first
        if (n <= NPAR_ARTH) then
          do k = 2,n
            arth_d(k)=arth_d(k-1)+increment
          enddo
        else
          do k = 2,NPAR2_ARTH
            arth_d(k)=arth_d(k-1)+increment
          enddo
          temp = increment*NPAR2_ARTH
          k = NPAR2_ARTH
          do
            if (k >= n) exit
            k2 = k+k
            arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
            temp = temp+temp
            k = k2
          enddo
        endif
        end function arth_d
      !BL
        function arth_i(first,increment,n)
        INTEGER(I4B), INTENT(IN) :: first,increment,n
        INTEGER(I4B), DIMENSION(n) :: arth_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
          do k = 2,n
            arth_i(k)=arth_i(k-1)+increment
          enddo
        else
          do k = 2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
          enddo
          temp = increment*NPAR2_ARTH
          k = NPAR2_ARTH
          do
            if (k >= n) exit
            k2 = k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp = temp+temp
            k = k2
          enddo
        endif
        end function arth_i
      !BL
      !BL
        function geop_r(first,factor,n)
        REAL(SP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: geop_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        if (n > 0) geop_r(1)=first
        if (n <= NPAR_GEOP) then
          do k = 2,n
            geop_r(k)=geop_r(k-1)*factor
          enddo
        else
          do k = 2,NPAR2_GEOP
            geop_r(k)=geop_r(k-1)*factor
          enddo
          temp = factor**NPAR2_GEOP
          k = NPAR2_GEOP
          do
            if (k >= n) exit
            k2 = k+k
            geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
            temp = temp*temp
            k = k2
          enddo
        endif
        end function geop_r
      !BL
        function geop_d(first,factor,n)
        REAL(DP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: geop_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        if (n > 0) geop_d(1)=first
        if (n <= NPAR_GEOP) then
          do k = 2,n
            geop_d(k)=geop_d(k-1)*factor
          enddo
        else
          do k = 2,NPAR2_GEOP
            geop_d(k)=geop_d(k-1)*factor
          enddo
          temp = factor**NPAR2_GEOP
          k = NPAR2_GEOP
          do
            if (k >= n) exit
            k2 = k+k
            geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
            temp = temp*temp
            k = k2
          enddo
        endif
        end function geop_d
      !BL
        function geop_i(first,factor,n)
        INTEGER(I4B), INTENT(IN) :: first,factor,n
        INTEGER(I4B), DIMENSION(n) :: geop_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) geop_i(1)=first
        if (n <= NPAR_GEOP) then
          do k = 2,n
            geop_i(k)=geop_i(k-1)*factor
          enddo
        else
          do k = 2,NPAR2_GEOP
            geop_i(k)=geop_i(k-1)*factor
          enddo
          temp = factor**NPAR2_GEOP
          k = NPAR2_GEOP
          do
            if (k >= n) exit
            k2 = k+k
            geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
            temp = temp*temp
            k = k2
          enddo
        endif
        end function geop_i
      !BL
        function geop_c(first,factor,n)
        COMPLEX(SP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        COMPLEX(SP), DIMENSION(n) :: geop_c
        INTEGER(I4B) :: k,k2
        COMPLEX(SP) :: temp
        if (n > 0) geop_c(1)=first
        if (n <= NPAR_GEOP) then
          do k = 2,n
            geop_c(k)=geop_c(k-1)*factor
          enddo
        else
          do k = 2,NPAR2_GEOP
            geop_c(k)=geop_c(k-1)*factor
          enddo
          temp = factor**NPAR2_GEOP
          k = NPAR2_GEOP
          do
            if (k >= n) exit
            k2 = k+k
            geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
            temp = temp*temp
            k = k2
          enddo
        endif
        end function geop_c
      !BL
        function geop_dv(first,factor,n)
        REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(size(first),n) :: geop_dv
        INTEGER(I4B) :: k,k2
        REAL(DP), DIMENSION(size(first)) :: temp
        if (n > 0) geop_dv(:,1)=first(:)
        if (n <= NPAR_GEOP) then
          do k = 2,n
            geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
          enddo
        else
          do k = 2,NPAR2_GEOP
            geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
          enddo
          temp = factor**NPAR2_GEOP
          k = NPAR2_GEOP
          do
            if (k >= n) exit
            k2 = k+k
            geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
              spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
            temp = temp*temp
            k = k2
          enddo
        endif
        end function geop_dv
      !BL
      !BL
        RECURSIVE function cumsum_r(arr,seed) RESULT(ans)
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        REAL(SP), OPTIONAL, INTENT(IN) :: seed
        REAL(SP), DIMENSION(size(arr)) :: ans
        INTEGER(I4B) :: n,j
        REAL(SP) :: sd
        n=size(arr)
        if (n == 0_i4b) return
        sd = 0.0_sp
        if (present(seed)) sd = seed
        ans(1)=arr(1)+sd
        if (n < NPAR_CUMSUM) then
          do j = 2,n
            ans(j)=ans(j-1)+arr(j)
          enddo
        else
          ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
          ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
        endif
        end function cumsum_r
      !BL
        RECURSIVE function cumsum_i(arr,seed) RESULT(ans)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
        INTEGER(I4B), DIMENSION(size(arr)) :: ans
        INTEGER(I4B) :: n,j,sd
        n=size(arr)
        if (n == 0_i4b) return
        sd = 0_i4b
        if (present(seed)) sd = seed
        ans(1)=arr(1)+sd
        if (n < NPAR_CUMSUM) then
          do j = 2,n
            ans(j)=ans(j-1)+arr(j)
          enddo
        else
          ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
          ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
        endif
        end function cumsum_i
      !BL
      !BL
        RECURSIVE function cumprod(arr,seed) RESULT(ans)
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        REAL(SP), OPTIONAL, INTENT(IN) :: seed
        REAL(SP), DIMENSION(size(arr)) :: ans
        INTEGER(I4B) :: n,j
        REAL(SP) :: sd
        n=size(arr)
        if (n == 0_i4b) return
        sd = 1.0_sp
        if (present(seed)) sd = seed
        ans(1)=arr(1)*sd
        if (n < NPAR_CUMPROD) then
          do j = 2,n
            ans(j)=ans(j-1)*arr(j)
          enddo
        else
          ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
          ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
        endif
        end function cumprod
      !BL
      !BL
        function poly_rr(x,coeffs)
        REAL(SP), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(SP) :: poly_rr
        REAL(SP) :: pow
        REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
          poly_rr = 0.0_sp
        else if (n < NPAR_POLY) then
          poly_rr=coeffs(n)
          do i = n-1,1,-1
            poly_rr=x*poly_rr+coeffs(i)
          enddo
        else
          allocate(vec(n+1))
          pow = x
          vec(1:n)=coeffs
          do
            vec(n+1)=0.0_sp
            nn=ishft(n+1,-1)
            vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
            if (nn == 1) exit
            pow = pow*pow
            n = nn
          enddo
          poly_rr=vec(1)
          deallocate(vec)
        endif
        end function poly_rr
      !BL
        function poly_dd(x,coeffs)
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(DP) :: poly_dd
        REAL(DP) :: pow
        REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
          poly_dd = 0.0_dp
        else if (n < NPAR_POLY) then
          poly_dd=coeffs(n)
          do i = n-1,1,-1
            poly_dd=x*poly_dd+coeffs(i)
          enddo
        else
          allocate(vec(n+1))
          pow = x
          vec(1:n)=coeffs
          do
            vec(n+1)=0.0_dp
            nn=ishft(n+1,-1)
            vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
            if (nn == 1) exit
            pow = pow*pow
            n = nn
          enddo
          poly_dd=vec(1)
          deallocate(vec)
        endif
        end function poly_dd
      !BL
        function poly_rc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_rc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
          poly_rc = 0.0_sp
        else if (n < NPAR_POLY) then
          poly_rc=coeffs(n)
          do i = n-1,1,-1
            poly_rc=x*poly_rc+coeffs(i)
          enddo
        else
          allocate(vec(n+1))
          pow = x
          vec(1:n)=coeffs
          do
            vec(n+1)=0.0_sp
            nn=ishft(n+1,-1)
            vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
            if (nn == 1) exit
            pow = pow*pow
            n = nn
          enddo
          poly_rc=vec(1)
          deallocate(vec)
        endif
        end function poly_rc
      !BL
        function poly_cc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_cc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
          poly_cc = 0.0_sp
        else if (n < NPAR_POLY) then
          poly_cc=coeffs(n)
          do i = n-1,1,-1
            poly_cc=x*poly_cc+coeffs(i)
          enddo
        else
          allocate(vec(n+1))
          pow = x
          vec(1:n)=coeffs
          do
            vec(n+1)=0.0_sp
            nn=ishft(n+1,-1)
            vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
            if (nn == 1) exit
            pow = pow*pow
            n = nn
          enddo
          poly_cc=vec(1)
          deallocate(vec)
        endif
        end function poly_cc
      !BL
        function poly_rrv(x,coeffs)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(SP), DIMENSION(size(x)) :: poly_rrv
        INTEGER(I4B) :: i,n,m
        m=size(coeffs)
        n=size(x)
        if (m <= 0) then
          poly_rrv = 0.0_sp
        else if (m < n .or. m < NPAR_POLY) then
          poly_rrv=coeffs(m)
          do i = m-1,1,-1
            poly_rrv=x*poly_rrv+coeffs(i)
          enddo
        else
          do i = 1,n
            poly_rrv(i)=poly_rr(x(i),coeffs)
          enddo
        endif
        end function poly_rrv
      !BL
        function poly_ddv(x,coeffs)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(DP), DIMENSION(size(x)) :: poly_ddv
        INTEGER(I4B) :: i,n,m
        m=size(coeffs)
        n=size(x)
        if (m <= 0) then
          poly_ddv = 0.0_dp
        else if (m < n .or. m < NPAR_POLY) then
          poly_ddv=coeffs(m)
          do i = m-1,1,-1
            poly_ddv=x*poly_ddv+coeffs(i)
          enddo
        else
          do i = 1,n
            poly_ddv(i)=poly_dd(x(i),coeffs)
          enddo
        endif
        end function poly_ddv
      !BL
        function poly_msk_rrv(x,coeffs,mask)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
        poly_msk_rrv = unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
        end function poly_msk_rrv
      !BL
        function poly_msk_ddv(x,coeffs,mask)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
        poly_msk_ddv = unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
        end function poly_msk_ddv
      !BL
      !BL
        RECURSIVE function poly_term_rr(a,b) RESULT(u)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a
        REAL(SP), INTENT(IN) :: b
        REAL(SP), DIMENSION(size(a)) :: u
        INTEGER(I4B) :: n,j
        n=size(a)
        if (n <= 0) return
        u(1)=a(1)
        if (n < NPAR_POLYTERM) then
          do j = 2,n
            u(j)=a(j)+b*u(j-1)
          enddo
        else
          u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
          u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
        endif
        end function poly_term_rr
      !BL
        RECURSIVE function poly_term_cc(a,b) RESULT(u)
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(SPC), INTENT(IN) :: b
        COMPLEX(SPC), DIMENSION(size(a)) :: u
        INTEGER(I4B) :: n,j
        n=size(a)
        if (n <= 0) return
        u(1)=a(1)
        if (n < NPAR_POLYTERM) then
          do j = 2,n
            u(j)=a(j)+b*u(j-1)
          enddo
        else
          u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
          u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
        endif
        end function poly_term_cc
      !BL
      !BL
        function zroots_unity(n,nn)
        INTEGER(I4B), INTENT(IN) :: n,nn
        COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
        INTEGER(I4B) :: k
        REAL(SP) :: theta
        zroots_unity(1)=1.0
        theta = TWOPI/n
        k = 1
        do
          if (k >= nn) exit
          zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
          zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
          k = 2*k
        enddo
        end function zroots_unity
      !BL
        function outerprod_r(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
        outerprod_r = spread(a,dim = 2,ncopies = size(b)) * &
          spread(b,dim = 1,ncopies = size(a))
        end function outerprod_r
      !BL
        function outerprod_d(a,b)
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
        outerprod_d = spread(a,dim = 2,ncopies = size(b)) * &
          spread(b,dim = 1,ncopies = size(a))
        end function outerprod_d
      !BL
        function outerdiv(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
        outerdiv = spread(a,dim = 2,ncopies = size(b)) / &
          spread(b,dim = 1,ncopies = size(a))
        end function outerdiv
      !BL
        function outersum(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(size(a),size(b)) :: outersum
        outersum = spread(a,dim = 2,ncopies = size(b)) + &
          spread(b,dim = 1,ncopies = size(a))
        end function outersum
      !BL
        function outerdiff_r(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
        outerdiff_r = spread(a,dim = 2,ncopies = size(b)) - &
          spread(b,dim = 1,ncopies = size(a))
        end function outerdiff_r
      !BL
        function outerdiff_d(a,b)
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
        outerdiff_d = spread(a,dim = 2,ncopies = size(b)) - &
          spread(b,dim = 1,ncopies = size(a))
        end function outerdiff_d
      !BL
        function outerdiff_i(a,b)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
        INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
        outerdiff_i = spread(a,dim = 2,ncopies = size(b)) - &
          spread(b,dim = 1,ncopies = size(a))
        end function outerdiff_i
      !BL
        function outerand(a,b)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
        LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
        outerand = spread(a,dim = 2,ncopies = size(b)) .and. &
          spread(b,dim = 1,ncopies = size(a))
        end function outerand
      !BL
        subroutine scatter_add_r(dest,source,dest_index)
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(SP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n = assert_eq2(size(source),size(dest_index),'scatter_add_r')
        m=size(dest)
        do j = 1,n
          i=dest_index(j)
          if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
        enddo
        end subroutine scatter_add_r
        subroutine scatter_add_d(dest,source,dest_index)
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n = assert_eq2(size(source),size(dest_index),'scatter_add_d')
        m=size(dest)
        do j = 1,n
          i=dest_index(j)
          if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
        enddo
        end subroutine scatter_add_d
        subroutine scatter_max_r(dest,source,dest_index)
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(SP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n = assert_eq2(size(source),size(dest_index),'scatter_max_r')
        m=size(dest)
        do j = 1,n
          i=dest_index(j)
          if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
        enddo
        end subroutine scatter_max_r
        subroutine scatter_max_d(dest,source,dest_index)
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n = assert_eq2(size(source),size(dest_index),'scatter_max_d')
        m=size(dest)
        do j = 1,n
          i=dest_index(j)
          if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
        enddo
        end subroutine scatter_max_d
      !BL
        subroutine diagadd_rv(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), DIMENSION(:), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
        do j = 1,n
          mat(j,j)=mat(j,j)+diag(j)
        enddo
        end subroutine diagadd_rv
      !BL
        subroutine diagadd_r(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j = 1,n
          mat(j,j)=mat(j,j)+diag
        enddo
        end subroutine diagadd_r
      !BL
        subroutine diagmult_rv(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), DIMENSION(:), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
        do j = 1,n
          mat(j,j)=mat(j,j)*diag(j)
        enddo
        end subroutine diagmult_rv
      !BL
        subroutine diagmult_r(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j = 1,n
          mat(j,j)=mat(j,j)*diag
        enddo
        end subroutine diagmult_r
      !BL
        function get_diag_rv(mat)
        REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
        INTEGER(I4B) :: j
        j = assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
        do j = 1,size(mat,1)
          get_diag_rv(j)=mat(j,j)
        enddo
        end function get_diag_rv
      !BL
        function get_diag_dv(mat)
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
        INTEGER(I4B) :: j
        j = assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
        do j = 1,size(mat,1)
          get_diag_dv(j)=mat(j,j)
        enddo
        end function get_diag_dv
      !BL
        subroutine put_diag_rv(diagv,mat)
        REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        INTEGER(I4B) :: j,n
        n = assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
        do j = 1,n
          mat(j,j)=diagv(j)
        enddo
        end subroutine put_diag_rv
      !BL
        subroutine put_diag_r(scal,mat)
        REAL(SP), INTENT(IN) :: scal
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j = 1,n
          mat(j,j)=scal
        enddo
        end subroutine put_diag_r
      !BL
        subroutine unit_matrix(mat)
        REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
        INTEGER(I4B) :: i,n
        n = min(size(mat,1),size(mat,2))
        mat(:,:)=0.0_sp
        do i = 1,n
          mat(i,i)=1.0_sp
        enddo
        end subroutine unit_matrix
      !BL
        function upper_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
        INTEGER(I4B) :: n
        n = 0
        if (present(extra)) n = extra
        upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
        end function upper_triangle
      !BL
        function lower_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
        INTEGER(I4B) :: n
        n = 0
        if (present(extra)) n = extra
        lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
        end function lower_triangle
      !BL
        function vabs(v)
        REAL(SP), DIMENSION(:), INTENT(IN) :: v
        REAL(SP) :: vabs
        vabs = sqrt(dot_product(v,v))
        end function vabs
      !BL
      END MODULE nrutil


