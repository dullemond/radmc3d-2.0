module mathroutines_module
use constants_module

double precision, allocatable :: exp_x(:),exp_y(:)
double precision :: exp_fact
integer :: exp_nx

contains

!-------------------------------------------------------------------
!                   REMAP FUNCTION ONTO NEW GRID
!
! This is a general purpose routine for remapping onto a new grid.
! Interpolation is done in log(f) space. Care is taken that in case
! of going to lower resolution, the average of the high-res function 
! is taken as the low-res value.
! 
! ARGUMENTS:
!  nold        Nr of grid points of old function
!  xold        Array with old x coordinates 
!  fold        Array with old function values
!  nnew        Nr of grid points of new function
!  xnew        Array with new x coordinates
!  av          =0  Do not use the averaging trick in case of lowering resolution
!              =1  Use averaging trick: Make sure that the value fnew in case
!                  of lower resolution equals the average of fold over the
!                  region of x belonging to that grid point.
!  elow        Way to treat if xnew has lower lowest value than xold
!              =0  Do not allow out of bound
!              =1  Take boundary value of fold
!              =2  Logarithmic extrapolation
!              =3  Out of bound -> 0
!              =4  Smooth version of out of bound -> 0
!  eup         Way to treat if xnew has higher highest value than xold
!              =0  Do not allow out of bound
!              =1  Take boundary value of fold
!              =2  Logarithmic extrapolation
!              =3  Out of bound -> 0
!              =4  Smooth version of out of bound -> 0
!
! RETURNS:
!  fnew        Array with new function values
!  ierror      =0  means everything went fine. Else: error:
!              =-1 xnew has lower lowest value than xold (out of bound),
!                  Aborted.
!              =+1 xnew has higher highest value than xold (out of bound)
!                  Aborted.
!-------------------------------------------------------------------
subroutine remap_function(nold,xold,fold,nnew,xnew,fnew,av,elow,eup,ierror)
  implicit none
  integer :: nold,nnew,elow,eup,ierror,av
  double precision :: xold(1:nold),fold(1:nold)
  double precision :: xnew(1:nnew),fnew(1:nnew)
  integer :: iold,inew,sgnold,sgnnew,ioldstart,ioldend,inewstart,inewend
  double precision :: eps
  !
  ! Reset ierror
  !
  ierror = 0
  !
  ! Check if x is ascending or descending, and figure out the 
  ! minimum and maximum values of xold and xnew
  !
  if(xold(nold).ge.xold(1)) then
     sgnold    = 1
     ioldstart = 1
     ioldend   = nold
  else
     sgnold    = -1
     ioldstart = nold
     ioldend   = 1
  endif
  if(xnew(nnew).ge.xnew(1)) then
     sgnnew    = 1
     inewstart = 1
     inewend   = nnew
  else
     sgnnew    = -1
     inewstart = nnew
     inewend   = 1
  endif
  !
  ! Now do a loop over the new x values, from low to high
  !
  do inew=inewstart,inewend,sgnnew
     if(xnew(inew).lt.xold(ioldstart)) then
        !
        ! Lower than lowest boundary 
        !
        if(elow.eq.0) then
           ierror = -1
           return
        elseif(elow.eq.1) then
           fnew(inew) = fold(ioldstart)
        elseif(elow.eq.2) then
           if(nold.gt.1) then
              if((fold(ioldstart+sgnold).gt.0.d0).and. &
                   (fold(ioldstart).gt.0.d0)) then
                 fnew(inew) = fold(ioldstart) *                               &
                        (fold(ioldstart+sgnold)/fold(ioldstart))**            &
                        ((log(xnew(inew))-log(xold(ioldstart)))/              &
                         (log(xold(ioldstart+sgnold))-log(xold(ioldstart))))
              else
                 fnew(inew) = 0.d0
              endif
           else
              fnew(inew) = fold(ioldstart)
           endif
        elseif(elow.eq.3) then
           fnew(inew) = 0.d0
        elseif(elow.eq.4) then
           if(nold.gt.1) then
              eps = ((xnew(inew)-xold(ioldstart))/                     &
                         (xold(ioldstart+sgnold)-xold(ioldstart)))
              if(eps.gt.-1.d0) then
                 fnew(inew) = (1-abs(eps))*fold(ioldstart)
              else
                 fnew(inew) = 0.d0
              endif
           else
              fnew(inew) = fold(ioldstart)
           endif
        else
           stop 1349
        endif
     elseif(xnew(inew).gt.xold(ioldend)) then
        !
        ! higher than highest boundary 
        !
        if(eup.eq.0) then
           ierror = -1
           return
        elseif(eup.eq.1) then
           fnew(inew) = fold(ioldend)
        elseif(eup.eq.2) then
           if(nold.gt.1) then
              if((fold(ioldend-sgnold).gt.0.d0).and. &
                   (fold(ioldend).gt.0.d0)) then
                 fnew(inew) = fold(ioldend) *                               &
                        (fold(ioldend-sgnold)/fold(ioldend))**              &
                        ((log(xnew(inew))-log(xold(ioldend)))/              &
                         (log(xold(ioldend-sgnold))-log(xold(ioldend))))
              else
                 fnew(inew) = 0.d0
              endif
           else
              fnew(inew) = fold(ioldend)
           endif
        elseif(eup.eq.3) then
           fnew(inew) = 0.d0
        elseif(eup.eq.4) then
           if(nold.gt.1) then
              eps = ((xnew(inew)-xold(ioldend))/                     &
                         (xold(ioldend-sgnold)-xold(ioldend)))
              if(eps.gt.-1.d0) then
                 fnew(inew) = (1-abs(eps))*fold(ioldend)
              else
                 fnew(inew) = 0.d0
              endif
           else
              fnew(inew) = fold(ioldend)
           endif
        else
           stop 1349
        endif
     else
        !
        ! Within the domain of xold
        !
        if(av.eq.0) then
           !
           ! Simple remapping: just interpolation
           !
           call hunt(xold,nold,xnew(inew),iold)
           if(xnew(inew).eq.xold(1)) then
              fnew(inew) = fold(1)
           elseif(xnew(inew).eq.xold(nold)) then
              fnew(inew) = fold(nold)
           else
              if((iold.le.0).or.(iold.ge.nold)) stop 1350
              eps = (xnew(inew)-xold(iold)) / (xold(iold+1)-xold(iold))
              if((fold(iold).gt.0.d0).and.(fold(iold+1).gt.0.d0)) then
                 fnew(inew) = exp((1.d0-eps)*log(fold(iold))+eps*log(fold(iold+1)))
              else
                 fnew(inew) = (1.d0-eps)*fold(iold)+eps*fold(iold+1)
              endif
           endif
        elseif(av.eq.1) then
           !
           ! Remapping with assurance that the integral over the function
           ! remains roughly OK even if much lower resolution is used
           !
           write(stdo,*) 'SORRY: The integral-conserving remapping of functions'
           write(stdo,*) '       is not yet implemented...'
           stop 4722
        else
           stop 4721
        endif
     endif
  enddo
  !
end subroutine remap_function


!--------------------------------------------------------------------------
!                     PREPARE TABULATED EXPONENTIAL
!--------------------------------------------------------------------------
subroutine prep_tabulated_exp(x0,x1,nx)
  implicit none
  double precision :: x0,x1
  integer :: nx,ix,ierr
  exp_nx = nx
  if(allocated(exp_x)) deallocate(exp_x)
  if(allocated(exp_y)) deallocate(exp_y)
  allocate(exp_x(nx),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate tabulated exponential array'
     stop
  endif
  allocate(exp_y(nx),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR: Could not allocate tabulated exponential array'
     stop
  endif
  do ix=1,nx
     exp_x(ix) = x0 - (x1-x0)*((ix-1.d0)/(nx-1.d0))
     exp_y(ix) = exp(exp_x(ix))
  enddo
  exp_fact = (exp_nx-1)/(exp_x(exp_nx)-exp_x(1))
end subroutine prep_tabulated_exp

!--------------------------------------------------------------------------
!                   THE EXPONENTIAL WITH LOOKUP
!--------------------------------------------------------------------------
function tabexp(x)
  implicit none
  double precision :: x,tabexp,eps
  integer :: ix
  if(exp_nx.eq.0) then
     tabexp = exp(x)
     return
  else
     if(x.le.exp_x(1)) then
        tabexp = exp(x)
        return
     endif
     if(x.ge.exp_x(exp_nx)) then
        tabexp = exp(x)
        return
     endif
     ix  = exp_fact*(x-exp_x(1)) + 1
     if((ix.lt.1).or.(ix.ge.exp_nx)) stop 3234
     eps = (x-exp_x(ix))/(exp_x(ix+1)-exp_x(ix))
     if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 3235
     tabexp = (1.d0-eps)*exp_y(ix)+eps*exp_y(ix+1)
     return
  endif
end function tabexp

!--------------------------------------------------------------------------
!                     DEALLOCATE TABULATED EXPONENTIALS
!--------------------------------------------------------------------------
subroutine tabulated_exp_cleanup()
  implicit none
  if(allocated(exp_x)) deallocate(exp_x)
  if(allocated(exp_y)) deallocate(exp_y)
end subroutine tabulated_exp_cleanup



!--------------------------------------------------------------------------
!                       NUMERICAL RECIPES ROUTINES
!--------------------------------------------------------------------------
FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  doubleprecision ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  !$OMP THREADPRIVATE(iv,iy,idum2)
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
END FUNCTION ran2
! C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE hunt(xx,n,x,jlo)
  INTEGER jlo,n
  doubleprecision x,xx(n)
  INTEGER inc,jhi,jm
  LOGICAL ascnd
  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt
! C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  double precision a(np,np),b(n)
  INTEGER i,ii,j,ll
  double precision sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        enddo
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  enddo
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     enddo
     b(i)=sum/a(i,i)
  enddo
  return
END SUBROUTINE lubksb
! (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE ludcmp(a,n,np,indx,d,success)
  INTEGER n,np,indx(n),NMAX
  double precision d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  double precision aamax,dum,sum,vv(NMAX)
  logical,optional :: success
  if(present(success)) success = .true.
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     enddo
     if (aamax.eq.0.) then
        if(present(success)) then
           success=.false.
           return
        else
           stop 'singular matrix in ludcmp'
        endif
     endif
     vv(i)=1./aamax
  enddo
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
     enddo
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     enddo
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        enddo
     endif
  enddo
  return
END SUBROUTINE ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

end module mathroutines_module
