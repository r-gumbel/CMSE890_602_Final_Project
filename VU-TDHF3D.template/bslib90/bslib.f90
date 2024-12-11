module spline_grid
! subroutine grid is created as a module since various spline dimensions are
! not known in the main program. use of pointers take care of this issue.
      contains
      subroutine grid(coefl, coefr, xkv, xcl, xin, xout, nnot, ncol,&
                      nord, nord1, nspl, igrd, ianlyz, iper)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: nnot
      integer , intent(out) :: ncol
      integer , intent(out) :: nord
      integer , intent(out) :: nord1
      integer , intent(out) :: nspl
      integer , intent(out) :: igrd
      integer , intent(out) :: ianlyz
      integer , intent(out) :: iper
      real(wp) , intent(inout) :: xin
      real(wp) , intent(inout) :: xout
      real(wp) , dimension(:,:), pointer :: coefl
      real(wp) , dimension(:,:), pointer :: coefr
      real(wp) , dimension(:), pointer   :: xkv
      real(wp) , dimension(:), pointer   :: xcl
      real(wp) , dimension(:), allocatable   :: xkv_local
      real(wp) , dimension(:), allocatable   :: xg
!
      integer , parameter :: ncolmax = 256
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, nbc, n, npsum, np, mode, iu
      real(wp) :: xl, xu, fn, powdi, dx, dxtoi, xxin, xxout, delx, a, b, yi
!----------------------------------------------------------------------
!    calculate the knots xkv
!
!        xkv(nord)      = in boundary point
!        xkv(nord+ncol) = out boundary point
!
!        equally spaced knots for periodic boundary conditions.
!        multiple knots at the boundary for fixed bc.
!
!        nord : order of splines
!        if igrd < 0 print out the knots,collocation points and stop
!        if igrd > 0 print out the knots,collocation points
!        if ianlyz=1 print out all derivative matrices etc.
!        iper = 0 fixed boundary conditions
!        iper = 1 periodic boundary conditions
!----------------------------------------------------------------------
      read (3, *) nord, igrd, ianlyz, iper
!
      if (associated(xkv)) nullify(xkv)
      if (associated(xcl)) nullify(xcl)
      if (associated(coefl)) nullify(coefl)
      if (associated(coefr)) nullify(coefr)
!
      allocate ( coefl(nord,nord), coefr(nord,nord) )
!
      do i = 1, nord
         do j = 1, nord
            coefl(i,j) = 0.0_wp
            coefr(i,j) = 0.0_wp
         end do
      end do
!
      if (iper == 0) then
         nbc = (nord - 1)/2
         do i = 1, nbc
            read (3, *) (coefl(i,j),j=1,nord - 1)
         end do
         do i = 1, nbc
            read (3, *) (coefr(i,j),j=1,nord - 1)
         end do
      endif
!
      nord1 = nord + 1
      if (mod(nord,2) == 0) then
         write (8, '(2/,A,2/)') ' spline order must be odd'
         stop
      endif
      if (iper/=0 .and. iper/=1) then
         write (8, '(2/,A,2/)') ' incorrect value of iper entered'
         stop
      endif
!----------------------------------------------------------------------
!       read the coefficients determining the type of boundary
!       condition for iper=0 (fixed bc) option. for periodic bc
!       we still read them in but they are not used. the total
!       number of boundary conditions must be (nord-1). these are
!       divided equally between the left and right boundaries. the
!       first (nord-1)/2 lines read from input represent the bc for
!       the left boundary and the next (nord-1)/2 lines for the
!       right boundary. each line of input should contain (nord-1)
!       flags (0 or 1) indicating the type of bc:
!
!       example:  nord=5 and we would like to impose the condition
!                 that function vanishes on left boundary and its
!                 derivative vanishes on the right boundary.
!                 the input will look like:
!
!                 1 0 0 0
!                 0 0 0 1
!                 0 1 0 0
!                 0 0 0 1
!
!                 note that on each line a nonzero value will enforce
!                 the bc that either the function or some of its
!                 derivatives are zero. the first flag is for the
!                 function itself, the second flag for the first
!                 derivative etc. the second and fourth lines make
!                 the third derivative zero. these are just added so
!                 that the total number of conditions is (nord-1).
!                 they have almost no influence on the result.
!----------------------------------------------------------------------
!
!       read  parameters defining the spline lattice.
!       input is four numbers per line, each line generates a
!       region of knots as follows:
!
!       xl(i),xu(i),np(i),mode(i) in free format
!
!                     xl(i)   -  lower limit of region i
!                     xu(i)   -  upper limit of region i
!                     np(i)   -  no. of collocation points in region i
!                     mode(i) -  mode for generating the
!                                points in the  region i
!                                mode = 0  arith.  series
!                                mode = 1  geom.   series
!                                mode < 0  terminate the input
!
!        note 1: when we enter multiple knot intervals the upper limit
!                of the preceeding interval should be the same as the
!                lower limit of the next interval.
!        note 2: np is actually the number of collocation points. so
!                the total number of np's should be less than the
!                parameter ncolmax defined here (check is provided)
!------------------------------------------------------------------------
      allocate (xkv_local(ncolmax), xg(ncolmax) )
!
      n = nord
      npsum = 0
   30 continue
      read (3, *) xl, xu, np, mode
!
      npsum = npsum + np
      if (npsum > ncolmax) then
         write (8, '(2/,A,i4,2/)') ' ncol exceeds maximum allowed value of',ncolmax
         stop
      endif
      if (mode < 0) go to 60
      xout = xu
      if (n == nord) xin = xl
      if (iper==1 .and. mode/=0) then
         write (8, '(2/,A,2/)') ' periodic bc only for equal mesh'
         stop
      endif
!
      select case (mode)
!
      case (0)
         dx = (xu - xl)/np
         do i = 0, np
            xkv_local(n+i) = xl + dx*i
         end do
         n = n + np
!
      case (1)
         fn = np
         powdi = 1.0_wp/fn
         dx = (xu/xl)**powdi
         dxtoi = 1.0_wp/dx
         do i = 0, np
            dxtoi = dxtoi*dx
            xkv_local(n+i) = xl*dxtoi
         end do
         n = n + np
!
      case (2)
         xxin = 0.0_wp
         xxout = 1.0_wp
         delx = (xxout - xxin)/np
         do i = 0, np
            xg(i+1) = xxin + i*delx
         end do
         a = 1.0_wp
         b = 3.3_wp
         do i = 0, np
            yi = xg(i+1)
            xkv_local(n+i) = a*(exp(b*yi) - 1.0_wp)
         end do
         n = n + np
!
      case (3)
         xxin = 0.0_wp
         xxout = 1.0_wp
         delx = 2.0_wp*(xxout - xxin)/np
         do i = 1, np/2
            xg(i) = xxin + i*delx
         end do
         b = log(xu)
         do i = 1, np/2
            yi = xg(i)
            xkv_local(n+np/2+i) = exp(b*yi)
         end do
         do i = 1, np/2
            xkv_local(n+np/2-i) = - xkv_local(n+np/2+i)
         end do
         xkv_local(n+np/2) = 0.0_wp
         n = n + np
!
      end select
!
      go to 30
   60 continue
      ncol = n - nord
      nnot = ncol + 2*nord - 1
      nspl = ncol + nord - 1

      allocate (xkv(nnot), xcl(ncol) )
      xkv = xkv_local(1:nnot)
!----------------------------------------------------------------------
!    calculate the additional knots needed for the generation of
!    splines. for fixed boundary conditions we need multiple knots
!    at the boundary. this is required for the proper calculation of
!    the areas under the splines at the boundaries. for periodic
!    splines (equally spaced mesh) we simply extend the splines
!    using the constant spacings of the first and last knot intervals
!----------------------------------------------------------------------
      if (iper == 0) then
         do i = 1, nord
            xkv(i) = xin
            xkv(ncol+nord+i-1) = xout
         end do
      else
         dx = xkv(nord+1) - xkv(nord)
         do i = nord - 1, 1, -1
            xkv(i) = xkv(i+1) - dx
         end do
         dx = xkv(nord+ncol) - xkv(nord+ncol-1)
         do i = 1, nord - 1
            xkv(nord+ncol+i) = xkv(nord+ncol+i-1) + dx
         end do
      endif
!----------------------------------------------------------------------
!    calculate the collocation points xcl, each falling between two
!    data points starting at nord'th knot point
!----------------------------------------------------------------------
      do i = 1, ncol
         iu = nord - 1 + i
         xcl(i) = (xkv(iu)+xkv(iu+1))/2.0_wp
      end do
!
      deallocate ( xkv_local, xg)
!
      return
      end subroutine grid
!
end module spline_grid
      subroutine spline(der1, der2, binvp, xcl, xkv, coefl, coefr, wx, xin, xout, ncol, nord, nnot,&
                        nspl, igrd, ianlyz, iper, cname)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)   :: ncol
      integer, intent(in)   :: nord
      integer, intent(in)   :: nnot
      integer, intent(in)   :: nspl
      integer, intent(in)   :: igrd
      integer, intent(in)   :: ianlyz
      integer, intent(in)   :: iper
      real(wp), intent(in)  :: xcl(ncol)
      real(wp), intent(in)  :: xkv(nnot)
      real(wp), intent(in)  :: coefl(nord,nord), coefr(nord,nord)
      real(wp), intent(in)  :: xin
      real(wp), intent(in)  :: xout
      real(wp), intent(out) :: der1(ncol,ncol)
      real(wp), intent(out) :: der2(ncol,ncol)
      real(wp), intent(out) :: wx(ncol)
      character, intent(in) :: cname
      Real(wp) :: binvp(ncol,ncol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(:),   allocatable  :: indx0, iwk
      real(wp), dimension(:),   allocatable  :: ekin, ws1, wkr
      real(wp), dimension(:,:), allocatable  :: btil, btilm
      real(wp), dimension(:,:), allocatable  :: beta
      real(wp), dimension(:,:), allocatable  :: tmp
      real(wp), dimension(:,:), allocatable  :: xkbl, xkbr
      real(wp), dimension(:,:), allocatable  :: zkin
      real(wp), dimension(:,:), allocatable  :: cgwk, wsk1, wsk2
      real(wp), dimension(:,:), allocatable  :: bwrk
      integer   :: nord1
      integer   :: i
! set some extra dimensions
      nord1 = nord + 1
!-----------------------------------------------------------------------
!        allocate arrays for the spline routines
!-----------------------------------------------------------------------
      allocate( indx0(ncol), iwk(ncol) )
      allocate( ekin(ncol), ws1(ncol), wkr(ncol) )
      allocate( btil(nspl,nspl), btilm(nspl,nspl) )
      allocate( beta(nord,nord), xkbl(nord,nord), xkbr(nord,nord) )
      allocate( zkin(ncol,ncol) )
      allocate( cgwk(ncol,ncol) )
      allocate( wsk1(ncol,ncol), wsk2(ncol,ncol) )
      allocate( bwrk(nord,ncol) )
      allocate( tmp(nord,nord1) )
!----------------------------------------------------------------------
!       this routine calculates all spline quantities. the details
!       of calculations are described in the subroutines the returned
!       quantities are (for dimensions see the statement):
!
!        wx     :   integration weigts
!        der1   :   first derivative matrix
!        der2   :   second derivative matrix (direct d2 or d*d)
!        zkin   :   eigenfunctions of -0.5*d2
!        ekin   :   eigenvalues of -0.5*d2
!
!       active input quantities are:
!
!        nnot   :   number of knots
!        ncol   :   number of collocation points
!        nspl   :   number of splines for fixed boundary conditions
!        xkv    :   the spline knots
!        xcl    :   the collocation points
!        and all flags
!
!       everything else is for work space
!----------------------------------------------------------------------
      if (iper == 1) then
         call setspp (bwrk, cgwk, der1, der2, zkin, binvp, ekin, xkv, xcl, tmp,&
                      wsk1, wsk2, ws1, wkr, iwk, wx, indx0, nnot, ncol, nspl, nord,&
                      nord1)
      else if (iper == 0) then
         call setspb (bwrk, cgwk, btil, btilm, der1, der2, zkin, ekin, wsk1, wsk2, &
                      ws1, wkr, iwk, xkv, xcl, coefr, coefl, beta, tmp, xkbl, xkbr,&
                      wx, indx0, xin, xout, nnot, ncol, nspl, nord, nord1)
      else
         write (8, '(2/,A,2/)') ' wrong value of iper entered'
         stop
      endif
!----------------------------------------------------------------------
!        print grid information and weights
!        if detailed analsis is required print all matrices (ianlyz=1)
!----------------------------------------------------------------------
      if (cname /= ' ') then
         if (igrd /= 0) then
            write (8, '(2/,1X,A,A,2/)') cname, '-lattice'
            call putgrd (xkv, xcl, wx, nnot, ncol, nord, igrd)
         endif
!
         if (ianlyz == 1) then
            write (8, '(A,/)') ' first derivative matrix:'
            call putmat (der1, ncol)
            write (8, '(/,A,/)') ' second derivative matrix:'
            call putmat (der2, ncol)
            write (8, '(/,A,/)') ' kinetic energy (-0.5*d2) eigenvalues:'
            write (8, '(1P,10E12.4)') (ekin(i),i=1,ncol)
         endif
      endif
!----------------------------------------------------------------------
!        test overlaps of eigenfunctions
!----------------------------------------------------------------------
!       do 45 j=1,ncol
!       do 42 k=1,ncol
!          sumo=0.0
!       do 41 l=1,ncol
!          sumo=sumo+wx(l)*zkin(l,j)*zkin(l,k)
!41     continue
!          write(*,*) 'states, overlap:',j,k,sumo
!42     continue
!45     continue
!
      deallocate( indx0, iwk, ekin, ws1, wkr, btil, btilm, beta, xkbl, xkbr )
      deallocate( zkin, cgwk, wsk1, wsk2, bwrk, tmp )
!
      return
      end subroutine spline
      subroutine bcset(beta, xkbl, xkbr, coefr, coefl, xkv, xin, xout, tmp,&
                       indxl, indxr, nord, nord1, nnot, nspl)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: indxl
      integer , intent(out) :: indxr
      integer  :: nord
      integer  :: nord1
      integer  :: nnot
      integer  :: nspl
      real(wp)  :: xin
      real(wp)  :: xout
      real(wp) , intent(inout) :: beta(nord,nord)
      real(wp)  :: xkbl(nord,nord)
      real(wp)  :: xkbr(nord,nord)
      real(wp) , intent(in) :: coefr(nord,nord)
      real(wp) , intent(in) :: coefl(nord,nord)
      real(wp)  :: xkv(nnot)
      real(wp)  :: tmp(nord,nord1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ii, nbcl, nbcr, ik, ir, index, irb
!----------------------------------------------------------------------
!      beta : boundary sub-matrix (output)
!      xkbl : work array
!      xkbr : work array
!      coefr: coefficients for selecting the right boundary cond.
!      coefl: coefficients for selecting the left boundary cond.
!      xkv  : array containing the knot values
!      xin  : left boundary point value
!      xout : right boundary point value
!      tmp  : work array used in splgen
!      nord : order of splines
!      nord1: nord+1
!      nnot : number of knots (ncol+2*nord-1)
!      nspl : number of splines (ncol+nord-1)
!
!      comments:
!      routine to calculate the boundary sub-matrix which is appended
!      to the b-matrix before inversion. the form of this matrix
!      is:
!              beta(is,ir)=sum(k) coef(ir,ik)*b**(ik) (is)
!      where 'is' is the number of splines (ncol+nord-1) and 'b' is
!      evaluted at boundaries specified by 'xin' and 'xout'.
!      the values of ik 1,..,nord-1. ik=1 is for setting the function
!      itself to zero, ik=2 is for setting the first derivative to
!      zero etc. we have to specify a total of nord-1 boundary
!      conditions (left and right summed up).
!
!      the selection of particular boundary conditions is done by
!      the arrays coefl(ir,ik) and coefr(ir,ik). in each case
!      ir=1,... counts the number of boundary conditions and
!      ik=1,...,nord-1 specifies the type of boundary condition as
!      mentioned above. for examples coefl(1,1)=1.0 indicates a
!      a vanishing function at the left boundary, coefl(1,2)=1.0
!      indicates a vanishing derivative at the left boundary.
!      coefl(1,1)=1.0 and coefl(2,2)=1.0 indicates two boundary
!      conditions; vanishing function and vanishing derivative
!      on the left boundary. the total number of boundary conditions
!      on left and right must be nord-1. so to impose more number of
!      boundary conditions we need to increase the spline order.
!----------------------------------------------------------------------
      do i = 1, nord
         do ii = 1, nord
            xkbl(ii,i) = 0.0_wp
            xkbr(ii,i) = 0.0_wp
         end do
      end do
!----------------------------------------------------------------------
!      construct splines (or derivatives) evaluated at left and
!      right boundary points using also the coefficient arrays
!      to determine the non-vanishing coefficients in the expression
!      for the boundary sub-matrix beta.
!      nbcl and nbcr are the number of left and right boundary cond.
!----------------------------------------------------------------------
      nbcl = (nord - 1)/2
      nbcr = (nord - 1)/2
!
      l3: do ik = 1, nord
         do ir = 1, nord
            if (coefl(ir,ik) == 0.0_wp) cycle
            call bndc (ik - 1, xkv, xin, xkbl(1,ik), tmp, index, nord, nord1, nnot, nspl)
            cycle  l3
         end do
      end do l3
      indxl = index
!
      l5: do ik = 1, nord
         do ir = 1, nord
            if (coefr(ir,ik) == 0.0_wp) cycle
            call bndc (ik - 1, xkv, xout, xkbr(1,ik), tmp, index, nord, nord1, nnot, nspl)
            cycle  l5
         end do
      end do l5
      indxr = index
!----------------------------------------------------------------------
!      construct the boundary sub-matrix beta by summing over all k
!      values. the left and right results are added on top of one
!      another. total number of boundary conditions (nbcl+nbcr) must
!      be (nord-1) in order to create a square matrix when added to
!      the non-square b-matrix.
!-----------------------------------------------------------------------
      do ir = 1, nbcl
         do ii = 1, nord
            beta(ii,ir) = 0.0_wp
            do ik = 1, nord
               beta(ii,ir) = beta(ii,ir) + coefl(ir,ik)*xkbl(ii,ik)
            end do
         end do
      end do
!
      do ir = 1, nbcr
         irb = ir + nbcl
         do ii = 1, nord
            beta(ii,irb) = 0.0_wp
            do ik = 1, nord
               beta(ii,irb) = beta(ii,irb) + coefr(ir,ik)*xkbr(ii,ik)
            end do
         end do
      end do
!
      return
      end subroutine bcset
      subroutine bndc(ip, xkv, xbnd, bbnd, tmp, index, nord, nord1, nnot, nspl)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer  :: ip
      integer  :: index
      integer  :: nord
      integer  :: nord1
      integer  :: nnot
      integer  :: nspl
      real(wp)  :: xbnd
      real(wp)  :: xkv(nnot)
      real(wp)  :: bbnd(nord)
      real(wp)  :: tmp(nord,nord1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii
!----------------------------------------------------------------------
!      routine to call for the calculation of b-splines at the
!      boundary points (xbnd). ip has the same meaning as in splgen
!      and calulates derivatives of b-splines. what is returned are
!      the non-zero spline values at point x=xbnd. these points
!      are used in calculating the submatrix for boundary conditions
!----------------------------------------------------------------------
      do ii = 1, nord
         bbnd(ii) = 0.0_wp
      end do
!
      call splgens (ip, xkv, xbnd, bbnd, tmp, index, nord, nord1, nnot, nspl)
!
      return
      end subroutine bndc
      subroutine findi(xkv, nnot, nord, nspl, tcol, i)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nnot
      integer , intent(in) :: nord
      integer , intent(in) :: nspl
      integer , intent(out) :: i
      real(wp) , intent(in) :: tcol
      real(wp) , intent(in) :: xkv(nnot)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij, i1, i2
!----------------------------------------------------------------------
!      calculates the lower index of the knot interval containing the
!      the point tcol. the algorithm is primitive but it also works
!      for multiply defined boundary knots.
!----------------------------------------------------------------------
      i = 0
      do ij = 1, nnot
         if (xkv(ij) <= tcol) cycle
         i = ij - 1
         exit
      end do
      i1 = i
      i = 0
      do ij = nnot, 1, -1
         if (xkv(ij) >= tcol) cycle
         i = ij
         exit
      end do
      i2 = i
      i = max0(i1,i2)
      i = max0(nord,i)
      i = min0(nspl,i)
      return
      end subroutine findi
      subroutine invtst(a, b, n)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      real(wp) , intent(in) :: a(n,n)
      real(wp) , intent(in) :: b(n,n)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(wp), parameter :: epsinv = 0.00000001_wp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k
      real(wp) :: sumt
!----------------------------------------------------------------------
!      tests the inverse of a matrix by comparing the product to an
!      identity matrix.
!----------------------------------------------------------------------
      sumt = 0.0_wp
      do i = 1, n
         do j = 1, n
            do k = 1, n
               sumt = sumt + a(i,k)*b(k,j)
            end do
         end do
      end do
!
      sumt = sumt - 1.0_wp*n
      sumt = sumt/(n*n)
      if (abs(sumt) > epsinv) then
         write (8, '(2/,A,/)') ' invtst: inversion accuracy problem'
         write (8, '(A,1P,E12.5,A,E12.5,/)') ' accuracy=', sumt, '  eps=', &
                                                           epsinv
      endif
!
      return
      end subroutine invtst
      subroutine jacobi(a, n, np, d, v)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: np
      real(wp) , intent(inout) :: a(np,np)
      real(wp) , intent(inout) :: d(np)
      real(wp) , intent(inout) :: v(np,np)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nmax = 1000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ip, iq, i, j
      real(wp), dimension(nmax) :: b, z
      real(wp) :: sm, tresh, g, h, t, theta, c, s, tau
!-----------------------------------------------
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0_wp
         end do
         v(ip,ip) = 1.0_wp
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0_wp
      end do
      do i = 1, 50
         sm = 0.0_wp
         do ip = 1, n - 1
            do iq = ip + 1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm == 0.0_wp) return
         if (i < 4) then
            tresh = 0.2_wp*sm/n**2
         else
            tresh = 0.0_wp
         endif
         do ip = 1, n - 1
            do iq = ip + 1, n
               g = 100.0_wp*abs(a(ip,iq))
               if (i>4 .and. abs(d(ip))+g==abs(d(ip)) .and. abs(d(iq))+g==abs(d(iq))) then
                  a(ip,iq) = 0.0_wp
               else if (abs(a(ip,iq)) > tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h) + g == abs(h)) then
                     t = a(ip,iq)/h
                  else
                     theta = 0.5_wp*h/a(ip,iq)
                     t = 1.0_wp/(abs(theta) + sqrt(1.0_wp + theta**2))
                     if (theta < 0.0_wp) t = -t
                  endif
                  c = 1.0_wp/sqrt(1.0_wp + t**2)
                  s = t*c
                  tau = s/(1.0_wp + c)
                  h = t*a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0_wp
                  do j = 1, ip - 1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h + g*tau)
                     a(j,iq) = h + s*(g - h*tau)
                  end do
                  do j = ip + 1, iq - 1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h + g*tau)
                     a(j,iq) = h + s*(g - h*tau)
                  end do
                  do j = iq + 1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h + g*tau)
                     a(iq,j) = h + s*(g - h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h + g*tau)
                     v(j,iq) = h + s*(g - h*tau)
                  end do
               endif
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0_wp
         end do
      end do
      write (*, *) 'JACOBI: 50 iterations should never happen'
      return
      end subroutine jacobi
      subroutine matq(a, x, nr, nv, na, nx)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nr
      integer , intent(in) :: nv
      integer , intent(in) :: na
      integer , intent(in) :: nx
      real(wp) , intent(inout) :: a(na*na)
      real(wp) , intent(inout) :: x(na*nx)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, kkna, mx, mr, nr1, kr1, nvx, nvxr, ma, na1, nrnr1, nrnr, &
                 nrna1, nrnak, nvxk, kk, ik, iprj, iprjpr, kj, ij, nrj, ii, il, lj
      real(wp) :: pivot, z, sumt, dtest
!----------------------------------------------------------------------
!     a subroutine for solving systems of linear equations or for
!     inverting a matrix (real quantities only).
!     the input parameters are as follows:
!         a  =  real na x na matrix containing the coefficients of the
!               linear system. in the case of the matrix inversion
!               a is the matrix to be inverted.
!         x  =  for solving a single set of linear equations x is a
!               vector of length na. for solving more than one set
!               of equations x is a matrix dimensioned na x nx, and
!               for matrix inversion it should be set to a nr x nr
!               unit matrix. on output x contains either the coeffs.
!               of the linear system or the inverse matrix.
!         nr =  number of rows or the number of equations to be
!               solved.
!         nv =  number of right hand side vectors. nv=1 for a single
!               linear system and nv=nr for inverting a nr x nr matrix
!         na =  the actual dimension of a and x in the calling program
!               i.e. na.ge.nr where nr is the dimension for the
!               subsystem.
!         nx =  1 for a single system of linear equations or =na for
!               matrix inversion.
!     examples:
!         a) linear system of equations
!            call matq(a,x,nr,1,na,1)
!         b) matrix inversion
!            call matq(a,x,nr,nr,na,na)
!
!         more details can be realized by examining the dimension
!         statement in the program. on output a is destroyed.
!-----------------------------------------------------------------------
      k = 0
      kkna = k
      mx = nx
      mr = nr
      nr1 = mr - 1
      kr1 = nr1
      nvx = (nv - 1)*mx
      nvxr = nvx + mr
      ma = na
      na1 = ma + 1
      nrnr1 = na1*nr1
      nrnr = nrnr1 + 1
      nrna1 = nrnr - nr1
      do nrnak = nrna1, nrnr1
         pivot = 0.0_wp
         k = k + 1
         nvxk = k + nvx
         kkna = kkna + na1
         kk = kkna - ma
         do ik = kk, kr1
            z = -a(ik+1)
            if (z >= 0.0_wp) then
               if (z == 0.0_wp) cycle
               z = -z
            endif
            if ((-z) - pivot <= 0.0_wp) cycle
            pivot = -z
            iprj = ik
         end do
         sumt = -a(kk)
         if (sumt + pivot >= 0.0_wp) then
            if (sumt + pivot <= 0.0_wp) then
               if (sumt /= 0.0_wp) go to 11
               return
            endif
            if (pivot - sumt > 0.0_wp) then
               iprjpr = iprj - kk + k
               sumt = -a(iprj+1)
               do kj = kk, nrnak, ma
                  z = a(iprj+1)
                  a(iprj+1) = a(kj)
                  a(kj) = z
                  iprj = iprj + ma
               end do
               do kj = k, nvxk, mx
                  z = x(iprjpr+1)
                  x(iprjpr+1) = x(kj)
                  x(kj) = z
                  iprjpr = iprjpr + mx
               end do
            endif
         endif
   11    continue
         pivot = (-1.0_wp)/sumt
         do kj = k, nvxk, mx
            ij = kj
            x(kj) = pivot*x(kj)
            do ik = kk, kr1
               x(ij+1) = (-a(ik+1)*x(kj)) + x(ij+1)
               ij = ij + 1
            end do
         end do
         do kj = kkna, nrnak, ma
            ij = kj
            a(kj) = pivot*a(kj)
            do ik = kk, kr1
               a(ij+1) = (-a(kj))*a(ik+1) + a(ij+1)
               ij = ij + 1
            end do
         end do
         kr1 = kr1 + ma
      end do
      dtest = abs(a(nrnr))
      if (dtest < 1.0E-29_wp) write (8, '(2/,A,2/)') &
         ' in matq: determinant almost zero'
      pivot = 1.0_wp/a(nrnr)
      do nrj = mr, nvxr, mx
         x(nrj) = pivot*x(nrj)
         ij = nrj + 1
         ii = nrnr + ma
         do k = 1, nr1
            ii = ii - na1
            il = ii
            ij = ij - 1
            do lj = ij, nrj
               x(ij-1) = (-a(il)*x(lj)) + x(ij-1)
               il = il + ma
            end do
         end do
      end do
      return
      end subroutine matq
      subroutine ordeig(nm, n, x, wr, wi, zr, zi)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nm
      integer , intent(in) :: n
      real(wp) , intent(inout) :: x(n)
      real(wp) , intent(inout) :: wr(n)
      real(wp) , intent(inout) :: wi(n)
      real(wp) , intent(inout) :: zr(nm,n)
      real(wp) , intent(inout) :: zi(nm,n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, is, i, j
      real(wp) :: s
!-----------------------------------------------
!...   orders w,z according to array  x
      do m = 1, n - 1
         s = x(m)
         is = m
         do i = m + 1, n
            if (x(i) > s) cycle
            s = x(i)
            is = i
         end do
         if (is == m) cycle
         x(is) = x(m)
         x(m) = s
         s = wr(is)
         wr(is) = wr(m)
         wr(m) = s
         s = wi(is)
         wi(is) = wi(m)
         wi(m) = s
         do j = 1, n
            s = zr(j,is)
            zr(j,is) = zr(j,m)
            zr(j,m) = s
            s = zi(j,is)
            zi(j,is) = zi(j,m)
            zi(j,m) = s
         end do
      end do
      return
      end subroutine ordeig
      subroutine putgrd(xkv, xcl, wq, nnot, ncol, nord, igrd)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nnot
      integer , intent(in) :: ncol
      integer , intent(in) :: nord
      integer , intent(in) :: igrd
      real(wp) , intent(in) :: xkv(nnot)
      real(wp) , intent(in) :: xcl(ncol)
      real(wp) , intent(in) :: wq(ncol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!----------------------------------------------------------------------
!      write grid information and weights
!----------------------------------------------------------------------
      if (igrd == 0) return
!
      write (8, '(/,A,I6,A,/,A,I6,A,/,A,I6,A,/)') &
         ' nord  =', nord, '    (order  of b-splines)',&
         ' ncol  =', ncol, '    (number of collocation points)', &
         ' nnot  =', nnot, '    (number of knots)'
      write (8, '(A,/)') '   knots:'
      write (8, '(7(1X,F10.3))') (xkv(i),i=1,nnot)
      write (8, '(/,A,/)') '   collocation points:'
      write (8, '(7(1X,F10.3))') (xcl(i),i=1,ncol)
      write (8, '(/,A,/)') '   weights:'
      write (8, '(7(1X,F10.3))') (wq(i),i=1,ncol)
      write (8, '(2/)')
!
      if (igrd < 0) then
         write (8, '(2/,A,2/)') ' end of grid examination'
         stop
      endif
      return
      end subroutine putgrd
      subroutine putmat(dmat, ncol)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ncol
      real(wp) , intent(in) :: dmat(ncol,ncol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!----------------------------------------------------------------------
!      write the matrix dmat in collocation space
!----------------------------------------------------------------------
      do i = 1, ncol
         write (8, '(1P,10E12.4)') (dmat(i,j),j=1,ncol)
      end do
      return
      end subroutine putmat
      subroutine setspb(bwrk, cgwk, btil, btilm, der1, der2, zkin, ekin, wsk1, wsk2, ws1,&
                        wkr, iwk, xkv, xdv, coefr, coefl, beta, tmp, xkbl, xkbr, wq,&
                        indx0, xin, xout, nnot, ncol, nspl, nord, nord1)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer  :: nnot
      integer  :: ncol
      integer  :: nspl
      integer  :: nord
      integer  :: nord1
      real(wp)  :: xin
      real(wp)  :: xout
      integer , intent(out) :: iwk(ncol)
      integer  :: indx0(ncol)
      real(wp)  :: bwrk(nord,ncol)
      real(wp)  :: cgwk(ncol,ncol)
      real(wp)  :: btil(nspl,nspl)
      real(wp)  :: btilm(nspl,nspl)
      real(wp) , intent(inout) :: der1(ncol,ncol)
      real(wp) , intent(inout) :: der2(ncol,ncol)
      real(wp)  :: zkin(ncol,ncol)
      real(wp)  :: ekin(ncol)
      real(wp)  :: wsk1(ncol,ncol)
      real(wp) , intent(out) :: wsk2(ncol,ncol)
      real(wp)  :: ws1(ncol)
      real(wp)  :: wkr(ncol)
      real(wp)  :: xkv(nnot)
      real(wp)  :: xdv(ncol)
      real(wp)  :: coefr(nord,nord)
      real(wp)  :: coefl(nord,nord)
      real(wp)  :: beta(nord,nord)
      real(wp)  :: tmp(nord,nord1)
      real(wp)  :: xkbl(nord,nord)
      real(wp)  :: xkbr(nord,nord)
      real(wp) , intent(inout) :: wq(ncol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nbcl, nbcr, i, j, indxl, indxr, ia, jj, k
      real(wp) :: dx, sumt
!----------------------------------------------------------------------
!       subroutine to calculate the following quantities for fixed
!       boundary conditions:
!
!      der1   :   first derivative operator in collocation space
!      der2   :   second derivative operator in collocation space
!      zkin   :   eigenvectors of the -0.5*der2
!      ekin   :   eigenvalues of -0.5*der2
!      wq     :   integration weights in collocation space
!      btilm  :   inverse b-tilde matrix. this could be used to
!                 calculate expansion coefficients for interpolation
!      bwrk   :   array containing the non-zero splines at knots in
!                 the physical boundary
!
!      fixed boundary conditions are written as an additional set
!      of linear equations in terms of splines evaluated at the mesh
!      boundaries xin and xout. these equations when appended to the
!      non-square b-matrix result in a square matrix which is then
!      inverted and used in the calculation of operators in  the
!      collocation space. the details of the boundary condition setup
!      is given in the comments of the subroutine bcset.
!----------------------------------------------------------------------
!    begin block to compute the tilde matrix for inversion
!    call splgen to evaluate all nonzero splines at all collocation
!    points.
!----------------------------------------------------------------------
      call splgen (0, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!----------------------------------------------------------------------
!    call bcset to get the additional splines (or derivatives)
!    evaluated at the boundary points. beta is then appended to
!    the collocation matrix to form a square tilde matrix.
!----------------------------------------------------------------------
      nbcl = 0
      nbcr = 0
      do i = 1, nord
         do j = 1, nord
            if (coefl(i,j) /= 0.0_wp) nbcl = nbcl + 1
            if (coefr(i,j) == 0.0_wp) cycle
            nbcr = nbcr + 1
         end do
      end do
!
      call bcset (beta, xkbl, xkbr, coefr, coefl, xkv, xin, xout, tmp, indxl,&
                  indxr, nord, nord1, nnot, nspl)
!
      do i = 1, nspl
         do j = 1, nspl
            btilm(i,j) = 0.0_wp
            btil(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            btil(ia,i-jj+1) = bwrk(jj,ia)
         end do
      end do
!
      do ia = ncol + 1, ncol + nbcl
         i = indxl
         do j = 1, nord
            btil(ia,i-j+1) = beta(j,ia-ncol)
         end do
      end do
!
      do ia = ncol + nbcl + 1, ncol + nbcl + nbcr
         i = indxr
         do j = 1, nord
            btil(ia,i-j+1) = beta(j,ia-ncol)
         end do
      end do
!----------------------------------------------------------------------
!    invert the b-tilde matrix and store in btilm
!----------------------------------------------------------------------
      do i = 1, nspl
         btilm(i,i) = 1.0_wp
      end do
!
      call matq (btil, btilm, nspl, nspl, nspl, nspl)
!----------------------------------------------------------------------
!    since the inversion destroys btil we need to generate it again
!    in order to be able to test the inverse. this is a duplication
!    of code done because there is no other work array (nspl x nspl)
!----------------------------------------------------------------------
      do i = 1, nspl
         do j = 1, nspl
            btil(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            btil(ia,i-jj+1) = bwrk(jj,ia)
         end do
      end do
!
      do ia = ncol + 1, ncol + nbcl
         i = indxl
         do j = 1, nord
            btil(ia,i-j+1) = beta(j,ia-ncol)
         end do
      end do
!
      do ia = ncol + nbcl + 1, ncol + nbcl + nbcr
         i = indxr
         do j = 1, nord
            btil(ia,i-j+1) = beta(j,ia-ncol)
         end do
      end do
!
      call invtst (btil, btilm, nspl)
!----------------------------------------------------------------------
!     calculate integration weights and all other quantities using
!     this inverse.
!----------------------------------------------------------------------
      do i = 1, ncol
         wq(i) = 0.0_wp
         do k = 1, nspl
            dx = (xkv(k+nord)-xkv(k))/nord
            wq(i) = wq(i) + btilm(k,i)*dx
         end do
      end do
!----------------------------------------------------------------------
!    begin the calculation of the first derivative operator.
!----------------------------------------------------------------------
      call splgen (1, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      do i = 1, nspl
         do j = 1, nspl
            btil(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            btil(ia,i-jj+1) = bwrk(jj,ia)
         end do
      end do
!
      do i = 1, ncol
         do j = 1, ncol
            der1(i,j) = 0.0_wp
            do k = 1, nspl
               der1(i,j) = der1(i,j) + btil(i,k)*btilm(k,j)
            end do
         end do
      end do
!----------------------------------------------------------------------
!    begin the calculation of the second derivative operator.
!
!    note: the second derivative can be calculated two ways: 1) a
!          direct calculation from second derivative splines.
!          2) as d1*d1.
!----------------------------------------------------------------------
      call splgen (2, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      do i = 1, nspl
         do j = 1, nspl
            btil(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            btil(ia,i-jj+1) = bwrk(jj,ia)
         end do
      end do
!
      do i = 1, ncol
         do j = 1, ncol
            der2(i,j) = 0.0_wp
            do k = 1, nspl
               der2(i,j) = der2(i,j) + btil(i,k)*btilm(k,j)
            end do
         end do
      end do
!
!       do 120 i=1,ncol
!       do 120 j=1,ncol
!       der2(i,j)=0.0d0
!       do 115 k=1,ncol
!       der2(i,j)=der2(i,j)+der1(i,k)*der1(k,j)
!115    continue
!120    continue
!----------------------------------------------------------------------
!    calculate the eigenvalues and eigenvectors of the kinetic
!    energy operator
!----------------------------------------------------------------------
      do i = 1, ncol
         ekin(i) = 0.0_wp
         ws1(i) = 0.0_wp
         iwk(i) = 0
         do j = 1, ncol
            zkin(i,j) = 0.0_wp
            cgwk(i,j) = 0.0_wp
            wsk1(i,j) = -0.5_wp*der2(i,j)
            wsk2(i,j) = 0.0_wp
         end do
      end do
!
      call jacobi (wsk1, ncol, ncol, ekin, zkin)
      do i = 1, ncol
         wkr(i) = ekin(i)
      end do
!
      do j = 1, ncol
         sumt = 0.0_wp
         do i = 1, ncol
            sumt = sumt + zkin(i,j)**2
         end do
         do i = 1, ncol
            zkin(i,j) = zkin(i,j)/sqrt(sumt)
            zkin(i,j) = zkin(i,j)/sqrt(wq(i))
         end do
      end do
      call ordeig (ncol, ncol, wkr, ekin, ws1, zkin, cgwk)
!----------------------------------------------------------------------
!    refill the spline array bwrk for output
!----------------------------------------------------------------------
      call splgen (0, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      return
      end subroutine setspb
      subroutine setspp(bwrk, cgwk, der1, der2, zkin, binvp, ekin, xkv, xdv, tmp,&
                        wsk1, wsk2, ws1, wkr, iwk, wq, indx0, nnot, ncol, nspl, nord,&
                        nord1)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer  :: nnot
      integer  :: ncol
      integer  :: nspl
      integer  :: nord
      integer  :: nord1
      integer  :: indx0(ncol)
      integer  , intent(out) :: iwk(ncol)
      real(wp) , intent(out) :: wsk2(ncol,ncol)
      real(wp) , intent(inout) :: wq(ncol)
      real(wp)  :: bwrk(nord,ncol)
      real(wp)  :: cgwk(ncol,ncol)
      real(wp)  :: der1(ncol,ncol)
      real(wp)  :: der2(ncol,ncol)
      real(wp)  :: zkin(ncol,ncol)
      real(wp)  :: binvp(ncol,ncol)
      real(wp)  :: ekin(ncol)
      real(wp)  :: xkv(nnot)
      real(wp)  :: xdv(ncol)
      real(wp)  :: tmp(nord,nord1)
      real(wp)  :: wsk1(ncol,ncol)
      real(wp)  :: ws1(ncol)
      real(wp)  :: wkr(ncol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ia, jj, ii, k
      real(wp) :: dx, sumt
!----------------------------------------------------------------------
!      this routine returns the following quantites for periodic
!      boundary conditions:
!
!      der1   :   first derivative operator in collocation space
!      der2   :   second derivative operator in collocation space
!      zkin   :   eigenvectors of the -0.5*der2
!      ekin   :   eigenvalues of -0.5*der2
!      wq     :   integration weights in collocation space
!      binvp  :   inverse b-matrix. this could be used to calculate
!                 expansion coefficients for interpolation
!      bwrk   :   array containing the non-zero splines at knots in
!                 the physical boundary
!
!      periodic boundary conditions are accomplished by closing the
!      mesh on both ends, i.e. the splines that are outside of our
!      physical left boundary are continued back into the right
!      boundary. so the first ncol splines are the same as always but
!      the remaining part of the mesh is covered by the extension of
!      the first nord-1 splines through the right boundary.
!----------------------------------------------------------------------
!    begin block to compute the b-matrix for inversion
!    use der2 as work space initially
!    call splgen to calculate all nonzero splines at all collocation
!    points.
!----------------------------------------------------------------------
      call splgen (0, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      do i = 1, ncol
         do j = 1, ncol
            der1(i,j) = 0.0_wp
            der2(i,j) = 0.0_wp
            zkin(i,j) = 0.0_wp
            cgwk(i,j) = 0.0_wp
            binvp(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            ii = i - jj + 1
            if (ii > ncol) ii = ii - ncol
            der2(ia,ii) = bwrk(jj,ia)
            der1(ia,ii) = der2(ia,ii)
         end do
      end do
!----------------------------------------------------------------------
!    invert the matrix and then test the inverse
!----------------------------------------------------------------------
      do i = 1, ncol
         binvp(i,i) = 1.0_wp
      end do
      call matq (der2, binvp, ncol, ncol, ncol, ncol)
      call invtst (der1, binvp, ncol)
!----------------------------------------------------------------------
!    calculate weights using the inverse stored in zkin. note that
!    we have an analytic expression for the integral of the splines.
!    this expression is valid for non-equal mesh spacings and also
!    for multiple knots.
!----------------------------------------------------------------------
      do i = 1, ncol
         wq(i) = 0.0_wp
         do k = 1, ncol
            dx = (xkv(k+nord)-xkv(k))/nord
            wq(i) = wq(i) + binvp(k,i)*dx
         end do
      end do
!----------------------------------------------------------------------
!    begin the calculation of the first derivative operator. the same
!    method of imposing periodic boundary conditions is used here.
!----------------------------------------------------------------------
      call splgen (1, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      do i = 1, ncol
         do j = 1, ncol
            der2(i,j) = 0.0_wp
            der1(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            ii = i - jj + 1
            if (ii > ncol) ii = ii - ncol
            der2(ia,ii) = bwrk(jj,ia)
         end do
      end do
!
      do i = 1, ncol
         do j = 1, ncol
            der1(i,j) = 0.0_wp
            do k = 1, ncol
               der1(i,j) = der1(i,j) + der2(i,k)*binvp(k,j)
            end do
         end do
      end do
!----------------------------------------------------------------------
!    begin the calculation of the second derivative operator. the
!    same method of imposing periodic boundary conditions is used
!
!    note: the second derivative can be calculated two ways: 1) a
!          direct calculation from second derivative splines.
!          2) as d1*d1.
!----------------------------------------------------------------------
      call splgen (2, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      do i = 1, ncol
         do j = 1, ncol
            cgwk(i,j) = 0.0_wp
         end do
      end do
!
      do ia = 1, ncol
         i = indx0(ia)
         do jj = 1, nord
            ii = i - jj + 1
            if (ii > ncol) ii = ii - ncol
            cgwk(ia,ii) = bwrk(jj,ia)
         end do
      end do
!
      do i = 1, ncol
         do j = 1, ncol
            der2(i,j) = 0.0_wp
            do k = 1, ncol
!       der2(i,j)=der2(i,j)+der1(i,k)*der1(k,j)
               der2(i,j) = der2(i,j) + cgwk(i,k)*binvp(k,j)
            end do
         end do
      end do
!----------------------------------------------------------------------
!    calculate the eigenvalues and eigenvectors of the kinetic
!    energy operator
!----------------------------------------------------------------------
      do i = 1, ncol
         ekin(i) = 0.0_wp
         ws1(i) = 0.0_wp
         iwk(i) = 0
         do j = 1, ncol
            zkin(i,j) = 0.0_wp
            cgwk(i,j) = 0.0_wp
            wsk1(i,j) = -0.5_wp*der2(i,j)
            wsk2(i,j) = 0.0_wp
         end do
      end do
!
      call jacobi (wsk1, ncol, ncol, ekin, zkin)
      do i = 1, ncol
         wkr(i) = ekin(i)
      end do
!
      do j = 1, ncol
         sumt = 0.0_wp
         do i = 1, ncol
            sumt = sumt + zkin(i,j)**2
         end do
         do i = 1, ncol
            zkin(i,j) = zkin(i,j)/sqrt(sumt)
            zkin(i,j) = zkin(i,j)/sqrt(wq(i))
         end do
      end do
      call ordeig (ncol, ncol, wkr, ekin, ws1, zkin, cgwk)
!----------------------------------------------------------------------
!    refill the spline array bwrk for output
!----------------------------------------------------------------------
      call splgen (0, xkv, xdv, bwrk, tmp, indx0, nord, nord1, nnot, ncol, nspl)
!
      return
      end subroutine setspp
      subroutine splgens(ip, xkv, xdv, bx, tmp, indx0, nord, nord1, nnot, nspl)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ip
      integer  :: nord
      integer , intent(in) :: nord1
      integer  :: nnot
      integer  :: nspl
      integer , intent(out) :: indx0
      real(wp)  :: xkv(nnot)
      real(wp) , intent(in) :: xdv
      real(wp) , intent(out) :: bx(nord)
      real(wp) , intent(inout) :: tmp(nord,nord1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, i, iz, jz, no, nom1, jp, ij, ip1, ip2
      real(wp) :: tcol, xnum1, di, gi, coef
!----------------------------------------------------------------------
!    subroutine to generate the collocation matrix b or derivatives
!    of b and return it in bx. the splines
!    are generated at the collocation points xdv using the knots
!    xkv. for ip=0 it calculates the matrix b, for ip=number it
!    calculates the derivative of order number. the calculation of
!    the derivative matrix uses essentially the same kind of
!    recurrance relation used for b (see notes)
!-----------------------------------------------------------------------
      do ii = 1, nord
            bx(ii) = 0.0_wp
         end do
!
      if (ip >= nord) return
      if (ip == 0) then
!-----------------------------------------------------------------------
!    branch for the calculation of the collocation matrix b. below
!    the detailed explanation of the coding is given (see also notes)
!
!    loop 40 is over the collocation points tcol=xdv .
!    i is the lower limit of the knot interval containing
!    the current collocation point, tcol.
!    zero a work array tmp, and set the initial condition for the
!    recursion.
!    loop 30 starts the iteration upto the maximum order required
!    loop 20 uses the triangular relation to calculate only the
!    non-zero splines for a given order (no). what we are doing is
!    to calculate the triangular form column by column until we reach
!    nord. the results are than stored in bx
!-----------------------------------------------------------------------
            tcol = xdv
            call findi (xkv, nnot, nord, nspl, tcol, i)
            indx0 = i
            do iz = 1, nord
               do jz = 1, nord1
                  tmp(iz,jz) = 0.0_wp
               end do
            end do
            tmp(1,2) = 1.0_wp
            do no = 2, nord
               nom1 = no - 1
               do ii = 1, no
                  jp = i - ii + 1
                  xnum1 = xkv(nom1+jp) - xkv(jp)
                  if (abs(xnum1) < 1.0E-15_wp) then
                     di = 0.0_wp
                  else
                     di = (tcol - xkv(jp))/(xkv(nom1+jp)-xkv(jp))
                  endif
                  xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
                  if (abs(xnum1) < 1.0E-15_wp) then
                     gi = 0.0_wp
                  else
                     gi = (xkv(nom1+jp+1)-tcol)/(xkv(nom1+jp+1)-xkv(jp+1))
                  endif
                  tmp(no,ii+1) = di*tmp(nom1,ii+1) + gi*tmp(nom1,ii)
               end do
            end do
            do ij = 1, nord
               bx(ij) = tmp(nord,ij+1)
            end do
         return
      endif
!----------------------------------------------------------------------
!    branch for the calculation of the derivative matrix.
!    we have essentially a duplicate of the above coding except
!    here the initial condition is also a recursive expression.
!    the details are the same as above and they are explained in notes
!----------------------------------------------------------------------
      ip1 = ip + 1
      ip2 = ip + 2
         tcol = xdv
         call findi (xkv, nnot, nord, nspl, tcol, i)
         indx0 = i
!----------------------------------------------------------------------
!    calculate the initial condition for a given colloc. point
!----------------------------------------------------------------------
         do iz = 1, nord
            do jz = 1, nord1
               tmp(iz,jz) = 0.0_wp
            end do
         end do
         tmp(1,2) = 1.0_wp
         do no = 2, ip1
            nom1 = no - 1
            do ii = 1, no
               jp = i - ii + 1
               xnum1 = xkv(nom1+jp) - xkv(jp)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  di = 0.0_wp
               else
                  di = nom1/(xkv(nom1+jp)-xkv(jp))
               endif
               xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  gi = 0.0_wp
               else
                  gi = nom1/(xkv(nom1+jp+1)-xkv(jp+1))
               endif
               tmp(no,ii+1) = di*tmp(nom1,ii+1) - gi*tmp(nom1,ii)
            end do
         end do
!----------------------------------------------------------------------
!    use the initial condition which is in tmp(ip1,i) i=1..nord to
!    calculate the recursion relation for the derivative. the result
!    is than stored in bx.
!----------------------------------------------------------------------
         do no = ip2, nord
            nom1 = no - 1
            do ii = 1, no
               jp = i - ii + 1
               xnum1 = xkv(nom1+jp) - xkv(jp)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  di = 0.0_wp
               else
                  di = (tcol - xkv(jp))/(xkv(nom1+jp)-xkv(jp))
               endif
               xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  gi = 0.0_wp
               else
                  gi = (xkv(nom1+jp+1)-tcol)/(xkv(nom1+jp+1)-xkv(jp+1))
               endif
               coef = nom1
               coef = coef/(nom1 - ip)
               tmp(no,ii+1) = (di*tmp(nom1,ii+1)+gi*tmp(nom1,ii))*coef
            end do
         end do
         do ij = 1, nord
            bx(ij) = tmp(nord,ij+1)
         end do
!
      return
      end subroutine splgens
      subroutine splgen(ip, xkv, xdv, bx, tmp, indx0, nord, nord1, nnot, ncol, nspl)
      implicit none
      integer, parameter :: wp = kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      integer  :: nord
      integer  :: nnot
      integer  :: nspl
      integer , intent(in) :: ip
      integer , intent(in) :: nord1
      integer , intent(in) :: ncol
      integer , intent(out) :: indx0(ncol)
      real(wp) , intent(in) :: xdv(ncol)
      real(wp) , intent(out) :: bx(nord,ncol)
      real(wp) , intent(inout) :: tmp(nord,nord1)
      real(wp)  :: xkv(nnot)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, jj, n, i, iz, jz, no, nom1, jp, ij, ip1, ip2
      real(wp) :: tcol, xnum1, di, gi, coef
!----------------------------------------------------------------------
!    subroutine to generate the collocation matrix b or derivatives
!    of b and return it in bx. the splines
!    are generated at the collocation points xdv using the knots
!    xkv. for ip=0 it calculates the matrix b, for ip=number it
!    calculates the derivative of order number. the calculation of
!    the derivative matrix uses essentially the same kind of
!    recurrance relation used for b (see notes)
!----------------------------------------------------------------------
      do ii = 1, nord
         do jj = 1, ncol
            bx(ii,jj) = 0.0_wp
         end do
      end do
!
      if (ip >= nord) return
      if (ip == 0) then
!----------------------------------------------------------------------
!    branch for the calculation of the collocation matrix b. below
!    the detailed explanation of the coding is given (see also notes)
!
!    loop 40 is over the collocation points tcol=xdv .
!    i is the lower limit of the knot interval containing
!    the current collocation point, tcol.
!    zero a work array tmp, and set the initial condition for the
!    recursion.
!    loop 30 starts the iteration upto the maximum order required
!    loop 20 uses the triangular relation to calculate only the
!    non-zero splines for a given order (no). what we are doing is
!    to calculate the triangular form column by column until we reach
!    nord. the results are than stored in bx
!----------------------------------------------------------------------
         do n = 1, ncol
            tcol = xdv(n)
            call findi (xkv, nnot, nord, nspl, tcol, i)
            indx0(n) = i
            do iz = 1, nord
               do jz = 1, nord1
                  tmp(iz,jz) = 0.0_wp
               end do
            end do
            tmp(1,2) = 1.0_wp
            do no = 2, nord
               nom1 = no - 1
               do ii = 1, no
                  jp = i - ii + 1
                  xnum1 = xkv(nom1+jp) - xkv(jp)
                  if(abs(xnum1) < 1.0E-15_wp) then
                     di = 0.0_wp
                  else
                     di = (tcol - xkv(jp))/(xkv(nom1+jp)-xkv(jp))
                  endif
                  xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
                  if(abs(xnum1) < 1.0E-15_wp) then
                     gi = 0.0_wp
                  else
                     gi = (xkv(nom1+jp+1)-tcol)/(xkv(nom1+jp+1)-xkv(jp+1))
                  endif
                  tmp(no,ii+1) = di*tmp(nom1,ii+1) + gi*tmp(nom1,ii)
               end do
            end do
            do ij = 1, nord
               bx(ij,n) = tmp(nord,ij+1)
            end do
         end do
         return
      endif
!----------------------------------------------------------------------
!    branch for the calculation of the derivative matrix.
!    we have essentially a duplicate of the above coding except
!    here the initial condition is also a recursive expression.
!    the details are the same as above and they are explained in notes
!----------------------------------------------------------------------
      ip1 = ip + 1
      ip2 = ip + 2
      do n = 1, ncol
         tcol = xdv(n)
         call findi (xkv, nnot, nord, nspl, tcol, i)
         indx0(n) = i
!----------------------------------------------------------------------
!    calculate the initial condition for a given colloc. point
!----------------------------------------------------------------------
         do iz = 1, nord
            do jz = 1, nord1
               tmp(iz,jz) = 0.0_wp
            end do
         end do
         tmp(1,2) = 1.0_wp
         do no = 2, ip1
            nom1 = no - 1
            do ii = 1, no
               jp = i - ii + 1
               xnum1 = xkv(nom1+jp) - xkv(jp)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  di = 0.0_wp
               else
                  di = nom1/(xkv(nom1+jp)-xkv(jp))
               endif
               xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  gi = 0.0_wp
               else
                  gi = nom1/(xkv(nom1+jp+1)-xkv(jp+1))
               endif
               tmp(no,ii+1) = di*tmp(nom1,ii+1) - gi*tmp(nom1,ii)
            end do
         end do
!----------------------------------------------------------------------
!    use the initial condition which is in tmp(ip1,i) i=1..nord to
!    calculate the recursion relation for the derivative. the result
!    is than stored in bx.
!----------------------------------------------------------------------
         do no = ip2, nord
            nom1 = no - 1
            do ii = 1, no
               jp = i - ii + 1
               xnum1 = xkv(nom1+jp) - xkv(jp)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  di = 0.0_wp
               else
                  di = (tcol - xkv(jp))/(xkv(nom1+jp)-xkv(jp))
               endif
               xnum1 = xkv(nom1+jp+1) - xkv(jp+1)
               if(abs(xnum1) <= 1.0E-15_wp) then
                  gi = 0.0_wp
               else
                  gi = (xkv(nom1+jp+1)-tcol)/(xkv(nom1+jp+1)-xkv(jp+1))
               endif
               coef = nom1
               coef = coef/(nom1 - ip)
               tmp(no,ii+1) = (di*tmp(nom1,ii+1)+gi*tmp(nom1,ii))*coef
            end do
         end do
         do ij = 1, nord
            bx(ij,n) = tmp(nord,ij+1)
         end do
      end do
!
      return
      end subroutine splgen
