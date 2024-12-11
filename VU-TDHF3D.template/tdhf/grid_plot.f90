Subroutine grid_plot(d, header, xg, yg, nx, ny, unit, istyl)
!------------------------------------------------------------------------------
!   purpose
!   =======
!
!   grid_plot constructs an output file in the mtvdat format for 2-d
!   data stored in array d.
!
!   description
!   ===========
!
!   this routine provides optional control over the main title, selection of
!   the output file and the choice of plot styles.   control for the labeling,
!   extent, and location of plot points on each axis is also provided. the
!   maximum number of contours allowed (40) is always specified,
!   if fewer countours are desired the output file can be edited manually.
!
!   arguments
!   =========
!
!          d  - (required input) real(wp) 2-d array
!               data to be represented in the contour plot.  if no optional
!               arguments are specified the `graduate color' style will be
!               selected.
!     header  - character(len=*)
!               header is passed verbatim to plotmtv as the main title.
!      xg     - 1-d x array
!      yg     - 1-d y array
!      unit   - output unit
!      istyl  -
!                value            action
!                -----            ------
!                default   --->  graduated colors
!                  1       --->  normal contours
!                  2       --->  graduated colors
!                  3       --->  3-d surface mesh
!                 else     --->  graduated colors
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
      Integer, Intent(In) :: nx, ny, unit, istyl
      Real(wp), Intent(In) :: d(nx, ny)
      Real(wp), Intent(In) :: xg(nx), yg(ny)
      Character(Len=*) :: header
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: j, lu
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      Integer, Parameter :: nsteps = 40
!
      lu = unit
!  specify plot style ..
      Write(lu,*) '$data=contour'
      Write(lu,*) '%contstyle=', istyl
!  specify the grid ..
      Write(lu,*) '%nx=', nx, ' xgrid=true'
      Write(lu,*) xg
      Write(lu,*) '%ny=', ny, ' ygrid=true'
      Write(lu,*) yg
!  specify axis minima/maxima
      Write(lu,*) '%xmin=', xg(1)
      Write(lu,*) '%ymin=', yg(1)
      Write(lu,*) '%xmax=', xg(nx)
      Write(lu,*) '%ymax=', yg(ny)
!  specify number of contours ..
      Write(lu,*) '%nsteps=', nsteps
!  specify labels
      Write(lu,*) '%toplabel=', header
      Write(lu,*) '%xlabel=', 'x'
      Write(lu,*) '%ylabel=', 'z'
      Write(lu,*) '%zlabel=', 'f(x,z)'
!  write countour data ..
      Do j = 1, ny
         Write(lu,*) d(:, j)
      End Do
!
      Return
End Subroutine grid_plot
