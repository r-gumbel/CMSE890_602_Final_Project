Subroutine distance(ncolx,ncoly,ncolz,xclx,xcly,xclz,wx,wy,wz,rho,damp)
!.......................................................................
!     cutoff of the constraint depending on the density distribution and
!     profile as proposed by K. Rutz et al, Nucl. Phys. A590 (1995) 690
!
!     adopted from EV8 version II code
!
!     this subroutine calculates the minimum distance of each point on the
!     mesh to the surface which is defined by having the value sval in the
!     mesh function func and then defines the cutoff function damp from
!     this distance such that
!
!                       1
!          --------------------------
!           1 + exp[(dist-rcut)/acut]
!
!       sval    : switch value, recommended value is 1/10 of the density max
!       rcut    : additional distance from the surface
!       acut    : width of the cutoff
!       damp    : cutoff function
!
!     the values rcut = 4.0 fm and acut = 0.4 fm should be used.
!.......................................................................
      Implicit None
      Integer, Parameter  :: wp = Kind(1.0D0)
      Real(wp), Parameter :: rcut = 4.0_wp
      Real(wp), Parameter :: acut = 0.4_wp
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In)   :: ncolx
      Integer, Intent(In)   :: ncoly
      Integer, Intent(In)   :: ncolz
      Real(wp), Intent(In)  :: xclx(ncolx)
      Real(wp), Intent(In)  :: xcly(ncoly)
      Real(wp), Intent(In)  :: xclz(ncolz)
      Real(wp), Intent(In)  :: wx(ncolx)
      Real(wp), Intent(In)  :: wy(ncoly)
      Real(wp), Intent(In)  :: wz(ncolz)
      Real(wp), Intent(In)  :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(Out) :: damp(ncolx, ncoly, ncolz)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer  :: ix, iy, iz, ixx, iyy, izz, jx, is, ismax, mq
      Real(wp) :: func(ncolx, ncoly, ncolz), dist(ncolx, ncoly, ncolz)
      Real(wp) :: xi(ncolx, ncoly, ncolz), yi(ncolx, ncoly, ncolz), zi(ncolx, ncoly, ncolz)
      Real(wp) :: surf(ncolx*ncoly*ncolz,3), sig(ncolx, ncoly, ncolz)
      Real(wp) :: x, y, z, xx, yy, zz
      Real(wp) :: sval, fxyz, ffxyz, frac, earg, surfx, surfy, surfz, d, dmin
!
!........................................................... define the surface
!                         i[x,y,z]    are the coordinates of the current  point
!                         i[xx,yy,zz] are the coordinates of the adjacent point

      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               xi(ix, iy, iz) = xclx(ix)
               yi(ix, iy, iz) = xcly(iy)
               zi(ix, iy, iz) = xclz(iz)
               func(ix, iy, iz) = rho(ix, iy, iz, 1) + rho(ix, iy, iz, 2)
            End Do
         End Do
      End Do
      sval = Maxval(func) / 10.0_wp
!
      is = 0         ! index of field with surface points
      mq = ncolx*ncoly*ncolz
      surf = 0.0_wp
                     ! initialization - to be checked
      surfx = 75.0d0 ! sometimes, the values are not set below, which is not
      surfy = 75.0d0 ! a problem for the algorithm when they remain zero,
      surfz = 75.0d0 ! but might have some consequences on some platforms
!
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               x = xi(ix, iy, iz)
               y = yi(ix, iy, iz)
               z = zi(ix, iy, iz)
               fxyz = func(ix, iy, iz)
               If(fxyz .Lt. sval) Then
                  sig(ix, iy, iz) = 1.0_wp
               Else
                  sig(ix, iy, iz) = -1.0_wp
               End If
               ixx = ix
               iyy = iy
               izz = iz
!           ......................................................  x direction
               If(ix /= 1) ixx = ix - 1
               ffxyz = func(ixx, iyy, izz)
               If((fxyz < sval .And. ffxyz < sval) .Or. (fxyz >= sval .And. ffxyz >= sval)) Then
                  jx = ix ! both points inside surface or outside surface, do nothing!
               Else ! one point inside one point outside - surface in between
                  frac = (sval-fxyz) / (ffxyz-fxyz)
                  If(ix /= ixx) surfx = x - frac * wx(ix)
                  If(is >= mq) Then
                     STOP 'distance: is > mq'
                  Else
                     is = is + 1
                     surf(is, 1) = surfx
                     surf(is, 2) = y
                     surf(is, 3) = z
                  End If
               End If
!
               ixx = ix
!           ......................................................  y direction
               If(iy /= 1) iyy = iy - 1
               ffxyz = func(ixx, iyy, izz)
               If((fxyz < sval .And. ffxyz < sval) .Or. (fxyz >= sval .And. ffxyz >= sval)) Then
                  jx = ix
               Else
                  frac = (sval-fxyz) / (ffxyz-fxyz)
                  If(iy /= iyy) surfy = y - frac * wy(iy)
                  If(is >= mq) Then
                     STOP 'distance: is > mq'
                  Else
                     is = is + 1
                     surf(is, 1) = x
                     surf(is, 2) = surfy
                     surf(is, 3) = z
                  End If
               End If
!
               iyy = iy
!           ......................................................  z direction
               If(iz /= 1) izz = iz - 1
               ffxyz = func(ixx, iyy, izz)
               If((fxyz < sval .And. ffxyz < sval) .Or. (fxyz >= sval .And. ffxyz >= sval)) Then
                  jx = ix
               Else
                  frac = (sval-fxyz) / (ffxyz-fxyz)
                  If(iz /= izz) surfz = z - frac * wz(iz)
                  If(is >= mq) Then
                     STOP 'distance: is > mq'
                  Else
                     is = is + 1
                     surf(is, 1) = x
                     surf(is, 2) = y
                     surf(is, 3) = surfz
                  End If
               End If
               izz = iz
            End Do
         End Do
!
      End Do
!
      ismax = is
!     ...................................................... calculate distance
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               dmin = 1.e+12_wp
               x = xi(ix, iy, iz)
               y = yi(ix, iy, iz)
               z = zi(ix, iy, iz)
               Do is = 1, ismax
                  xx = surf(is, 1)
                  yy = surf(is, 2)
                  zz = surf(is, 3)
                  d = (x-xx) * (x-xx) + (y-yy) * (y-yy) + (z-zz) * (z-zz)
                  dmin = Min(dmin, d)
               End Do
               dist(ix, iy, iz) = Sqrt(dmin) * sig(ix, iy, iz)
            End Do
         End Do
!
      End Do
!     ................................ calculate cutoff function for constraint
      Do iz = 1, ncolz
         Do iy = 1, ncoly
            Do ix = 1, ncolx
               earg = (dist(ix, iy, iz)-rcut) / acut
               damp(ix, iy, iz) = 1.0_wp / (1.0_wp + Exp(earg))
            End Do
         End Do
!
!
      End Do
      Return
End Subroutine distance
