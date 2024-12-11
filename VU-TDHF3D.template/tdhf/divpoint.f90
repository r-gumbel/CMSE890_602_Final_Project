Subroutine divpoint(ncolx, ncoly, ncolz, xclx, xclz, slope, bb, ixc1, ixc2, izc1, izc2, rho, xmin, zmin, ctime)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0d0)
!---------------------------------------------------------------------------
!     subroutine to calculate the location of the mid-point between the
!     two nuclei in the x-z plane. collision axis is the x-axis and for
!     non-zero impact parameters the nuclei are placed on the z-axis.
!
!     the calculation of the dividing plane starts by examining the line
!     along the axis of largest quadrupole moment found in getslope. We
!     integrate the density along a line perpindicular to the slope
!     at each point and find the minimum integrated density.
!
!     we interpolate the density in x-z plane along the principal axis
!     i.e., equation z=slope*x+bb.
!---------------------------------------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: ixc1, ixc2, izc1, izc2
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: slope
      Real(wp), Intent(In) :: bb
      Real(wp), Intent(Inout) :: xmin
      Real(wp), Intent(Inout) :: zmin
      Real(wp), Intent(Inout) :: ctime
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer, Parameter :: np = 500
      Real(wp), Parameter :: small = 1.0e-25_wp
      Real(wp), Parameter :: contour = 0.08_wp
      Integer :: iz1, iz2, ix1, ix2, il, ix, iz, locxz(1), npz, npx, fragments, npmax
      Real(wp) :: xx, zz, dx, dz, bplina, slopev, bbv, angle, xminold, zminold, rhoold
      Real(wp) :: rhob(ncolx, ncolz), rhol(np*np), rhopl(0:np*np), rhosum(np*np), xa(np), za(np), xap(np), zap(np)
!
      Do ix = 1, ncolx
         Do iz = 1, ncolz
            !rhob(ix,iz) = sum(wy(:)*(rho(ix,:,iz,1) + rho(ix,:,iz,2)))
            rhob(ix, iz) = (rho(ix, ncoly/2, iz, 1)+rho(ix, ncoly/2, iz, 2))
         End Do
      End Do
!
      xminold = xmin
      zminold = zmin
      rhol = 1.e+30_wp
      rhopl = 0.0_wp
! indices determining nuclear centers on numerical mesh (from tinfo)
      ix1 = Min(ixc1, ixc2)
      ix2 = Max(ixc1, ixc2)
      iz1 = Min(izc1, izc2)
      iz2 = Max(izc1, izc2)
      npz = Min(20*Abs(izc2-izc1), np)
      npx = Min(20*Abs(ixc2-ixc1), np)
      If(ix1 == ix2) npx = 2
      If(iz1 == iz2) npz = 2
! generate interpolation arrays between nuclear centers
      dx = (xclx(ix2)-xclx(ix1)) / (npx-1)
      dz = (xclz(iz2)-xclz(iz1)) / (npz-1)
      Do ix = 1, npx
         xa(ix) = xclx(ix1) + (ix-1) * dx
      End Do
      Do iz = 1, npz
         za(iz) = xclz(iz1) + (iz-1) * dz
      End Do
! interpolate the density between the two centers (along the slope)
      slopev = - 1.0_wp / (slope+small)
      angle = Atan(slopev)
      If(Abs(slope) < 1.0_wp) Then
         npmax = npx
         il = 0
         rhosum = 1.e+30_wp
         Do ix = 1, npx
            xx = xa(ix)
            zz = slope * xx + bb
            bbv = zz - slopev * xx
            zap(1:2) = zz
            xap(1:2) = xx
            rhosum(ix) = 0.0_wp
            rhopl(0) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(1), zap(1))
            rhosum(ix) = rhosum(ix) + rhopl(0)
            Do iz = 1, npx
               zap(1) = zap(1) + iz * dx * Sin(angle)
               xap(1) = xap(1) + iz * dx * Cos(angle)
               rhopl(iz) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(1), zap(1))
               rhosum(ix) = rhosum(ix) + rhopl(iz)
               zap(2) = zap(2) - iz * dx * Sin(angle)
               xap(2) = xap(2) - iz * dx * Cos(angle)
               rhopl(iz) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(2), zap(2))
               rhosum(ix) = rhosum(ix) + rhopl(iz)
            End Do
            rhol(ix) = rhopl(0)
         End Do
         locxz = minloc(rhosum)
         ix = 1
         xmin = xa(locxz(1))
         zmin = slope * xmin + bb
      ! Need to add some fail cases here and around the code
         Do While((locxz( 1) .Eq. ix) .Or. (locxz(1) .Eq. npx))
            If(ix .Eq. npx) Then
               xmin = xminold
               zmin = zminold
               Exit
            Else If(locxz(1) .Eq. ix) Then
               rhosum(ix) = 1.e+30_wp
               ix = ix + 1
            Else If(locxz(1) .Eq. npx) Then
               rhosum(npx) = 1.e+30_wp
               npx = npx - 1
            End If
            locxz = minloc(rhosum)
            xmin = xa(locxz(1))
            zmin = slope * xmin + bb
         End Do
! Case where the slope is greater than 1.0
      Else
         npmax = npz
         rhosum = 1.e+30_wp
         Do iz = 1, npz
            zz = za(iz)
            xx = (zz-bb) / (slope+small)
            bbv = zz - slopev * xx
            zap(1:2) = zz
            xap(1:2) = xx
            rhosum(iz) = 0.0_wp
            rhopl(0) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(1), zap(1))
            rhosum(iz) = rhosum(iz) + rhopl(0)
            Do ix = 1, npz
               zap(1) = zap(1) + ix * dz * Sin(angle)
               xap(1) = xap(1) + ix * dz * Cos(angle)
               rhopl(ix) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(1), zap(1))
               rhosum(iz) = rhosum(iz) + rhopl(ix)
               zap(2) = zap(2) - ix * dz * Sin(angle)
               xap(2) = xap(2) - ix * dz * Cos(angle)
               rhopl(ix) = bplina(ncolx, ncolz, xclx, xclz, rhob, xap(2), zap(2))
               rhosum(iz) = rhosum(iz) + rhopl(ix)
            End Do
            rhol(iz) = rhopl(0)
         End Do
         locxz = minloc(rhosum)
         iz = 1
         zmin = za(locxz(1))
         xmin = (zmin-bb) / (slope+small)
      ! Need to add some fail cases here and around the code
         Do While((locxz( 1) .Eq. iz) .Or. (locxz(1) .Eq. npz))
            If(iz .Eq. npz) Then
               xmin = xminold
               zmin = zminold
               Exit
            Else If(locxz(1) .Eq. iz) Then
               rhosum(iz) = 1.e+30_wp
               iz = iz + 1
            Else If(locxz(1) .Eq. npz) Then
               rhosum(npz) = 1.e+30_wp
               npz = npz - 1
            End If
            locxz = minloc(rhosum)
            zmin = za(locxz(1))
            xmin = (zmin-bb) / (slope+small)
         End Do
      End If
      rhoold = 0.0_wp
      Do ix = 2,npmax-1
         If(fragments==0 .And. rhol(ix)>contour) Then
            fragments = 1
         Else If (fragments==1 .And. rhol(ix)>contour .And. rhoold<contour) Then
            fragments = 2
         End If
         rhoold = rhol(ix)
      End Do
      If(fragments == 1) Then
         ctime = ctime + 1
      End If
!      If(xb <= 1.0e-2_wp) zmin = 0.0_wp
!
End Subroutine divpoint
