Subroutine getslope(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcmtot, ycmtot, zcmtot, rho, slope)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: xclx(ncolx)
      Real(wp), Intent(In) :: xcly(ncoly)
      Real(wp), Intent(In) :: xclz(ncolz)
      Real(wp), Intent(In) :: wx(ncolx)
      Real(wp), Intent(In) :: wy(ncoly)
      Real(wp), Intent(In) :: wz(ncolz)
      Real(wp), Intent(In) :: rho(ncolx, ncoly, ncolz, 2)
      Real(wp), Intent(In) :: xcmtot
      Real(wp), Intent(In) :: ycmtot
      Real(wp), Intent(In) :: zcmtot
      Real(wp), Intent(Out) :: slope
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: q2(3, 3), xx, yy, zz, vol, discr
      Integer :: ix, iy, iz
!-------------------------------------------------------
!   calculate the traceless-symmetric quadrupole tensor
!-------------------------------------------------------
      q2 = 0.0_wp
      Do iz = 1, ncolz
         zz = xclz(iz) - zcmtot
         Do iy = 1, ncoly
            yy = xcly(iy) - ycmtot
            Do ix = 1, ncolx
               xx = xclx(ix) - xcmtot
               vol = wx(ix) * wy(iy) * wz(iz) * (rho(ix, iy, iz, 1)+rho(ix, iy, iz, 2))
               q2(1, 1) = q2(1, 1) + (xx*xx+xx*xx-yy*yy-zz*zz) * vol
               q2(2, 2) = q2(2, 2) + (yy*yy+yy*yy-xx*xx-zz*zz) * vol
               q2(3, 3) = q2(3, 3) + (zz*zz+zz*zz-xx*xx-yy*yy) * vol
               q2(1, 2) = q2(1, 2) + 3.0_wp * xx * yy * vol
               q2(1, 3) = q2(1, 3) + 3.0_wp * xx * zz * vol
               q2(2, 3) = q2(2, 3) + 3.0_wp * yy * zz * vol
            End Do
         End Do
      End Do
      q2(2, 1) = q2(1, 2)
      q2(3, 1) = q2(1, 3)
      q2(3, 2) = q2(2, 3)
!-----------------------------------------------------------------------------------
!   if we assume that the only rotation takes place in the x-z plane, we can
!   set all elements of Q_ij with i=2 or j=2 to zero. In this case one can
!   analytically calculate the eigenvalues and eigenvectors as:
!   Eigenvalues:
!   0
!   1/2 (q11+q33-Sqrt[q11^2+4 q13^2-2 q11 q33+q33^2])
!   1/2 (q11+q33+Sqrt[q11^2+4 q13^2-2 q11 q33+q33^2]) <----largest eigenvalue
!   eigen_l=0.5_wp*(q2(1,1) + q2(3,3) + sqrt(q2(1,1)**2 - 2.0_wp*q2(1,1)*q2(3,3)&
!                           + q2(3,3)**2 + 4.0_wp*q2(3,1)**2))
!
!   Eigenvectors:
!   {0,1,0}
!   {-((-q11+q33+Sqrt[q11^2+4 q13^2-2 q11 q33+q33^2])/(2 q13)),0,1}
!   {-((-q11+q33-Sqrt[q11^2+4 q13^2-2 q11 q33+q33^2])/(2 q13)),0,1} <--- eigenvector
!   so, slope=(q11-q33+Sqrt[q11^2+4 q13^2-2 q11 q33+q33^2])/(2 q13)
!   which agrees with dz/dx coming from exact diagonalization.
!-----------------------------------------------------------------------------------
      discr = 0.5_wp * (q2(1, 1)-q2(3, 3)+Sqrt(Abs(q2(1, 1)**2-2.0_wp*q2(1, 1)*q2(3, 3)+q2(3, 3)**2+4.0_wp*q2(3, 1)**2)))
      If(Abs(discr) < 1.0e-4_wp) Then
         slope = 100.0_wp
      Else
         slope = q2(1, 3) / discr
      End If
End Subroutine getslope
