Subroutine moments(xclx, xcly, xclz, ncolx, ncoly, ncolz, wx, wy, wz, xcm, ycm, zcm, rho, q20, q22, q20tot, q22tot, dxq, dzq, &
& idp, etheta)
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
      Real(wp), Intent(In) :: xcm(2)
      Real(wp), Intent(In) :: ycm(2)
      Real(wp), Intent(In) :: zcm(2)
      Real(wp), Intent(Out) :: q20(2)
      Real(wp), Intent(Out) :: q22(2)
      Real(wp), Intent(Out) :: q20tot
      Real(wp), Intent(Out) :: q22tot
      Real(wp), Intent(Out) :: dxq, dzq
      Real(wp), Intent(Out) :: etheta(3)
      Integer, Intent(Out) :: idp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp) :: q2(3, 3, 2), q2mat(3, 3), xx, yy, zz, vol, dyq
!      Real(wp) :: q2(3, 3, 2), discr, slope
      Integer :: iq, ix, iy, iz
!-------------------------------------------------------
!   calculate the traceless-symmetric quadrupole tensor
!-------------------------------------------------------
      q2 = 0.0_wp
      Do iq = 1, 2
         Do iz = 1, ncolz
            zz = xclz(iz) - zcm(iq)
            Do iy = 1, ncoly
               yy = xcly(iy) - ycm(iq)
               Do ix = 1, ncolx
                  xx = xclx(ix) - xcm(iq)
                  vol = wx(ix) * wy(iy) * wz(iz) * rho(ix, iy, iz, iq)
                  q2(1, 1, iq) = q2(1, 1, iq) + (xx*xx+xx*xx-yy*yy-zz*zz) * vol
                  q2(2, 2, iq) = q2(2, 2, iq) + (yy*yy+yy*yy-xx*xx-zz*zz) * vol
                  q2(3, 3, iq) = q2(3, 3, iq) + (zz*zz+zz*zz-xx*xx-yy*yy) * vol
                  q2(1, 2, iq) = q2(1, 2, iq) + 3.0_wp * xx * yy * vol
                  q2(1, 3, iq) = q2(1, 3, iq) + 3.0_wp * xx * zz * vol
                  q2(2, 3, iq) = q2(2, 3, iq) + 3.0_wp * yy * zz * vol
               End Do
            End Do
         End Do
         q2(2, 1, iq) = q2(1, 2, iq)
         q2(3, 1, iq) = q2(1, 3, iq)
         q2(3, 2, iq) = q2(2, 3, iq)
      End Do
!
      ! q2mat = q2(:, :, 1) + q2(:, :, 2)
      ! discr = 0.5_wp * (q2mat(2, 2)-q2mat(3, 3)+Sqrt(Abs(q2mat(2, 2)**2-2.0_wp*q2mat(2, 2)*q2mat(3, 3)+q2mat(3, 3)**2+4.0_wp*q2mat(3, 2)**2)))
      ! If(Abs(discr) < 1.0e-4_wp) Then
      !    slope = 100.0_wp
      !    etheta(1) = 0.0_wp
      ! Else
      !    slope = q2mat(2, 3) / discr
      !    etheta(1) = Atan(slope)
      ! End If
      ! discr = 0.5_wp * (q2mat(1, 1)-q2mat(3, 3)+Sqrt(Abs(q2mat(1, 1)**2-2.0_wp*q2mat(1, 1)*q2mat(3, 3)+q2mat(3, 3)**2+4.0_wp*q2mat(3, 1)**2)))
      ! If(Abs(discr) < 1.0e-4_wp) Then
      !    slope = 100.0_wp
      !    etheta(2) = 0.0_wp
      ! Else
      !    slope = q2mat(1, 3) / discr
      !    etheta(2) = Atan(slope)
      ! End If
      ! discr = 0.5_wp * (q2mat(1, 1)-q2mat(2, 2)+Sqrt(Abs(q2mat(1, 1)**2-2.0_wp*q2mat(1, 1)*q2mat(2, 2)+q2mat(2, 2)**2+4.0_wp*q2mat(2, 1)**2)))
      ! If(Abs(discr) < 1.0e-4_wp) Then
      !    slope = 100.0_wp
      !    etheta(3) = 0.0_wp
      ! Else
      !    slope = q2mat(1, 2) / discr
      !    etheta(3) = Atan(slope)
      ! End If
      q2mat = q2(:, :, 1)
      Call q2diag(q2mat, q20(1), q22(1), dxq, dyq, dzq, idp)
      q2mat = q2(:, :, 2)
      Call q2diag(q2mat, q20(2), q22(2), dxq, dyq, dzq, idp)
! total should be called last because of id
      q2mat = q2(:, :, 1) + q2(:, :, 2)
      Call q2diag(q2mat, q20tot, q22tot, dxq, dyq, dzq, idp)
      etheta(1) = Atan(dyq/dxq)
      etheta(2) = Atan(dzq/dxq)
!
End Subroutine moments
Subroutine q2diag(q2, q20x, q22x, dxq, dyq, dzq, idp)
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!   A r g u m e n t s
!-----------------------------------------------
      Real(wp), Intent(In) :: q2(3, 3)
      Real(wp), Intent(Out) :: q20x, q22x
      Real(wp), Intent(Out) :: dxq, dyq, dzq
      Integer, Intent(Out) :: idp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Real(wp), Parameter :: pi = acos(-1.0_wp)
      Real(wp) :: w(3), z(3, 3)
      Integer :: loc1(1)
!      integer              :: i
!
      Call jacobi(q2, 3, 3, w, z)
!
! now calculate spherical moments
      loc1 = maxloc(w)
      idp = loc1(1)
      q20x = Sqrt(5.0_wp/(16.0_wp*pi)) * w(idp)
      dxq = z(1, idp)
      dyq = z(2, idp)
      dzq = z(3, idp)
!      WRITE(8,*) 'Cartesian quadrupole principal values and axes:',idp
!      WRITE(8,'(3(F10.2,''('',3F8.4,'')''))') (w(i),z(:,i),i=1,3)
! Re(q22) = SQRT(5.0_wp/(96.0_wp*pi))*(Q11 - Q33) (tensor components)
      If(idp == 1) Then
         q22x = Sqrt(5.0_wp/(96.0_wp*pi)) * (w(2)-w(3))
      Else If(idp == 2) Then
         q22x = Sqrt(5.0_wp/(96.0_wp*pi)) * (w(3)-w(1))
      Else If(idp == 3) Then
         q22x = Sqrt(5.0_wp/(96.0_wp*pi)) * (w(1)-w(2))
      End If
End Subroutine q2diag
