Subroutine setdmc(ncolx, ncoly, ncolz, unitx, unity, unitz, der2x, der2y, der2z, cdmpx, cdmpy, cdmpz, h2m, e0dmp)
!-----------------------------------------------------------------------
!       sets up the operator for damping with the non-relativistic
!       kinetic-energy:
!                                    1
!                cdmp(...) = -------------------
!                              (1.0 + t/e0dmp)
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Real(wp), Intent(In) :: e0dmp
      Real(wp), Intent(In) :: h2m(2)
      Real(wp), Intent(In) :: der2x(ncolx, ncolx)
      Real(wp), Intent(In) :: der2y(ncoly, ncoly)
      Real(wp), Intent(In) :: der2z(ncolz, ncolz)
      Real(wp) :: unitx(ncolx, ncolx)
      Real(wp) :: unity(ncoly, ncoly)
      Real(wp) :: unitz(ncolz, ncolz)
      Real(wp) :: cdmpx(ncolx, ncolx)
      Real(wp) :: cdmpy(ncoly, ncoly)
      Real(wp) :: cdmpz(ncolz, ncolz)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: i, j
      Real(wp) :: h2ma
!-----------------------------------------------
      h2ma = (h2m(1)+h2m(2)) / 2.0_wp
      cdmpx = 0.0_wp
      Do i = 1, ncolx
         Do j = 1, ncolx
            unitx(j, i) = - h2ma * der2x(j, i) / e0dmp
         End Do
         unitx(i, i) = 1.0_wp + unitx(i, i)
         cdmpx(i, i) = 1.0_wp
      End Do
      Call matq(unitx, cdmpx, ncolx, ncolx, ncolx, ncolx)
!
      cdmpy = 0.0_wp
      Do i = 1, ncoly
         Do j = 1, ncoly
            unity(j, i) = - h2ma * der2y(j, i) / e0dmp
         End Do
         unity(i, i) = 1.0_wp + unity(i, i)
         cdmpy(i, i) = 1.0_wp
      End Do
      Call matq(unity, cdmpy, ncoly, ncoly, ncoly, ncoly)
!
      cdmpz = 0.0_wp
      Do i = 1, ncolz
         Do j = 1, ncolz
            unitz(j, i) = - h2ma * der2z(j, i) / e0dmp
         End Do
         unitz(i, i) = 1.0_wp + unitz(i, i)
         cdmpz(i, i) = 1.0_wp
      End Do
      Call matq(unitz, cdmpz, ncolz, ncolz, ncolz, ncolz)
!
      Return
End Subroutine setdmc
