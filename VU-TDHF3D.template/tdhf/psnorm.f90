Subroutine psnorm(ncolx, ncoly, ncolz, pl, pr, anorm, irn, wxyz)
!-----------------------------------------------------------------------
!     pnorm - pgm to compute the complex norm of two vectors
!     if irn is 0 the vectors are renormalized and assumes for this
!     case the vectors are self adjoint.
!     otherwise the norm is returned
!-----------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: ncolx
      Integer, Intent(In) :: ncoly
      Integer, Intent(In) :: ncolz
      Integer, Intent(In) :: irn
      Real(wp), Intent(In) :: wxyz(ncolx, ncoly, ncolz)
      Complex(wp), Intent(Out) :: anorm
      Complex(wp) :: pl(ncolx, ncoly, ncolz, 2)
      Complex(wp) :: pr(ncolx, ncoly, ncolz, 2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: is
!-----------------------------------------------
!
      anorm = (0.0_wp, 0.0_wp)
      Do is = 1, 2
         anorm = anorm + Sum(wxyz(:, :, :)*conjg(pl(:, :, :, is))*pr(:, :, :, is))
      End Do
!
      If(irn /= 0) Return
!
      pr = pr / Sqrt(anorm)
!
      Return
End Subroutine psnorm
