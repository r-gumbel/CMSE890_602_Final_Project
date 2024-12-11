Subroutine eorder(spe, iorder, n)
!----------------------------------------------------
! create an index array containing ordered energies
!----------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-----------------------------------------------
!      A r g u m e n t s
!-----------------------------------------------
      Integer, Intent(In) :: n
      Integer, Intent(Inout) :: iorder(n)
      Real(wp), Intent(Inout) :: spe(n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      Integer :: i, k, loc(1)
      Logical, Dimension(n) :: selected
!
      selected = .False.
      k = 0
      Do i = 1, n
         k = k + 1
         loc = minloc(spe, mask= .Not. selected)
         selected(loc(1)) = .True.
         iorder(k) = loc(1)
      End Do
!
      Return
End Subroutine eorder
