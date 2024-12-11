Subroutine createdir(directname)
      Implicit None
      Character(Len=*) directname
!
      Call system('mkdir -p '//TRIM(directname))
End Subroutine createdir
