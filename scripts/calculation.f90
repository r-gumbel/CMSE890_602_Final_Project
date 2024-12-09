
    program calculation
    implict none

    ! Placeholder calculation
    real :: x, result
    x = acos(-1.0)
    result = x**2

    ! Write result to output file 
    open(unit=12, file = 'outputs/fortran_results.txt', status='replace')
    write(12, *) 'Calculation result:', result
    close(12)

    end program calculation   
        