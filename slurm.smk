#Snakefile
import os

# Configuration
EXECUTABLE_NAME = 'hello_world'

# Rule to generate some test fortran source code
rule create_fortran_source:
    output:
        "src/main.f90"
    run:
        # Ensure the src directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # Write the fortran source code 
        with open(output[0], 'w') as f:
            f.write('''program hello_world
    print *, "Hello World!"
end program hello_world''')
