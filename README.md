![ForDiff](media/logo.png)
============

**ForDiff**: Numerical differentiation

-----


## Table of Contents

- [](#)
  - [Table of Contents](#table-of-contents)
  - [TO DO](#to-do)
  - [Installation](#installation)
    - [fpm](#fpm)
  - [Usage](#usage)
  - [Tests](#tests)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
-----

## TO DO
- [x] Complex-step: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Complex-step: f(x) f is a scalar-valued function and x is a vector variable
- [x] Complex-step: f(x) f is a vector-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a vector-valued function and x is a vector variable
- [ ] Automatic Differentiation
-----
## Installation

### fpm
fordiff can be cloned and then built using [fpm](https://github.com/fortran-lang/fpm), following the instructions provided in the documentation available on Fortran Package Manager.

```bash
git clone https://github.com/gha3mi/fordiff.git
cd fordiff
fpm install --perfix .
```

Or you can easily include this package as a dependency in your `fpm.toml` file.

```toml
[dependencies]
fordiff = {git="https://github.com/gha3mi/fordiff.git"}
```

-----

## Usage
Here is an example of how to use the `fordiff` module in your Fortran code:
```fortran
module mod_func1

   use kinds
   implicit none

contains

   function func1(x) result(f)
      complex(rk), intent(in)  :: x
      complex(rk)              :: f

      f = x**2 + 2.0_rk*x

   end function func1

end module mod_func1

program test1

   use kinds
   use mod_func1
   use fordiff

   implicit none

   real(rk) :: dfdx

   dfdx = derivative(f=func1, x=1.0_rk, h=1e-100_rk)

   print*,dfdx

end program test1
```
-----

## Tests

The `tests` directory contains test programs to verify the functionality of the `fordiff` module. To run the tests using `fpm`, you can use response files for specific compilers:

- For Intel Fortran Compiler (ifort):
```bash
fpm @ifort
```

- For Intel Fortran Compiler (ifx):
```bash
fpm @ifx
```

- For NVIDIA Compiler (nvfortran):
```bash
fpm @nvidia
```

- For GNU Fortran Compiler (gfortran):
```bash
fpm @gfortran
```

-----

## Documentation
To generate the documentation for the `fordiff` module using [ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following command:
```bash
ford ford.yml
```

-----

## Contributing

Contributions to fordiff are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.