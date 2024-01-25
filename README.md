[![GitHub](https://img.shields.io/badge/GitHub-fordiff-blue.svg?style=social&logo=github)](https://github.com/gha3mi/fordiff)
[![Version](https://img.shields.io/github/release/gha3mi/fordiff.svg)](https://github.com/gha3mi/fordiff/releases/latest)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/fordiff/)
[![License](https://img.shields.io/github/license/gha3mi/fordiff?color=green)](https://github.com/gha3mi/fordiff/blob/main/LICENSE)
[![Build](https://github.com/gha3mi/fordiff/actions/workflows/CI_test.yml/badge.svg)](https://github.com/gha3mi/fordiff/actions/workflows/CI_test.yml)

**ForDiff**: A Fortran library for numericall differentiation


## Table of Contents

- [Table of Contents](#table-of-contents)
- [fpm dependency](#fpm-dependency)
- [Usage](#usage)
- [Run Tests](#run-tests)
- [TO DO](#to-do)
- [API documentation](#api-documentation)
- [Contributing](#contributing)

## fpm dependency

To use `ForDiff` as a dependency in your fpm project, include the following line in your `fpm.toml` file:

```toml
[dependencies]
fordiff = {git="https://github.com/gha3mi/fordiff.git"}
```

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

## Run Tests

The `tests` directory contains test programs to verify the functionality of the `fordiff` module. To run the tests using `fpm`, you can use response files for specific compilers:

```bash
fpm @test-<compiler>
```

`<compiler>: gfortran, ifx, ifort, nvfortran`

## TO DO
- [x] Complex-step: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Complex-step: f(x) f is a scalar-valued function and x is a vector variable
- [x] Complex-step: f(x) f is a vector-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a vector-valued function and x is a vector variable
- [ ] Automatic Differentiation

## API documentation

The most up-to-date API documentation for the main branch is available
[here](https://gha3mi.github.io/fordiff/).
To generate the API documentation for `ForDiff` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## Contributing

Contributions to fordiff are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.