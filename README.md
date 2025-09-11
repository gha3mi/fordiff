[![GitHub](https://img.shields.io/badge/GitHub-fordiff-blue.svg?style=social&logo=github)](https://github.com/gha3mi/fordiff)
[![Version](https://img.shields.io/github/v/tag/gha3mi/fordiff?label=version&sort=semver)](https://github.com/gha3mi/fordiff/releases)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/fordiff/)
[![Setup Fortran Conda CI/CD](https://github.com/gha3mi/fordiff/actions/workflows/CI-CD.yml/badge.svg?branch=main)](https://github.com/gha3mi/fordiff/actions/workflows/CI-CD.yml)
[![License](https://img.shields.io/github/license/gha3mi/fordiff?color=green)](https://github.com/gha3mi/fordiff/blob/main/LICENSE)


**ForDiff**: A Fortran library for numerical differentiation


- [fpm dependency](#fpm-dependency)
- [Usage](#usage)
- [Run Demos](#run-demos)
- [Run Tests](#run-tests)
- [TODO](#todo)
- [CI Status](#ci-status)
- [API documentation](#api-documentation)
- [Contributing](#contributing)

## fpm dependency

To use `ForDiff` as a dependency in your fpm project, include the following line in your `fpm.toml` file:

```toml
[dependencies]
fordiff = {git="https://github.com/gha3mi/fordiff.git"}
```

## Usage

```fortran
use fordiff
dfdx = derivative(f, x, h, method)
```

## Run Demos

`example/demo1.f90` demonstrates how to compute the derivative of a scalar-valued function w.r.t to a scalar variable using complex-step and finite-difference methods.

```bash
fpm run --example demo1
```

`example/demo2.f90` demonstrates how to compute the derivative of a scalar-valued function w.r.t to a vector variable using complex-step and finite-difference methods.

```bash
fpm run --example demo2
```

`example/demo3.f90` demonstrates how to compute the derivative of a vector-valued function w.r.t to a vector variable using complex-step and finite-difference methods.

```bash
fpm run --example demo3
```

## Run Tests

The `tests` directory contains test programs to verify the functionality of the `fordiff` module. To run the tests using `fpm`, you can use:

```bash
fpm test
```

## TODO
- [x] Complex-step: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Complex-step: f(x) f is a scalar-valued function and x is a vector variable
- [x] Complex-step: f(x) f is a vector-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a scalar variable
- [x] Finite Difference: f(x) f is a scalar-valued function and x is a vector variable
- [x] Finite Difference: f(x) f is a vector-valued function and x is a vector variable
- [ ] Automatic Differentiation

## CI Status

<!-- STATUS:setup-fortran-conda:START -->
<!-- STATUS:setup-fortran-conda:END -->

## API documentation

The most up-to-date API documentation for the main branch is available
[here](https://gha3mi.github.io/fordiff/).
To generate the API documentation for `ForDiff` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford README.md
```

## Contributing

Contributions to fordiff are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.