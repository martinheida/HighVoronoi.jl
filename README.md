## Documentation
There currently is a problem with my github pages documentation. You are forwarded to my home institution.
# HighVoronoi.jl
A Julia Package for setting up high dimensional (i.e. any dimension >= 2) Finite Volume problems on Voronoi Meshes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.wias-berlin.de/people/heida/HighVoronoi.jl/index.html)
[![Stable(broken)](https://img.shields.io/badge/docs-stable-blue.svg)](https://martinheida.github.io/HighVoronoi.jl/stable/)
[![Dev(broken)](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinheida.github.io/HighVoronoi.jl/dev/)
[![Build Status](https://github.com/martinheida/HighVoronoi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/martinheida/HighVoronoi.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/martinheida/HighVoronoi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/martinheida/HighVoronoi.jl)
[![Coverage](https://coveralls.io/repos/github/martinheida/HighVoronoi.jl/badge.svg?branch=main)](https://coveralls.io/github/martinheida/HighVoronoi.jl?branch=main)

This document provides important information about the project's code coverage which are not mentioned in the manual. These are due to experimental features, and certain parts of the code that are not meant to be accessed under normal circumstances.

## Test Coverage

The project aims to maintain a high level of code coverage to ensure the quality and reliability of the codebase. However, you might notice that the test coverage percentage is not at its maximum. This is intentional and is a result of several factors:

1. **Experimental Features**: I have marked certain parts of the codebase as "experimental", also in the manual. These features are in a testing/developing phase and might not be fully stable, but I include them in order to enable the user to see where it might lead to. As a result, the coverage for these experimental areas is only 20% or lower.

2. **Unreachable Code**: Some sections of the code are intentionally not meant to be invoked under typical usage scenarios. These sections are designed to handle unexpected situations that should ideally never occur during normal execution and I am indeed not capable to invoke these parts by tests.

finally, there is less than 1% of code that is currently not invoked but this may change in a later version. 
