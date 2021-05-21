# ValidatedTransfer

[![Build Status](https://travis-ci.org/wormell/ValidatedTransfer.jl.svg?branch=master)](https://travis-ci.org/wormell/ValidatedTransfer.jl)

[![Coverage Status](https://coveralls.io/repos/wormell/ValidatedTransfer.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wormell/ValidatedTransfer.jl?branch=master)

Specialised package for validated transfer operator algorithms that work:

(a) with ```ValidatedNumerics``` intervals

(b) on parallel processes

(c) in Julia 0.4 (i.e. not using ```DistributedArrays```, which would be better).

See the following paper:

* Caroline L. Wormell (2017), Spectral Galerkin methods for transfer operators in uniformly expanding dynamics ([preprint](https://arxiv.org/abs/1705.04431))
