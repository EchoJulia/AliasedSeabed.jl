# AliasedSeabed

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)

[![Build Status](https://travis-ci.org/robblackwell/AliasedSeabed.jl.svg?branch=master)](https://travis-ci.org/EchoJulia/AliasedSeabed.jl)

[![Coverage Status](https://coveralls.io/repos/robblackwell/AliasedSeabed.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/EchoJulia/AliasedSeabed.jl?branch=master)

[![codecov.io](http://codecov.io/github/robblackwell/AliasedSeabed.jl/coverage.svg?branch=master)](http://codecov.io/github/EchoJulia/AliasedSeabed.jl?branch=master)

## Introduction

A Julia implementation of the aliased seabed detection algorithm described in 
[Blackwell et al. (2019)](https://arxiv.org/abs/1904.10736).

## Example

```
using AliasedSeabed
using SimradEK60

# Load some 38 kHz data
filename = "JR280_-D20121202-T170437.raw"
ps = SimradEK60.load(filename)
ps38 = [p for p in ps if p.frequency == 38000]
_Sv = Sv(ps38) # Volume backscatter
ntheta = alongshipangle(ps38) # Split beam angle
nphi = athwartshipangle(ps38)

# Find aliased seabed
a = blackwell_asbmask(_Sv, ntheta, nphi)

```

