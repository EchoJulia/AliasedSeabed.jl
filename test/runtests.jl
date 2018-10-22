using AliasedSeabed
using Test
using SimradEK60

filename = joinpath(dirname(@__FILE__), "data/JR280_-D20121202-T170437.raw")

ps = SimradEK60.load(filename)
ps38 = [p for p in ps if p.frequency == 38000]
_Sv = Sv(ps38) # Volume backscatter
ntheta = alongshipangle(ps38) # Split beam angle
nphi = athwartshipangle(ps38)

_m1 = m1(ntheta)
_m2 = m2(nphi)
a = blackwell_asbmask(_Sv, ntheta, nphi)

@test sum(_m1) == 66165
@test sum(_m2) == 39929
@test sum(a) == 178399
