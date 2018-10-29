module AliasedSeabed

export blackwell_asbmask, m1, m2

using Statistics
using EchogramUtils
using EchogramProcessing

function m1(ntheta; Ttheta=702, dtheta=28)
    meanfilter(ntheta, dtheta, dtheta).^2 .> Ttheta
end

function m2(nphi; Tphi=282, dphi=52)
    meanfilter(nphi, dphi, dphi).^2  .> Tphi
end

function myfindn(x)
    I = findall(!iszero, x)
    (getindex.(I, 1), getindex.(I, 2))
end

"""

Implements aliased seabed detection according to the author's
forthcoming paper.

"""
function blackwell_asbmask(Sv, ntheta, nphi; Ttheta=702, Tphi=282, dtheta=28, dphi=52, minSv=nothing)

    _m1 = m1(ntheta, Ttheta=Ttheta, dtheta=dtheta)
    _m2 = m2(nphi, Tphi=Tphi, dphi=dphi)
             
    m = (_m1 .| _m2)
    
    a = Sv[m]

    b = filter(x -> x !=-999, a) # remove 'missing' cells

    if length(b) < 1
        # No ASB found
        return map(x -> false, Sv)
    end
    
    md = median(b)

    # info("md $md")

    if minSv == nothing
        T = md
    else
        T = max(md, minSv)
    end

    # info("T $T")
          
    aboveT = Sv.> T # Points above lower threshold. 
    rows, cols = myfindn(m)

    bw = bwselect(aboveT, cols, rows)

    bw .| m

end




end # module
