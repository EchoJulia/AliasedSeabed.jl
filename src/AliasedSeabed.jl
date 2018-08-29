module AliasedSeabed

export asbmask

using EchogramUtils
using EchogramProcessing

function asbmask(Sv, ntheta, nphi; Ttheta=702, Tphi=282, dtheta=28, dphi=52, minSv=nothing)
    
    m1 = meanfilter(ntheta, dtheta, dtheta).^2 .> Ttheta  
    m2  = meanfilter(nphi, dphi, dphi).^2  .> Tphi

    m = (m1 .| m2)

    a = Sv[m]

    b = filter(x -> x !=-999, a) # remove 'missing' cells

    md = median(b)

    info("md $md")

    if minSv == nothing
        T = md
    else
        T = max(md, minSv)
    end

    info("T $T")
          
    aboveT = Sv.> T # Points above lower threshold. 
    rows, cols = findn(m)

    bw = bwselect(aboveT, cols, rows)

    bw .| m

end



end # module
