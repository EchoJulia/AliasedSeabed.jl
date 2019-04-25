module AliasedSeabed

export blackwell_asbmask, m1, m2

using Statistics
using ImageFiltering

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

function meanfilter(A, height, width)
    kernel  = Array{Float64}(undef, height, width)
    kernel .= 1 / (width * height)
    kernel = centered(kernel)
    imfilter(A, kernel, Fill(0))
end

"""
    bwselect(BW, c , r)
Select objects in a binary image.
Similar to the MATLAB function of the same name.
"""
function bwselect(BW, c, r)
    # constants
    north = CartesianIndex(-1,  0)
    south = CartesianIndex( 1,  0)
    east  = CartesianIndex( 0,  1)
    west  = CartesianIndex( 0, -1)
    
    queue = CartesianIndex.(r,c)
    
    m,n = size(BW)
    out = falses(m,n)
    
    while !isempty(queue)
        node = pop!(queue)
        
        if BW[node]
            wnode = node
            enode = node + east
        
             # Move west until node is false
            while checkbounds(Bool, BW, wnode) && BW[wnode] 
                out[wnode] = true
                if checkbounds(Bool, BW, wnode + north) && !out[wnode + north]
                    push!(queue, wnode + north)
                end
                if checkbounds(Bool, BW, wnode + south) && !out[wnode + south]
                    push!(queue, wnode + south)
                end
                wnode += west
            end
        
            # Move east until node is false
            while checkbounds(Bool, BW, enode) && BW[enode]
                out[enode] = true
                if checkbounds(Bool, BW, enode + north) && !out[enode + north]
                    push!(queue, enode + north)
                end
                if checkbounds(Bool, BW, enode + south) && !out[enode + south]
                    push!(queue, enode + south)
                end
                enode += east
            end
        end
        
    end
    return out
end

"""

Aliased seabed detection according to [Blackwell et al. (2019)](https://arxiv.org/abs/1904.10736)

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

    @info("medianSv = $md")

    if minSv == nothing
        T = md
    else
        T = max(md, minSv)
    end

    @info("T = $T")
          
    aboveT = Sv.> T # Points above lower threshold. 
    rows, cols = myfindn(m)

    bw = bwselect(aboveT, cols, rows)

    bw .| m

end

end # module
