### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ b9ac6188-f68b-11ea-22a8-7dc92f1ff65e
using AliasedSeabed

# ╔═╡ f9aacfae-f68b-11ea-0384-19793b3c9ccb
using SimradEK60

# ╔═╡ be96809c-f68c-11ea-1488-bd84b262e45f
using EchogramPlots

# ╔═╡ bfd4a4e6-f690-11ea-0a5a-65af8e787377
using Colors

# ╔═╡ 12c55614-f68f-11ea-2a5e-ff590f917e88
using ColorSchemes

# ╔═╡ d1b07fe0-f694-11ea-191c-bb6141fb5416
md"# Aliased seabed detection demonstration"

# ╔═╡ 0b13af5a-f69f-11ea-0932-919b9d3afedf
md"Specify your Simrad EK60 RAW file here ..."

# ╔═╡ 5271d012-f68c-11ea-2310-e74aa00699c2
filename = "/home/reb/.julia/dev/AliasedSeabed/test/data/JR280_-D20121202-T170437.raw"
#filename = "/media/reb/gold/projects/fb-analysis/renfree/data/1301SH-D20130116-T232359.raw"
#filename = "/media/reb/gold/projects/fb-analysis/renfree/data/1301SH-D20130117-T000237.raw"
#filename = "/media/reb/gold/projects/fb-analysis/renfree/data/1301SH-D20130116-T234032.raw"

# ╔═╡ 64e44b58-f68c-11ea-0b14-d1391aa7cb35
ps = SimradEK60.load(filename);

# ╔═╡ 6a307878-f68c-11ea-2111-9b31ad9cbb8b
ps38 = [p for p in ps if p.frequency == 38000];

# ╔═╡ d93f7792-f68b-11ea-316a-d5388cd0e8e1
_Sv = Sv(ps38); # Volume backscatter

# ╔═╡ bf856e1c-f68e-11ea-1794-3d005b9cf095
ntheta = alongshipangle(ps38);

# ╔═╡ e21b213a-f68f-11ea-31ed-e9b46c3a7a73
nphi = athwartshipangle(ps38);

# ╔═╡ e324f26a-f69e-11ea-1871-81465679f79a
md"You might need to tweak the thresholds and window sizes below"

# ╔═╡ 697f3a7c-f696-11ea-13c1-11ac7c12041c
#Ttheta = 702
Ttheta = 702

# ╔═╡ 55320a56-f694-11ea-2f05-ebc235fb89b4
#Tphi=282
Tphi=282

# ╔═╡ 63464a8a-f694-11ea-1b3a-21da2e46541c
#dtheta=28
dtheta=28

# ╔═╡ 69634418-f694-11ea-23ad-036518d6e532
#dphi=52
dphi=52

# ╔═╡ 677daa64-f697-11ea-07df-a57cc0790621
md"The settings are Ttheta= $Ttheta, Tphi= $Tphi, dtheta= $dtheta, dphi= $dphi"

# ╔═╡ 1738661e-f694-11ea-0491-5dd5dd2da1b7
function showasb(Sv, a)
	foo = copy(Sv)
	foo[a] .=NaN
	echogram(foo, color=ColorSchemes.viridis.colors,title="ASB removed")
end

# ╔═╡ 8033d052-f691-11ea-2ee7-59d07f44288a
begin
	a = blackwell_asbmask(_Sv, ntheta, nphi, Ttheta=Ttheta, Tphi=Tphi,dtheta=dtheta,dphi=dphi)
	showasb(_Sv,a)
end

# ╔═╡ 23eeebc2-f69d-11ea-0ea4-eb5fb6692dec
echogram(_Sv, color=ColorSchemes.viridis.colors, vmin=-90,vmax=-50)

# ╔═╡ Cell order:
# ╟─d1b07fe0-f694-11ea-191c-bb6141fb5416
# ╠═b9ac6188-f68b-11ea-22a8-7dc92f1ff65e
# ╠═f9aacfae-f68b-11ea-0384-19793b3c9ccb
# ╠═be96809c-f68c-11ea-1488-bd84b262e45f
# ╠═bfd4a4e6-f690-11ea-0a5a-65af8e787377
# ╠═12c55614-f68f-11ea-2a5e-ff590f917e88
# ╟─0b13af5a-f69f-11ea-0932-919b9d3afedf
# ╠═5271d012-f68c-11ea-2310-e74aa00699c2
# ╠═64e44b58-f68c-11ea-0b14-d1391aa7cb35
# ╠═6a307878-f68c-11ea-2111-9b31ad9cbb8b
# ╠═d93f7792-f68b-11ea-316a-d5388cd0e8e1
# ╠═bf856e1c-f68e-11ea-1794-3d005b9cf095
# ╠═e21b213a-f68f-11ea-31ed-e9b46c3a7a73
# ╟─e324f26a-f69e-11ea-1871-81465679f79a
# ╠═697f3a7c-f696-11ea-13c1-11ac7c12041c
# ╠═55320a56-f694-11ea-2f05-ebc235fb89b4
# ╠═63464a8a-f694-11ea-1b3a-21da2e46541c
# ╠═69634418-f694-11ea-23ad-036518d6e532
# ╟─677daa64-f697-11ea-07df-a57cc0790621
# ╟─8033d052-f691-11ea-2ee7-59d07f44288a
# ╠═1738661e-f694-11ea-0491-5dd5dd2da1b7
# ╠═23eeebc2-f69d-11ea-0ea4-eb5fb6692dec
