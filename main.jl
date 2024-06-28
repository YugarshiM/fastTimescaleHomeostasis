using DifferentialEquations, Random, MAT, NaNMath, Statistics
include("integrator.jl")

# 1 is restarting
# use this during self assembly to create burster, then turn off homeostatic mechanism
# when homeostatic mechanism if off, this burster may change
restarter = 1;


#==========================================================================================#
# ICs for Models 1-5 in the paper
#==========================================================================================#

#ICs = [-50.0 0.24 0.20 0.20 0.29 0.30 0.20 0.30 0.23 0.21 0.21 0.28 0.40 0.54 0.49 0.74 0.60 0.74 0.60 0.33 0.42 0.42 -0.15 -0.08 -0.23 0.14 -0.48 -0.10 0.46 0.16 0.17 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

#ICs = [-50.0 0.21 0.21 0.27 0.20 0.29 0.24 0.29 0.27 0.23 0.23 0.30 0.40 0.41 0.45 0.73 0.56 0.57 0.72 0.42 -0.23 0.01 0.44 0.40 0.25 0.25 0.20 -0.25 0.27 0.07 0.08 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

#ICs = [-50.0 0.20 0.21 0.25 0.29 0.28 0.29 0.24 0.26 0.20 0.28 0.21 0.40 0.83 0.34 0.87 0.65 0.51 0.77 0.63 -0.36 0.32 -0.38 0.39 -0.03 -0.37 -0.34 0.33 -0.11 0.32 -0.03 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

#ICs = [-50.0 0.29 0.25 0.25 0.22 0.30 0.30 0.22 0.23 0.21 0.22 0.23 0.40 0.58 0.32 0.61 0.62 0.56 0.62 0.49 0.26 0.37 0.16 0.43 -0.08 -0.08 -0.14 -0.40 0.37 -0.03 0.46 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

#ICs = [-50.0 0.21 0.24 0.22 0.30 0.30 0.23 0.21 0.21 0.23 0.23 0.29 0.40 0.75 0.35 0.57 0.43 0.61 0.73 0.45 -0.25 0.22 0.42 0.27 -0.11 -0.48 0.16 -0.41 -0.42 -0.10 -0.38 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

#==========================================================================================#


tgt1 = .25
tgt2 = .03 
tgt3 = .02 
gf = 53.0
gs = 3.0

for i = 1:1:3

	chunkSize = 3000.0;  #this is in seconds! (not milliseconds)
	tEnd = chunkSize*30;
	hsh = randstring(['A':'Z'; '0':'9'], 6);
	path = "./simulations/"
	fname = path*"tgt1."*string(tgt1)*".tgt2."*string(tgt2)*".tgt3."*string(tgt3)*"."*hsh;
	mkpath(fname)

	#ICs = [-50.0 round.(rand(Float64, (1, 11))*.1 .+ .2,digits = 2) .4 round.(rand(Float64, (1, 7))*.6 .+ .3,digits = 2) round.(-0.5*ones(1,11)+rand(Float64, (1, 11)),digits=2) 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

	# Or Insert Here Another IC from above
	ICs = [-50.0 0.24 0.20 0.20 0.29 0.30 0.20 0.30 0.23 0.21 0.21 0.28 0.40 0.54 0.49 0.74 0.60 0.74 0.60 0.33 0.42 0.42 -0.15 -0.08 -0.23 0.14 -0.48 -0.10 0.46 0.16 0.17 0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 1.0];

	print("\n RUN NUMBER: ",i, "\n")
	alphaVal = integrator(fname,tgt1,tgt2,tgt3,gf,gs,ICs,tEnd,0.0001,chunkSize,[600000.0, 6000.0, 2000.0, 2000.0])

	## restart to capture burster w/o homeostasis
	if restarter == 1
		chunkSize = 300.0; tEnd = chunkSize*30;
		## LOAD ICs previously run simulations
		file = matopen(fname*"/"*hsh*"_030.mat");
		ICs = read(file,"lastValue")
		close(file)
		##
		## rename for new run
		fname = path*"tgt1."*string(tgt1)*".tgt2."*string(tgt2)*".tgt3."*string(tgt3)*".RESTART."*hsh;
		mkdir(fname)
		ICs = round.(ICs,digits = 4)
		alphaVal = integrator(fname,tgt1,tgt2,tgt3,gf,gs,ICs,tEnd,0.0001,chunkSize,[Inf, Inf, 2000.0, 2000.0])
	end
end
