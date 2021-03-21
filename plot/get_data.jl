using JLD
include("closeloop_est_error.jl")
MSE_data=zeros(2,15)
global i=0
for γ in [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500]
    global i
    i=i+1
    @time MSE_data[:,i]=get_MSE_data( γ )'
end
save("MSE_data_noattack.jld", "MSE_data",MSE_data)
