#=
with oracle feedback, used to quantify the estimation error


=#

using LinearAlgebra, GaussianDistributions, Random
using Plots; #pgfplotsx()  PGFPlotsX;
include("utility.jl")

m=1;M=1;l=1
g=9.8
b=0

Acon=[0 1 0 0 # A for continuous time
    0 -b/M -m*g/M 0
    0 0 0 1
    0 b/M*l g*(M+m)/M*l 0]
Bcon=[0; 1/M; 0; -1/M/l]
Ts=0.1
A_in=exp(Acon*Ts) # sampling every 0.1 s

B_in=[0.100000000000000	0.00500000000000000	-0.00164941490148451	-4.11010468697955e-05
0	0.100000000000000	-0.0498055805186480	-0.00164941490148451
0	0	0.103298829802969	0.00508220209373959
0	0	0.0996111610372960	0.103298829802969]*Bcon

C_in=[1 0 0 0
      1 0 0 0
      1 0 0 0
      0 0 1 0]

n=size(A_in,1)
m=size(C_in,1)
Q_in=Ts*Diagonal([0.01; 0.01; 0.01; 0.01])
R_in=Ts*Diagonal([0.01; 0.01; 0.01; 0.01])
Σ_in=Q_in

K_lqr= [-0.603524586303569	-1.67840806923021	-39.5143128608477	-9.72077388461956]
function LQGcontrol(X)
    u=-K_lqr*X # u is a scalar
    return u[1]
end
Λ, K_km, C, T = preprocess(A_in, B_in, C_in, Q_in, R_in, Σ_in)
# T is the transformation matrix from non-Diagonal to Diagonal
# K is the Kalman fixed gain
x0=[ 0 ; 1 ; 0 ; 1 ]

## get data
time_scale = 101

X=load("X.jld")["X"]
Y=load("Y.jld")["Y"]
w=load("w.jld")["w"]
v=load("v.jld")["v"]
K_km=load("K_km.jld")["K_km"]
Ya=Y#load("Y_attack4.jld")["Ya"]


function get_MSE_data(γ)

## original kalman
Xkm_hat = zeros(n,time_scale)
for k=1:time_scale
    if k==1
        Xkm_hat[:,k]=zeros(n,1)
    else
        Xkm_hat[:,k]=(A_in-K_km*C_in*A_in)*Xkm_hat[:,k-1]+(I-K_km*C_in)*B_in*LQGcontrol(X[:,k-1])+K_km*Ya[:,k]
    end

end

## lasso
ζ=complex(zeros(mn,time_scale))
Xls = zeros(n,time_scale)

# F=complex(zeros(n,n,m));
# for i=1:m
#     F[:,:,i]=V*Diagonal(V^(-1)*K[:,i])
# end

for k=1:time_scale
    # println("\n=========================================================")
    # println("solving LASSO at k = ", k)
    # zeta estimator
    if k==1
        ζ[:,k]=zeros(mn,1)
    # elseif k<=5
    #     ζ[:,k]=update_ζ( ζ[:,k-1], Yla[:,k], LQGcontrol(Xkm_hat[:,k-1]) )
    else
        ζ[:,k]=update_ζ( ζ[:,k-1], Ya[:,k], LQGcontrol(X[:,k-1]) )
    end
    # test_zero=zeros(n,1)
    # for i=1:m
    #     test_zero=test_zero+F[:,:,i]*ζ[(i-1)*n+1:i*n, k]
    # end
    # @show norm(test_zero-Xkm_hat[:,k], Inf)

    # solve opt problem
    x, μ, ν = solve_opt(ζ[:,k], γ, 0)
    # @show maximum(broadcast(abs, ν))
    # if k<=3
    #     Xls[:,k] = Xkm_hat[:,k]
    # else
    Xls[:,k] = T*x
    # end

end

MSEkm=norm(X[:,1:time_scale]-Xkm_hat[:,1:time_scale],2)^2/time_scale
MSEls=norm(X[:,1:time_scale]-Xls[:,1:time_scale],2)^2/time_scale
return [MSEkm, MSEls]
end

# i=3
# time_axis=[0:time_scale-1].*Ts
# # compare error
# # plot(time_axis, X[i,1:time_scale]-Xkm_hat[i,1:time_scale], label = "kalman est error", linecolor = "black", line = (:solid, 1))
# # plot!(time_axis, X[i,1:time_scale]-Xls[i,1:time_scale], label = "our est error", linecolor = "blue", line = (:dot, 2))
# # compare state under control
# plot(time_axis, X[i,1:time_scale], label = "State", linecolor = "black", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "my Estimation", linecolor = "blue", line = (:dot, 2))
# plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman State", linecolor = "red", line = (:solid, 1))
