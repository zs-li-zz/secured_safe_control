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
Q_in=Ts*Diagonal([0.1; 0.1; 0.01; 0.01])
R_in=Ts*Diagonal([0.1; 0.1; 0.01; 0.01])
Σ_in=Q_in

K_lqr= [-0.603524586303569	-1.67840806923021	-39.5143128608477	-9.72077388461956]
function LQGcontrol(X)
    u=-K_lqr*X # u is a scalar
    return u[1]
end
Λ, K, C, T = preprocess(A_in, B_in, C_in, Q_in, R_in, Σ_in)
# T is the transformation matrix from non-Diagonal to Diagonal
# K is the Kalman fixed gain
x0=[ 2 ; -1 ; 0.3 ; 0 ]

## get data
MAX_TIME = 1000

w=zeros(n,MAX_TIME)
v=zeros(m,MAX_TIME)

X=zeros(n,MAX_TIME)
Y=zeros(m,MAX_TIME)
Ya=zeros(m,MAX_TIME)

a=0*(rand(MAX_TIME).-0.5)

for k=1:MAX_TIME
    # original data
    w[:,k]=rand(Gaussian(zeros(n),Q_in))
    v[:,k]=rand(Gaussian(zeros(m),R_in))
    if k==1
        X[:,k]=x0
    else
        X[:,k]=A_in*X[:,k-1]+w[:,k-1]#+B_in*LQGcontrol(X[:,k-1])
    end
    Y[:,k]=C_in*X[:,k]+v[:,k]

end
# i=3
# time_axis=[0:MAX_TIME-1].*0.1
# state
# plot(time_axis, X[i,1:time_scale], label = "Oracle State", linecolor = "black", line = (:solid, 1))
# plot(time_axis, X[i,1:MAX_TIME], label = "my State", linecolor = "blue", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "my Estimation", linecolor = "blue", line = (:dot, 2))
# plot!(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "red", line = (:solid, 1))
# plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:dot, 2))


## get original kalman gain K_km

Xhat=zeros(n,MAX_TIME) # Xhat(:,k) means \hat{x} (k)
Xhatpost=zeros(n,MAX_TIME) # Xhatpost(:,k) means \hat{x} (k|k-1)
P=zeros(n,n,MAX_TIME) # P(:,:,k) means P(k)
Ppost=zeros(n,n,MAX_TIME) # Ppost(:,k) means \hat{x} (k|k-1)
Kmat=zeros(n,m,MAX_TIME)

for k in 1:MAX_TIME
    if k==1
        Ppost[:,:,k]=Σ_in
        Xhatpost[:,k]=zeros(n,1)
    else
        Ppost[:,:,k]=A_in*P[:,:,k-1]*A_in'+Q_in
        Xhatpost[:,k]=A_in*Xhat[:,k-1]#+B_in*LQGcontrol(Xhat[:,k-1])
    end
    Kmat[:,:,k]=Ppost[:,:,k]*C_in'*inv(C_in*Ppost[:,:,k]*C_in'+R_in)
    P[:,:,k]=Ppost[:,:,k]-Kmat[:,:,k]*C_in*Ppost[:,:,k]

    Xhat[:,k]=Xhatpost[:,k]+Kmat[:,:,k]*(Y[:,k]-C_in*Xhatpost[:,k])
end

Plim=P[:,:,MAX_TIME]
Pplus=A_in*Plim*A_in'+Q_in
K_km=Pplus*C_in'*inv(C_in*Pplus*C_in'+R_in)


## original kalman
time_scale=51
Xkm = zeros(n,time_scale) # real state
Xkm_hat = zeros(n,time_scale) # estimation
Ykm = zeros(m,time_scale) # real state
Ykma = zeros(m,time_scale)
for k=1:time_scale
    # real state
    if k==1
        Xkm[:,k]=x0
    else
        Xkm[:,k]=A_in*Xkm[:,k-1]+w[:,k-1]+B_in*LQGcontrol(Xkm_hat[:,k-1])
    end
    Ykm[:,k]=C_in*Xkm[:,k]+v[:,k]
    Ykma[:,k]=Ykm[:,k]
    Ykma[3,k]=Ykma[3,k]+a[k]
    if k==1
        Xkm_hat[:,k]=zeros(n,1)
    else
        Xkm_hat[:,k]=(A_in-K_km*C_in*A_in)*Xkm_hat[:,k-1]+(I-K_km*C_in)*B_in*LQGcontrol(Xkm_hat[:,k-1])+K_km*Ykma[:,k]
    end

end

## lasso
ζ=complex(zeros(mn,time_scale))
Ctrl=complex(zeros(mn,time_scale))
Xls = zeros(n,time_scale) # estimation
Xl = zeros(n,time_scale) # real state
Yl = zeros(m,time_scale) # real state
Yla = zeros(m,time_scale)

F=complex(zeros(n,n,m));
for i=1:m
    F[:,:,i]=V*Diagonal(V^(-1)*K[:,i])
end

for k=1:time_scale
    # real state
    if k==1
        Xl[:,k]=x0
    # elseif k<=5
    #     Xl[:,k]=A_in*Xl[:,k-1]+w[:,k-1]+B_in*LQGcontrol(Xkm_hat[:,k-1])
    else
        Xl[:,k]=A_in*Xl[:,k-1]+w[:,k-1]+B_in*LQGcontrol(Xls[:,k-1])
    end
    Yl[:,k]=C_in*Xl[:,k]+v[:,k]
    Yla[:,k]=Yl[:,k]
    # Yla[3,k]=Yla[3,k]+a[k]

    println("\n=========================================================")
    println("solving LASSO at k = ", k)
    # zeta estimator
    if k==1
        # ζ[:,k]=update_ζ( zeros(mn,1), Yla[:,k], 0 )
        ζ[:,k]=zeros(mn,1)
    # elseif k<=5
    #     ζ[:,k]=update_ζ( ζ[:,k-1], Yla[:,k], LQGcontrol(Xkm_hat[:,k-1]) )
    else
        ζ[:,k]=update_ζ( ζ[:,k-1], Yla[:,k], LQGcontrol(Xls[:,k-1]) )
    end
    # test_zero=zeros(n,1)
    # for i=1:m
    #     test_zero=test_zero+F[:,:,i]*ζ[(i-1)*n+1:i*n, k]
    # end
    # @show norm(test_zero-Xkm_hat[:,k], Inf)

    # solve opt problem
    γ = 10
        x, μ, ν = solve_opt(ζ[:,k], γ, 0)
    @show maximum(broadcast(abs, ν))
    if k<=5
        Xls[:,k] = Xkm_hat[:,k]
    else
        Xls[:,k] = T*x
    end

end


## transformed kalman
# Xtkm = zeros(n,time_scale)
# for k=1:time_scale
#     if k==1
#         Xtkm[:,k]=zeros(n,1)  # 初值随机生成
#     else
#         Xtkm[:,k]=(Λ-K*C*Λ)*Xtkm[:,k-1]+K*Ya[:,k]
#     end
# end
# Xtkm=T*Xtkm

##
# X=T_out*X
# Xl=T_out*Xl
# Xls=T_out*Xls
# Xkm=T_out*Xkm
# Xkm_hat=T_out*Xkm_hat

i=1
time_axis=[0:time_scale-1].*Ts
# compare state
# plot(time_axis, X[i,1:time_scale], label = "Oracle State", linecolor = "black", line = (:solid, 1))
# plot!(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "blue", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:dot, 2))
# compare error
plot(time_axis, Xkm[i,1:time_scale]-Xkm_hat[i,1:time_scale], label = "kalman est error", linecolor = "black", line = (:solid, 1))
plot!(time_axis, Xl[i,1:time_scale]-Xls[i,1:time_scale], label = "our est error", linecolor = "blue", line = (:dot, 2))
# compare state under control
# plot(time_axis, Xl[i,1:time_scale], label = "my State", linecolor = "blue", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "my Estimation", linecolor = "blue", line = (:dot, 2))
# plot!(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "red", line = (:solid, 1))
# plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:dot, 2))
