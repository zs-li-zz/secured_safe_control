#= coding: UTF-8
Julia 1.5.3
Apr 26th, 2021
@author: leo
This script simulates the linearlized inverted pendulum in discrete time
add noise, without the estimator
=#

using LinearAlgebra, GaussianDistributions, Random
using Plots;
include("CBF_control.jl")

m=1;M=1;l=1
In=m*l^2/3 # moment of inertia
g=9.8
b=0

Acon=[0 1 0 0 # A for continuous time
    0 -(In+m*l^2)*b/(In*(M+m)+M*m*l^2) m^2*g*l^2/(In*(M+m)+M*m*l^2) 0
    0 0 0 1
    0 -m*l*b/(In*(M+m)+M*m*l^2) m*g*l*(M+m)/(In*(M+m)+M*m*l^2) 0]

A_in=exp(Acon*0.1) # sampling every 0.1 s

B_in=[0; (In+m*l^2)/(In*(M+m)+M*m*l^2); 0; m*l/(In*(M+m)+M*m*l^2)]

C_in=[1 0 0 0
      1 0 0 0
      0 0 1 0
      0 0 1 0]

n=size(A_in,1)
m=size(C_in,1)
Q_in=0.0001*Matrix(1.0I,n,n)#+0.01*rand(n,n)
R_in= 0.001*Matrix(1.0I,m,m)#+0.01*rand(m,m)
Σ_in=Q_in

## get data
MAX_TIME = 100

w=zeros(n,MAX_TIME)
v=zeros(m,MAX_TIME)

X=zeros(n,MAX_TIME)
Y=zeros(m,MAX_TIME)
x0=[ 0 ; 0 ; 0 ; 0 ] # initial state

for k=1:MAX_TIME
    # original data
    w[:,k]=rand(Gaussian(zeros(n),Q_in)) # (第一个w没用到)
    v[:,k]=rand(Gaussian(zeros(m),R_in))
    if k==1
        X[:,k]=x0
    else
        X[:,k]=A_in*X[:,k-1]+w[:,k-1]+B_in*CBFControl(X[:,k-1])
    end
    Y[:,k]=C_in*X[:,k]+v[:,k]

end


time_axis=[0:MAX_TIME-1].*0.1
plot(time_axis, X[1,:], label = "position", linecolor = "blue", line = (:solid, 1))
plot!(time_axis, X[2,:], label = "velocity", linecolor = "blue", line = (:dot, 2))
plot!(time_axis, X[3,:], label = "angle", linecolor = "red", line = (:solid, 1))
plot!(time_axis, X[4,:], label = "angle velocity", linecolor = "red", line = (:dot, 2))
