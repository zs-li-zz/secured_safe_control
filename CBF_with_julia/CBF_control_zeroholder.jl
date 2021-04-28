using LinearAlgebra, GaussianDistributions, Random
using Plots, DifferentialEquations, PyCall

ENV["PYTHON"]="/home/zs/miniconda3/bin/python" # your python exectuable file path
qpsolvers = pyimport("qpsolvers")
np = pyimport("numpy")

## PARAMETERS
# [-3.28001894978812	-3.56225366617832	22.2826696463037	7.71624059731720]
K_gain = [0.447213595500539, 3.10711058638360, 128.125101649409, 127.003953329490] .* [10, 7, 1.5, 1.25]
epsilon = 1e-6
gamma = 65
zeta = 40

g = 9.8
M = 1
m = 0.1
l = 5
bx = 0.05
bt = 0.1


barrWeights = [3, 5, 3, 5]

x_barr = 15
v_barr = 10
t_barr = 3*π/18
w_barr = 5

P = zeros(4, 4)
P[1, 1] = -2/(x_barr^2)*barrWeights[1]
P[2, 2] = -2/(v_barr^2)*barrWeights[2]
P[3, 3] = -2/(t_barr^2)*barrWeights[3]
P[4, 4] = -2/(w_barr^2)*barrWeights[4]

function fs(X)
    A=[1 0 0 0
    0 m+M 0 m*l*cos(X[3])
    0 0 1 0
    0 m*l*cos(X[3]) 0 m*l^2]
    B=[X[2]; m*l*sin(X[3])*X[4]^2 - bx*X[2]; X[4]; m*g*l*sin(X[3]) - bt*l*X[4]]
    return A^(-1)*B
end

function gs(X)
    A=[1 0 0 0
    0 m+M 0 m*l*cos(X[3])
    0 0 1 0
    0 m*l*cos(X[3]) 0 m*l^2]
    return A^(-1)*[0; 1; 0; 0]
end


function computeJacobian(X_op)
    N = length(X_op)
    J = zeros(N, N)
    Id = Matrix(1.0I, N, N)
    for i = 1:N
        J[i, :] = (fs(X_op + epsilon*Id[i, :]) - fs(X_op - epsilon*Id[i, :]))/2/epsilon
    end
    return J'
end


A_lin = computeJacobian(zeros(4))
B_lin = gs(zeros(4))

Lf_h(X) = (P*X)'*(A_lin*X)
Lg_h(X) = (P*X)'*B_lin
h(X) = (1 - (X[1]/x_barr)^2)*barrWeights[1] + (1 - (X[2]/v_barr)^2)*barrWeights[2] + (1 - (X[3]/t_barr)^2)*barrWeights[3] + (1 - (X[4]/w_barr)^2)*barrWeights[4]


function CBFControl(X)
    qp_P = np.array([1 + epsilon])
    qp_q = np.array([-K_gain'*X])
    qp_h = np.array([Lf_h(X)+gamma*h(X)-epsilon-zeta*Lg_h(X)*Lg_h(X)])
    qp_G = np.array([-Lg_h(X)-epsilon])
    ## using CVXOPT
    qp_result = qpsolvers.cvxopt_solve_qp(qp_P,qp_q,qp_G,qp_h)
    # print("\nqp_P, qp_q, qp_G, qp_h", qp_P, qp_q, qp_G, qp_h)
    print("\nqp_result: ", qp_result)

    return qp_result
end

function func(X, t, q)
    X_dot = zeros(8)
    i = int(t/Ts)
    U = CBFControl(X[5:end])
    X_dot[1:4] = fs(X[1:4])+gs(X[1:4])*U
    X_dot[5:end] = f_est( X[5:end], f_out(X[1:4]) +  0.01*q*v[i]*U )
    #U = K(X + 30*noise[i]*np.array([1, 1, 1, 1]))
    return X_dot
end

function invpen_dynamic(X, p, t)
    U = CBFControl(X)
    X_dot = fs(X)+gs(X).*U
    return X_dot
end

## begin running simulation
C =  [1 0 0 0
      0 0 1 0]
n=4
m=size(C,1)
Ts=0.02
Q=Ts*0.0001*Matrix(1.0I,n,n)#+0.01*rand(n,n)
R= Ts*0.001*Matrix(1.0I,m,m)#+0.01*rand(m,m)

## Calculate discrete system Matrix
A_dis=exp(Ts*A_lin)
B_dis=[0.02 0.0001999334 -0.0000013570 0.0000001377
0 0.0199900046 -0.0001956617 0.00019997213
0 0.0000000371 0.0200028716 0.0001997213
0 0.0000019965 0.0004305862 0.0199589342]*B_lin

L=[1.7 0.1
0.541 -0.832
0.01 1.53
0.0148 2.4454]

## get data
MAX_TIME = 500

w=zeros(n,MAX_TIME)
v=zeros(m,MAX_TIME)

X=zeros(n,MAX_TIME)
Y=zeros(m,MAX_TIME)
X_est=zeros(n,MAX_TIME)
U=zeros(MAX_TIME)
# initialization
X[:,1]=[ 0 ; 11.856 ; 0 ; 0 ]
Y[:,1]=C*X[:,1]+Ts*rand(Gaussian(zeros(m),R))
U[1]= CBFControl(X_est[:,1])[1]

for k=2:MAX_TIME
    # original data
    w[:,k]=Ts*rand(Gaussian(zeros(n),Q)) # (第一个w没用到)
    v[:,k]=Ts*rand(Gaussian(zeros(m),R))
    X[:,k]=A_dis*X[:,k-1]+B_dis*U[k-1] # now without noise w
    Y[:,k]=C*X[:,k]+v[:,k]
    X_est[:,k]=A_dis*X_est[:,k-1]+B_dis*U[k-1]+L*(Y[:,k]-C*X_est[:,k-1])
    U[k]=CBFControl(X_est[:,k])[1]
end


# X = np.zeros([N, 4])
# Y = np.zeros([N, 2])
# X_est = np.zeros([N, 4])
# X[0, :] = initial
# Y[0, :] = f_out(X[0, :]) + Ts * q * noise[0] * np.array([1, 1]) # add output noise
# U[0] = CBFControl(X[0, :])
# #X_est at time 0 is zero
# for i in range(1,N):
#       # calculate control with states at last time step
#     X[i, :] = np.matmul(A_dis, X[i-1, :]) + B_dis * U[i-1]
#     Y[i, :] = f_out(X[i, :]) + Ts * q * noise[i] * np.array([1, 1]) # add output noise
#     X_est[i, :] = f_est(X[i - 1, :], Y[i, :], U[i-1])
#     U[i] = CBFControl(X[i, :])
#     print("step num: ", i)
#     print("X_RESULT: ", X[i, :])


time_axis=[0:MAX_TIME-1].*Ts
plot(time_axis, X[1,:]-X_est[1,:], label = "position", linecolor = "blue", line = (:solid, 1))
plot!(time_axis, X[2,:]-X_est[2,:], label = "velocity", linecolor = "blue", line = (:dot, 2))
plot!(time_axis, X[3,:]-X_est[3,:], label = "angle", linecolor = "red", line = (:solid, 1))
plot!(time_axis, X[4,:]-X_est[4,:], label = "angle velocity", linecolor = "red", line = (:dot, 2))

# h_test = (1 .- (X[1,:] / x_barr).^2) * barrWeights[1] + (1 .- (X[2,:] / v_barr).^2) * barrWeights[2] + (
#             1 .- (X[3,:] / t_barr).^2) * barrWeights[3] + (1 .- (X[4,:] / w_barr).^2) * barrWeights[4]
# plot(time_axis,h_test)
