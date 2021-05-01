## Inverted pendulum used for simulation

We use an inverted pendulum for simulation, the physical parameters are shown in the following figure:

<img src="/figs/inv_pen.png" width="480">


Based on the knowledge of theoretical mechanics, one obtains the inverted pendulum system dynamic:

<a href="https://www.codecogs.com/eqnedit.php?latex=(M&plus;m)\ddot{x}&plus;b_x&space;\dot{x}&space;&plus;ml&space;\ddot{\theta}&space;\cos&space;\theta&space;-&space;ml&space;\dot{\theta}^2\sin\theta=F" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(M&plus;m)\ddot{x}&plus;b_x&space;\dot{x}&space;&plus;ml&space;\ddot{\theta}&space;\cos&space;\theta&space;-&space;ml&space;\dot{\theta}^2\sin\theta=F" title="(M+m)\ddot{x}+b_x \dot{x} +ml \ddot{\theta} \cos \theta - ml \dot{\theta}^2\sin\theta=F" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=ml^2&space;\ddot{\theta}&plus;ml\ddot{x}\cos\theta&space;&plus;b_\theta&space;l&space;\dot{\theta}&space;=&space;mgl\sin\theta&space;." target="_blank"><img src="https://latex.codecogs.com/gif.latex?ml^2&space;\ddot{\theta}&plus;ml\ddot{x}\cos\theta&space;&plus;b_\theta&space;l&space;\dot{\theta}&space;=&space;mgl\sin\theta&space;." title="ml^2 \ddot{\theta}+ml\ddot{x}\cos\theta +b_\theta l \dot{\theta} = mgl\sin\theta ." /></a>

where $b_x$ and $b_\theta$ are the coefficients of friction w.r.t. cart and pendulum.

Rewrite it as a non-linear state space model :

<a href="https://www.codecogs.com/eqnedit.php?latex=\left[\begin{array}{cccc}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;m&plus;M&space;&&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;&&space;0&space;&&space;m&space;l^{2}&space;\end{array}\right]&space;\frac{d}{d&space;t}\left[\begin{array}{c}&space;x&space;\\&space;\dot{x}&space;\\&space;\theta&space;\\&space;\dot{\theta}&space;\end{array}\right]=\left[\begin{array}{c}&space;\dot{x}&space;\\&space;m&space;l&space;\dot{\theta}^{2}&space;\sin&space;\theta-b_x&space;\dot{x}&space;\\&space;\dot{\theta}&space;\\&space;m&space;g&space;l&space;\sin&space;\theta-b_\theta&space;l\dot{\theta}&space;\end{array}\right]&plus;\left[\begin{array}{c}&space;0&space;\\&space;1&space;\\&space;0&space;\\&space;0&space;\end{array}\right]&space;F." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[\begin{array}{cccc}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;m&plus;M&space;&&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;m&space;l&space;\cos&space;\theta&space;&&space;0&space;&&space;m&space;l^{2}&space;\end{array}\right]&space;\frac{d}{d&space;t}\left[\begin{array}{c}&space;x&space;\\&space;\dot{x}&space;\\&space;\theta&space;\\&space;\dot{\theta}&space;\end{array}\right]=\left[\begin{array}{c}&space;\dot{x}&space;\\&space;m&space;l&space;\dot{\theta}^{2}&space;\sin&space;\theta-b_x&space;\dot{x}&space;\\&space;\dot{\theta}&space;\\&space;m&space;g&space;l&space;\sin&space;\theta-b_\theta&space;l\dot{\theta}&space;\end{array}\right]&plus;\left[\begin{array}{c}&space;0&space;\\&space;1&space;\\&space;0&space;\\&space;0&space;\end{array}\right]&space;F." title="\left[\begin{array}{cccc} 1 & 0 & 0 & 0 \\ 0 & m+M & 0 & m l \cos \theta \\ 0 & 0 & 1 & 0 \\ 0 & m l \cos \theta & 0 & m l^{2} \end{array}\right] \frac{d}{d t}\left[\begin{array}{c} x \\ \dot{x} \\ \theta \\ \dot{\theta} \end{array}\right]=\left[\begin{array}{c} \dot{x} \\ m l \dot{\theta}^{2} \sin \theta-b_x \dot{x} \\ \dot{\theta} \\ m g l \sin \theta-b_\theta l\dot{\theta} \end{array}\right]+\left[\begin{array}{c} 0 \\ 1 \\ 0 \\ 0 \end{array}\right] F." /></a>

Reformulate the system as <a href="https://www.codecogs.com/eqnedit.php?latex=\dot{x}=f(x)&plus;g(x)u" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dot{x}=f(x)&plus;g(x)u" title="\dot{x}=f(x)+g(x)u" /></a> and linearize system at <a href="https://www.codecogs.com/eqnedit.php?latex=\theta=0,\dot{\theta}=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta=0,\dot{\theta}=0" title="\theta=0,\dot{\theta}=0" /></a>.
Then for continuous time $\dot{x}(t)=A_c x(t) +B_c u(t)$.
Sample the continuous time system with period $T_s$, then the discrete system is
$x(k+1)=Ax(k)+Bu(k)$ where

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." title="A=\exp(T\cdot A_c),\ B=\left(\int_{0}^{T} \exp(\tau A_c) d\tau \right)B_c ." /></a>

The system output equation is
$$Y(k)=Cx(k)+v(k)+a(k),\ C=\begin{bmatrix}1&0&0&0\\ 1&0&0&0\\ 0&0&1&0\\ 0&0&1&0 \end{bmatrix}$$
where $v(k)$ is the output noise and $a(k)$ is the injected attack.
There are 2 sensors monitoring the cart position and 2 sensors monitoring the pendulum angle.

The system is simulated with sampling time $T_s=1/200 s$.

The system noise process has covariance $Q=T_s^2 diag[0.1 0.1 0.01 0.01]$
and the measurement noise is $Q=T_s^2 diag[1 1 0.1 0.1]$,
i.e., the noise is scaled according to sampling time $T_s$.


# Safe control
This section focus on the safety validation of the safe control scheme.

The zeroing barrier function is ($n$ is the number of states)
$$h(x)=\sum_{i=1}^{n} w_i\left(1-\left(\frac{x_i}{A_i}\right)^2\right)$$
where $A_i$ is the safe barrier of state $x_i$, and $w_i$ is barrier weight. The states are in the safety region if $h(x)\geq 0$.

The safety barrier $A_i$ is set as

```
x_1 cart position: [-15,15] m
x_2 cart position: [-10,10] m/s
x_3 pendulum angle: [-π/6,π/6] rad
x_4 pendulum angle velocity position: [-1,1] rad/s^2
```

and the barrier weights $w_i$ are `[3, 5, 3, 5]`.

The initial values of the states are `[0; 10; 0; 0]`.
We first show that fix gain feedback control with control gain `K_control=[4.472; 21.750; 192.188; 158.755]` is not safe.

In all the following simulation experiments, if without explicitly stated, the feedback is provided by a fix gain Kalman estimator shown in the following:

<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{x}(k&plus;1)=A\hat{x}(k)&plus;K_{km}\left(y(k&plus;1)-CA\hat{x}(k)\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{x}(k&plus;1)=A\hat{x}(k)&plus;K_{km}\left(y(k&plus;1)-CA\hat{x}(k)\right)" title="\hat{x}(k+1)=A\hat{x}(k)+K_{km}\left(y(k+1)-CA\hat{x}(k)\right)" /></a>

## 1. Fix-gain Feedback control (not safe)
The first figure shows the states under this fix gain feedback controller.

![](/figs/States_lin_lin.png)

The following figure shows that zeroing barrier function goes below 0, which means the safe barrier is exceeded.

![](/figs/ZBF_lin_lin.png)


## 2. Safe control with linear estimator

In order to show that the performance of our proposed CBF controller, we implement our proosed controller with the fix gain estimator. The states and ZBF(zering barrier function) are shown in the following:

![](/figs/States_cbf_lin.png)

![](/figs/ZBF_cbf_lin.png)

Notice that even though the there are states that may exceed the barrier, as long as the weighted zeroing barrier function is greater than 0, the states are still in safety region.


## 3. Safe control with secure estimator in absence of sensor attack
In order to show that the performance of proposed estimator in absence of sensor attack. The system is simulated with secure estimator and begin sensors.

The first figure shows the states with our proposed CBF controller with states estimated by the proposed secure estimator.

![](/figs/States_cbf_sec.png)

The following figure shows that zeroing barrier function holds above 0, which means the states are in the safe barriers.

![](/figs/ZBF_cbf_sec.png)

# Secure Estimation
In order to validate the security of our proposed estimator. We design an attack on sensor 3 (angle sensor) and inject to $y_3$ a random value uniformly distributed on `[-π/2,π/2]`. Considering that $x_3$ is the pendulum angle, the measurement $y_3$ will be deviated severely from its original value. We will the power of this attack in the following.


## 1. Safe control with linear estimator
Even though with CBF safe controller, the system is not in the safe region because of corrupted estimation.

![](/figs/States_cbf_lin_att.png)

![](/figs/esterr_lin_att.png)

## 2. Safe control with secure estimator

The following figures show the states and zeroing barrier function with our proposed CBF controller and secure estimator.

![](/figs/States_cbf_sec_att.png)

The zeroing barrier functions stays above zero despite the attack.

![](/figs/ZBF_cbf_sec_att.png)

The estimation error of angle and angle velocity are smaller than linear estimator, which means our estimator works well.

![](/figs/esterr_sec_att.png)
