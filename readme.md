## Inverted pendulum used for simulation

We use an inverted pendulum for simulation, the physical parameters are shown in the following figure:

<img src="/figs/inv_pen.png" width="480">

Notice that the pendulum length is 5m long, which is harder to stablize than shorter pendulums and also pose challenge for controlling the states in safety barriers.

Based on the knowledge of theoretical mechanics, one obtains the inverted pendulum system dynamic:
<!-- $$(M+m)\ddot{x}+b \dot{x} +ml \ddot{\theta} \cos \theta - ml \dot{\theta}^2\sin\theta=F,$$ -->
<!-- $$ml^2 \ddot{\theta} +mgl\sin\theta =ml\ddot{x}\cos\theta.$$ -->
Rewrite it as a non-linear state space model :
<!-- $$
\left[\begin{array}{cccc}
	1 & 0 & 0 & 0 \\
	0 & m+M & 0 & m l \cos \theta \\
	0 & 0 & 1 & 0 \\
	0 & m l \cos \theta & 0 & m l^{2}
\end{array}\right] \frac{d}{d t}\left[\begin{array}{c}
	x \\
	\dot{x} \\
	\theta \\
	\dot{\theta}
\end{array}\right]=\left[\begin{array}{c}
	\dot{x} \\
	m l \dot{\theta}^{2} \sin \theta-b \dot{x} \\
	\dot{\theta} \\
	m g l \sin \theta
\end{array}\right]+\left[\begin{array}{c}
	0 \\
	1 \\
	0 \\
	0
\end{array}\right] F.
$$ -->
Reformulate the system as $\dot{x}=f(x)+g(x)u$ and linearize system at $\theta=0,\dot{\theta}=0$.
Then for continuous time $\dot{x}(t)=A_c x(t) +B_c u(t)$

<a href="https://www.codecogs.com/eqnedit.php?latex=A_c=\partial&space;f(x)/&space;\partial&space;x&space;|_{x_3=x_4=0}&space;=&space;\left[\begin{array}{cccc}&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;-\frac{b}{M}&space;&&space;-\frac{mg}{M}&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\\&space;0&space;&&space;\frac{b}{Ml}&space;&&space;\frac{g(M&plus;m)}{Ml}&space;&&space;0&space;\end{array}\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A_c=\partial&space;f(x)/&space;\partial&space;x&space;|_{x_3=x_4=0}&space;=&space;\left[\begin{array}{cccc}&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;-\frac{b}{M}&space;&&space;-\frac{mg}{M}&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\\&space;0&space;&&space;\frac{b}{Ml}&space;&&space;\frac{g(M&plus;m)}{Ml}&space;&&space;0&space;\end{array}\right]" title="A_c=\partial f(x)/ \partial x |_{x_3=x_4=0} = \left[\begin{array}{cccc} 0 & 1 & 0 & 0 \\ 0 & -\frac{b}{M} & -\frac{mg}{M} & 0 \\ 0 & 0 & 0 & 1 \\ 0 & \frac{b}{Ml} & \frac{g(M+m)}{Ml} & 0 \end{array}\right]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=B_c=g(x)|_{x_3=x_4=0}&space;=&space;\left[\begin{array}{c}&space;0&space;\\&space;\frac{1}{M}&space;\\&space;0&space;\\&space;-\frac{1}{Ml}&space;\end{array}\right]&space;." target="_blank"><img src="https://latex.codecogs.com/gif.latex?B_c=g(x)|_{x_3=x_4=0}&space;=&space;\left[\begin{array}{c}&space;0&space;\\&space;\frac{1}{M}&space;\\&space;0&space;\\&space;-\frac{1}{Ml}&space;\end{array}\right]&space;." title="B_c=g(x)|_{x_3=x_4=0} = \left[\begin{array}{c} 0 \\ \frac{1}{M} \\ 0 \\ -\frac{1}{Ml} \end{array}\right] ." /></a>

Sample the continuous time system with period $T$, then the discrete system is
$x(k+1)=Ax(k)+Bu(k)$ where

<a href="https://www.codecogs.com/eqnedit.php?latex=A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." target="_blank"><img src="https://latex.codecogs.com/gif.latex?A=\exp(T\cdot&space;A_c),\&space;B=\left(\int_{0}^{T}&space;\exp(\tau&space;A_c)&space;d\tau&space;\right)B_c&space;." title="A=\exp(T\cdot A_c),\ B=\left(\int_{0}^{T} \exp(\tau A_c) d\tau \right)B_c ." /></a>

The system is simulated with sampling time `Ts=1/200 s`.

The system noise process has covariance `Q=Ts^2*diag[0.1 0.1 0.01 0.01]` and the measurement noise is `Q=Ts^2*diag[1 1 0.1 0.1]`, i.e., the noise is scaled according to sampling time `Ts`.


# Safe control
This section focus on the safety validation of the safe control scheme. The safety barrier is set as

```
x_1 cart position: [-15,15] m
x_2 cart position: [-10,10] m/s
x_3 pendulum angle: [-π/6,π/6] rad
x_4 pendulum angle velocity position: [-1,1] rad/s^2
```

and the barrier weights are `[3, 5, 3, 5]`.

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
In order to validate the security of our proposed estimator. We design an attack on sensor 3 (angle sensor) and inject to `y_3` a random value uniformly distributed on `[-5,5]`.


## 1. Safe control with linear estimator
Even though with CBF safe controller, the system is not in the safe region because of corrupted estimation.

![](/figs/States_cbf_lin_att.png)

![](/figs/esterr_lin_att.png)

## 2. Safe control with secure estimator

The following figures show the states and zeroing barrier function with our proposed CBF controller and secure estimator.

![](/figs/States_cbf_sec_att.png)

![](/figs/ZBF_cbf_sec_att.png)

The estimation error of angle and angle velocity are smaller than linear estimator.

![](/figs/esterr_sec_att.png)
