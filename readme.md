## Inverted pendulum used for simulation

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
# Safe control
This section focus on the safety validation of the safe control scheme. The safety barrier is set as

```
x_1 cart position: [-15,15] m
x_2 cart position: [-10,10] m/s
x_3 pendulum angle: [-π/6,π/6] rad
x_4 pendulum angle velocity position: [5,5] rad/s
```

and the barrier weights are `[3, 5, 3, 5]`.

The initial values of the states are `[0; 10.5; 0; 0]`.
We first show that fix gain feedback control with control gain `K_control=[4.472; 21.750; 192.188; 158.755]` is not safe.

## Fix-gain Feedback control (not safe)
The first figure shows the states under this fix gain feedback controller.

![States with fixgain](/figs/States_with_fixgain.png)

The following figure shows that zeroing barrier function goes beyond 0, which means the safe barrier is exceeded.

![ZBF with fixgain](/figs/ZBF_with_fixgain.png)



## Safe control with oracle Feedback

## Safe control with secure Estimation as Feedback

The first figure shows the states with our proposed CBF controller with states estimated by the proposed secure estimator.

![States with attack and secure est](/figs/States_with_attack_and_secure_est.png)

The following figure shows that zeroing barrier function holds above 0, which means the states are in the safe barriers.

![ZBF with attack and secure est](/figs/ZBF_with_attack_and_secure_est.png)
