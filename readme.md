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
