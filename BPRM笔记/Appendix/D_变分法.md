### Appendix D 变分法（Calculus of Variations）

对于单变量函数$y(x)$而言，在$x$处可以做泰勒展开
$$
y(x+\epsilon)=y(x)+\frac{dy}{dx}\epsilon+O(\epsilon^2) \tag{D.1}
$$
类似的可以对于多变量函数$y(x_1,...,x_D)$做泰勒展开
$$
y(x_1+\epsilon_1,...,x_D+\epsilon_D)=y(x_1,...,x_D)+\sum_{i=1}^{D}\frac{\part y}{\part x_i}\epsilon_i+O(\epsilon^2) \tag{D.2}
$$
对于泛函$F[y]$而言，可以看作无限维变量的函数在函数$y$定义域$\mathcal{D}$上做泰勒展开
$$
F[y+\epsilon\eta]=F[y]+\epsilon\int_\mathcal{D}\frac{\delta F}{\delta y}\eta(x)\mathrm{d}x+O(\epsilon^2) \tag{D.3}
$$
其中，泛函$F$的变分$\delta$定义为
$$
\delta F=F[y+\delta y]-F[y] \tag{D.4}
$$
令泛函的值在函数$y(x)$发生微小扰动$\epsilon\eta(x)$时，值不变，可得
$$
\int_{\mathcal{D}}\frac{\delta F}{\delta y}\eta(x)\mathrm{d}x=0 \tag{D.5}
$$
上式对任意的$\eta(x)$都成立，所以可得
$$
\frac{\delta F}{\delta y}\bigg|_\mathcal{D}=0 \tag{D.6}
$$


考虑下面的泛函
$$
F[y]=\int_\alpha^\beta G(y,y',x)\mathrm{d}x \tag{D.7}
$$
计算变分$\delta F$
$$
\begin{aligned}
\delta F
&=F[y+\delta y]-F[y] \\
&=\int_\alpha^\beta\{G(y+\delta y,(y+\delta y)',x)-G(y,y',x)\}\mathrm{d}x \\
&=\int_\alpha^\beta\{\frac{\part G}{\part y}\delta y+\frac{\part G}{\part y'}(\delta y)'\}\mathrm{d}x \\
&=\int_\alpha^\beta\frac{\part G}{\part y}\delta y\mathrm{d}x+\int_\alpha^\beta\frac{\part G}{\part y'}(\delta y)'\mathrm{d}x \\
&=\int_\alpha^\beta\frac{\part G}{\part y}\delta y\mathrm{d}x+\frac{\part G}{\part y'}\delta y\bigg|_{x=\alpha}^\beta-\int_\alpha^\beta\frac{\mathrm{d}}{\mathrm{d}x}\frac{\part G}{\part y'}\delta y\mathrm{d}x \\
&=\int_\alpha^\beta(\frac{\part G}{\part y}-\frac{\mathrm{d}}{\mathrm{d}x}\frac{\part G}{\part y'})\delta y\mathrm{d}x+\frac{\part G}{\part y'}\delta y\bigg|_{x=\alpha}^\beta \\
\end{aligned}
\tag{D.8}
$$
增加约束$y(x)$的两端固定，即$\delta y|_{x=\alpha}=\delta y|_{x=\beta}=0$，要使$F[y]$取得最大值，$y$的函数选取应该使得$\delta F$在$[\alpha,\beta]$间处处为0，即
$$
\delta F|_{x\in[\alpha,\beta]} = 0 \tag{D.9}
$$
解这个式子可以得到欧拉-拉格朗日方程
$$
(\frac{\part G}{\part y}-\frac{\mathrm{d}}{\mathrm{d}x}\frac{\part G}{\part y'})\bigg|_{x\in[\alpha,\beta]}=0 \tag{D.10}
$$
