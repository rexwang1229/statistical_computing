# Lecture 1 (chap 3)
## Abbreviation
- d: density
- p: cdf
- q: quantile
- r: random number from the distribution

## The Inverse Transform Method
Every random variable $U$ from a cdf $F_{x}(x)$ is going to follow the distribution $Uniform(0,1)$
**方法：**\
    **Continuous**
    1. 生成一組服從 $Unif(0,1)$的隨機變數$u$
    2. 找到想要生成分配cdf的反函數，再將這些前面生成的均勻分布隨機變數代入
    
    **Discrete**
    1. 生成一組服從 $Unif(0,1)$的隨機變數$u$
    2. 藉由$F_X(x_{i-1})< u ≤ F_X(x_i)$決定這個u在轉換後變成$x_i$，也就是我們想要的random sample

- Inverse transformation: $F_{X}^{-1}(u)=inf\{x:F_{X}(x)=u\},0<x<1$
    - Continuous: 直接對應，若$F(x_1)=u_1$，則反函數$F^{-1}(u_1)=x_1$
    - Discrete: 因為離散型的cdf是step function，所以可能會發生$x_k$到$x_{k+1}$之間所對應的cdf都是一樣的。因此，$F(X)=u$的反函數$F_{X}^{-1}(u)$被定義為能使$F(X)$到達u最小的那個x。
- The way to prove this method work -> 證明$P(F_X^{-1}(U)≤x)=F_X(x)$

**Example:**
- Ex2: 生成pdf為 $f_X(x) = 3x^2, 0 < x < 1$ 的random sample
- Ex3: 生成服從exponential distribution的random sample
- EX4: 生成服從geometric distribution的random sample
    - geometric distribution的cdf: $F(x)=1-q^{x+1}$
    - 利用$F_X(x_{i-1}) < u ≤ F_X(x_i)$的觀念來做
    - 可以得到$1-q^x < u ≤ 1 - q^{x+1}$的式子
    - 推導後可以得到$x + 1 = [log(1-u)/log(q)]$
    - R實作方法：
        - floor: 向下取整（公式不用再扣1）
        - ceiling: 向上取整（直接用x+1的那個下去代）
- Ex5: 生成服從Poisson distribution的random sample
    - 利用 $f(x+1) = \frac{\lambda f(x)}{x+1}; F(x+1) = F(x) + f(x+1)$來做
    - 先定義一個for迴圈包含所有的random sample
    - 再用一個while loop處理個別sample變換的狀況
- Ex6: Logarithmic Distribution （對數分佈）

## Acceptance-Rejection Method
- 使用時機：當inverse transformation不好做的時候（cdf反函數不好求）
- Target density $f(x)$: 我們想要生成random sample的分佈
- Instrumental density $g(x)$: 比較簡單的分佈，用其來幫助目標分布的生
- Concept: 
    - X and Y are random variables with pmf/pdf $f$ and $g$ 
    - Goal: 讓這兩個分布滿足 $\frac{f(t)}{g(t)} ≤ c$ 這樣的條件
- Procedure
    1. 生成服從instrumental density $g$ 的 random sample $y$
    2. 生成uniform distribution的random sample $u$
    3. 若 $u < \frac{f(y)}{cg(y)}$，則定 $x = y$，若不符合則回到第一步重新找一個
- A-R Method for discrete r.v.
    - 利用一個比較方便生成的Y，近一步生成X
    1. Generate X, with corresponding $q_i$
    2. Generate $u \sim Unif(0, 1)$
    3. if $u < \frac{p_j}{cq_j}$, set $X = x_j$
**Example:**
- Ex1: Generate triangular distribution
- Ex2: Generate gamma disribution


## Transformation Method
- 利用不同變數之間的變換，生成更多不同分布的變數
- $Z \sim N(0, 1)$ 系列
    - Chi-Squared: $V = Z^2 \sim \chi^2(1)$
    - F: $F = \frac{U / m}{V / n}$ if $U \sim \chi^2(m)$ and $V \sim \chi^2(n)$ are independent
    - t(n): $T = \frac{Z}{\sqrt{V / n}}$ if $Z \sim N(0, 1)$ and $V \sim \chi^2(n)$ are independent
    - Box-Muller algorithm: \
    If $U, V \sim Unif(0, 1)$ and independent, then

        $Z_1 = \sqrt{-2logU} cos(2 \pi V)$ and $Z_2 = \sqrt{-2logU} sin(2 \pi V)$ 
    
        are independent standard normal variables.
- Gamma系列
    - If $U \sim \Gamma(r, \lambda) $ and $V \sim \Gamma(s, \lambda)$ are independent, \
    then $X = \frac{U}{U + V} \sim Beta(r, s)$ 
    - If $U, V \sim Unif(0, 1)$ are independent, then\
    $X = [1 + \frac{log(V)}{log(1 - (1 - \theta)^U)}] \sim logarithemic(\theta)$
**Example:**
- Ex1: 生成Beta distribution
- Ex2: 生成logartihmic distribution
 
## Convolution and Mixture
**Convolution $F^*(n)_X$**
- $X_1, \dots, X_n$ is iid distributed
- $S = X_1 + \dots + X_n$ is called the n-fold convolution of X
- Example: Chi-square distribution
    1. 製作一個 $n * v$ 的矩陣，矩陣內填滿了服從N(0, 1)的r.v.
    2. 對矩陣內所有element取平方，並計算每個row的加總
    3. 每個row的加總形成的vector即為我們想要的chi-square random sample

**Mixture**
- 一個分佈是很多個分佈的組合
- $F_X(x) = \int_{-\infty}^{\infty} F_{X|Y=y}(x)f_Y(y) \, dy$ and $\int_{-\infty}^{\infty}f_Y(y)\,dy = 1$
- 如何在R生成一個mixture distribution
    1. 生成$k \in \{1, 2\}$，P(1) = P(2) = 0.5
    2. 若k = 1，生成一個x服從$f_{X_1}(x_1)$，反之則生成一個x服從$f_{X_2}(x_2)$
- Example:
    - 生成兩個Normal的mixture
    - 生成兩個gamma的mixture
    - 生成多個gamma的mixture 
        - $F_X = \sum_{i=1}^5 \theta_i F_{X_i}$
        - $X_i \sim \Gamma(r = 3, \lambda_i = 1/i)$, $\theta_i = i/15, i = 1,2,\dots,5$
            - $X_1 \sim \Gamma(r = 3, \lambda_i = 1)$, $\theta_1 = 1/15$
            - $X_2 \sim \Gamma(r = 3, \lambda_i = 1/2)$, $\theta_1 = 2/15$
            - $X_3 \sim \Gamma(r = 3, \lambda_i = 1/3)$, $\theta_1 = 3/15$
            - $X_4 \sim \Gamma(r = 3, \lambda_i = 1/4)$, $\theta_1 = 4/15$
            - $X_5 \sim \Gamma(r = 3, \lambda_i = 1/5)$, $\theta_1 = 5/15$
    - Poisson-Gamma mixture distribution

## Multivariate Normal Distribution

If $Z \sim N_d(0, I_d)$, then $CZ + b \sim N_d(b, CC^T)$ \
Supposed $\Sigma$ can be factored so that $\Sigma = CC^T$ for some matrix $C$. \
Then $CZ + \mu \sim N_d(\mu, \Sigma)$
- 目標: 想辦法對 $\Sigma$ 做factorization，進一步得到$\Sigma^{1/2}$
- 方法：
    1. spectral decomposition -> $\texttt{eigen()}$
    2. Choleski factorization -> $\texttt{chol()}$
    3. singular value decomposition (SVD) -> $\texttt{svd()}$

**Algorithm for generating random sample from $N_d(\mu, \Sigma)$**
1. 生成一個 $n*d$ 的矩陣$Z$，矩陣中都是$N(0, 1)$的隨機樣本
2. 對 $\Sigma$ 作factorization，得到 $\Sigma = Q^TQ$
3. 對Z進行轉換 $X = ZQ + 1_n\mu^T$
4. 轉換後X的每個row即為服從 $N_d(\mu, \Sigma)$ 的一個樣本

**Spectral Decomposition Method**

- $ \Sigma = P \Lambda P^T$
- $\Lambda$: 對角矩陣，矩陣內的元素為 $\Sigma$ 的eigenvalue
- $P$: eigenvector matrix
- 利用Spectral Decomposition Method將 $\Sigma$ 拆成 $Q^TQ$ 

    $\Sigma^{1/2} = P \Lambda^{1/2} P^T = Q \Rightarrow \Sigma = Q^TQ$

**Singular Value Decomposition (SVD)**

For any $m \times n$ matrix $X$
- There exist orthogonal matrix $U$ of size $m \times m$ and $V$ of size $n \times n$ such that $X = UDV^T$
- $D$: a $m \times n$ matrix with nondiagonal entries all zero

$\Rightarrow U = V = P$ and $\Sigma^{1/2} = UD^{1/2}V^T$ 

**Choleski Method**

For a symmetric positive-definite matrix X, there exists an upper triangular matrix $Q$ such that $X = QQ^T$

## Wishart Distribution
**General Method**

先生成n個服從 $N_d(\mu, \Sigma)$ 的random sample，再利用$M = X^TX$的轉換得到Wishart的random sample
1. $M = X^TX$, where $X$ is an $n \times d$ data matrix of a random sample form a $N_d(0, \Sigma)$ distribution
2. Then $M \sim W_d(\Sigma, n)$, a Wishart distribution with scale matrix $\Sigma$ and $n$ degree of freedom
3. When d = 1, $X \sim N(0, \sigma^2), W_1(\sigma^2, n) = \sigma^2\chi^2(n)$

**Bartlett's Decomposition**

Let $A = (A_{ij})$ be the lower triangular $d \times d$ random matirx with independent entries satisfying
- $A_{ij} \stackrel{iid}{\sim} N(0, 1), i>j$
- $A_{ii} \sim \sqrt{\chi^2(n-i+1)}$, for $i = 1,\dots, d$.

作法：
1. Obtian Choleski decomposition $\Sigma = LL^T$, where $L$ is 下三角矩陣且對角的值非負
2. Then $LAA^TL^T \sim W_d(\Sigma, n)$



# Leture 2
## Simulating Survival Data
### Cox Propotional Hazards Model

Survival analysis建立在hazard function $h_i(t)$ 以及survival function $S(t|X_i)$上 
- Hazard function: $h_i(t) = h_0(t)exp(X_i\beta)$
    - 概念：在t這一時間，瞬間發生事件的風險
        - $h(t) = \lim_{\Delta t \to \infty}\frac{P(t ≤ T ≤ t+\Delta t| T≥t)}{\Delta t}$
    - $h_0(t)$: baseline hazard
    - $X_i$: baseline covariate
    - $\beta$: vector of log hazard ratio
- Survival function: $S(t|X_i) = exp(-H_0(t)exp(X_i\beta))$
    - 概念: 個體存活超過t的機率 $S(t) = P(T > t)$
    - $H_0(t)$: Cumulative baseline hazard funciton
    - Proportional Hazard: 同一個風險函數下，survival rate只受$X_i$影響，和$t$無關

### Simulate the survival data
Simulate survival times $T_i$ under a Cox Propotional Hazard model
1. Set $S(t|X_i) = exp(-H_0(t)exp(X_i\beta)) = U_i$, where $U_i \sim Unif(0, 1)$
2. $H_0(t) = \frac{-ln(U_i)}{exp(X_i\beta)} \Rightarrow T_i = H_0^{-1}(\frac{-ln(U_i)}{exp(X_i\beta)})$ 

Hazard function $h(t)$ 是自己訂的，$T_i$對於幾個比較常見的hazard function具有close form
- **Exponential distribution**: constant baseline hazard 

    $h_0(t) = \lambda \Rightarrow H_0(t) = \lambda t$ 

    $T_i = \frac{-ln(U_i)}{\lambda exp(X_i \beta)}$

- **Weibull distribution**: monotonic baseline hazard 
    
    $h_0(t) = \lambda \nu t^{\nu-1}\Rightarrow H_0(t) = \lambda t^\nu$ 

    $T_i = (\frac{-ln(U_i)}{\lambda exp(X_i \beta)})^{1/\nu}$
- **Gompertz distribution**: exponentially changing baseline hazard
    
    $h_0(t) = a\exp(bt) \Rightarrow H_0(t) = \frac{a}{b}(\exp(bt) - 1)$ 

    $T_i = \frac{1}{b} \ ln(1 - \frac{b\ln(U_i)}{\lambda \exp(X_i \beta)})$

若是遇到hazard function太複雜沒有close-form的狀況\
$\Rightarrow$ 利用$\texttt{simsurv}$ to solve the root-finding problem numerically

## Censoring Data
受試者在實驗進行到一半時退出的狀況\
For each individual $i$
- $T_i$: True event time
- $C_i$: Censoring time
- $d_i$: Observed status indicator $d_i = I(T_i ≤ C_i)$
    - $d_i = 1$: Uncensored
    - $d_i = 0$: Censored
- Censor Rate = $P(T_i > C_i) = \frac{\sum_{i=1}^{n}I(d_i = 0)}{n} = 1-\bar{d}$

### Simulation of Right-Censoring Mechanisms
**1. Independent Censoring Time Generation**
- 對每個資料i個別模擬一個censoring time $C_i$
- $C_i \sim Uniform(0, \tau)$
- 透過調整upper bound $\tau$，可以控制資料中censored data的比例

**2. Constructing the Observed Variables**
- 每個survival dataset都會有兩個observed variables $(Y_i, d_i)$
- Observed Survival TIme $Y_i = min(T_i, C_i)$  
- Event Indicator $d_i = I(T_i ≤ C_i)$
    - 0: Event was actually observed
    - 1: Individual was censored

**如何控制censor rate**
- 透過調整$\tau$來控制，$\tau$越大則more censoring，反之亦然
- 找到最佳$\tau$的方式
    - Mathematical expection
    - Iterative simulation / Root-finding

### Models for Survival Data
#### Cox Proportional Hazards Model (Cox PH)




#### Accelerated Failure Time (AFT)
