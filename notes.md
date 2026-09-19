## Lecture 1 (chap 3)
**Abbreviation**
- d: density
- p: cdf
- q: quantile
- r: random number from the distribution

**The Inverse Transform Method**
- Every random variable $U$ from a cdf $F_{x}(x)$ is going to follow the distribution $Uniform(0,1)$
- 方法：\
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
- Example:
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

**Acceptance-Rejection Method**
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
- Example:
    - Ex1: Generate triangular distribution
    - Ex2: Generate gamma disribution
- A-R Method for discrete r.v.
    - 利用一個比較方便生成的Y，近一步生成X
    1. Generate X, with corresponding $q_i$
    2. Generate $u \sim Unif(0, 1)$
    3. if $u < \frac{p_j}{cq_j}$, set $X = x_j$

**Transformation Method**
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
- Example:
 - Ex1: 生成Beta distribution
 - Ex2: 生成logartihmic distribution
 






