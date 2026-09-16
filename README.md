# Statistical Computing

This is folder to store code files for the courese "Statistical Computing"

**Assignment**
- There are four assignment total in this semester

**Presentation**
- There are two presentation/report during the semester

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
    - Discrete: 因為離散型的pmf是step function，所以可能會發生$x_k$到$x_{k+1}$之間所對應的cdf都是一樣的。因此，$F(X)=u$的反函數$F_{X}^{-1}(u)$被定義為能使$F(X)$到達u最小的那個x。
- The way to prove this method work -> 證明$P(F_X^{-1}(U)≤x)=F_X(x)$
- Example:
    - Ex2: 生成pdf為 $f_X(x) = 3x^2, 0 < x < 1$ 的random sample
    - Ex3: 生成服從exponential distribution的random sample



