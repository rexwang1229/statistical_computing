# Lecture 1 (Chap 3)

## Abbreviation

| 前綴 | 意義 |
| --- | --- |
| `d` | density |
| `p` | cdf |
| `q` | quantile |
| `r` | random number from the distribution |

---

## The Inverse Transform Method

Every random variable $U$ from a cdf $F_{X}(x)$ is going to follow the distribution $Uniform(0,1)$.

### 方法

**Continuous**

1. 生成一組服從 $Unif(0,1)$ 的隨機變數 $u$
2. 找到想要生成分配 cdf 的反函數，再將這些前面生成的均勻分布隨機變數代入

**Discrete**

1. 生成一組服從 $Unif(0,1)$ 的隨機變數 $u$
2. 藉由 $F_X(x_{i-1}) < u \le F_X(x_i)$ 決定這個 $u$ 在轉換後變成 $x_i$，也就是我們想要的 random sample

### 觀念

- Inverse transformation：

  $$F_{X}^{-1}(u)=\inf\{x:F_{X}(x)=u\},\quad 0<u<1$$

  - **Continuous**：直接對應，若 $F(x_1)=u_1$，則反函數 $F^{-1}(u_1)=x_1$
  - **Discrete**：因為離散型的 cdf 是 step function，所以可能會發生 $x_k$ 到 $x_{k+1}$ 之間所對應的 cdf 都是一樣的。因此，$F(X)=u$ 的反函數 $F_{X}^{-1}(u)$ 被定義為能使 $F(X)$ 到達 $u$ 最小的那個 $x$。
- 證明此方法成立：證明 $P(F_X^{-1}(U)\le x)=F_X(x)$

### Examples

- **Ex2**：生成 pdf 為 $f_X(x) = 3x^2,\ 0 < x < 1$ 的 random sample
- **Ex3**：生成服從 exponential distribution 的 random sample
- **Ex4**：生成服從 geometric distribution 的 random sample
  - geometric distribution 的 cdf：$F(x)=1-q^{x+1}$
  - 利用 $F_X(x_{i-1}) < u \le F_X(x_i)$ 的觀念來做
  - 可以得到 $1-q^x < u \le 1 - q^{x+1}$ 的式子
  - 推導後可以得到 $x + 1 = \lceil \log(1-u)/\log(q) \rceil$
  - R 實作方法：
    - `floor`：向下取整（公式不用再扣 1）
    - `ceiling`：向上取整（直接用 $x+1$ 的那個下去代）
- **Ex5**：生成服從 Poisson distribution 的 random sample
  - 利用 $f(x+1) = \frac{\lambda f(x)}{x+1}$；$F(x+1) = F(x) + f(x+1)$ 來做
  - 先定義一個 `for` 迴圈包含所有的 random sample
  - 再用一個 `while` loop 處理個別 sample 變換的狀況
- **Ex6**：Logarithmic Distribution（對數分佈）

---

## Acceptance-Rejection Method

- **使用時機**：當 inverse transformation 不好做的時候（cdf 反函數不好求）
- **Target density** $f(x)$：我們想要生成 random sample 的分佈
- **Instrumental density** $g(x)$：比較簡單的分佈，用其來幫助目標分布的生成
- **Concept**：
  - $X$ and $Y$ are random variables with pmf/pdf $f$ and $g$
  - Goal：讓這兩個分布滿足 $\frac{f(t)}{g(t)} \le c$ 這樣的條件

### Procedure

1. 生成服從 instrumental density $g$ 的 random sample $y$
2. 生成 uniform distribution 的 random sample $u$
3. 若 $u < \frac{f(y)}{cg(y)}$，則定 $x = y$；若不符合則回到第一步重新找一個

### A-R Method for Discrete r.v.

利用一個比較方便生成的 $Y$，進一步生成 $X$：

1. Generate $Y$, with corresponding $q_j$
2. Generate $u \sim Unif(0, 1)$
3. If $u < \frac{p_j}{cq_j}$, set $X = x_j$

### Examples

- **Ex1**：Generate triangular distribution
- **Ex2**：Generate gamma distribution

---

## Transformation Method

利用不同變數之間的變換，生成更多不同分布的變數。

### $Z \sim N(0, 1)$ 系列

- **Chi-Squared**：$V = Z^2 \sim \chi^2(1)$
- **F**：$F = \frac{U / m}{V / n}$ if $U \sim \chi^2(m)$ and $V \sim \chi^2(n)$ are independent
- **t(n)**：$T = \frac{Z}{\sqrt{V / n}}$ if $Z \sim N(0, 1)$ and $V \sim \chi^2(n)$ are independent
- **Box-Muller algorithm**：If $U, V \sim Unif(0, 1)$ and independent, then

  $$Z_1 = \sqrt{-2\log U}\,\cos(2 \pi V), \qquad Z_2 = \sqrt{-2\log U}\,\sin(2 \pi V)$$

  are independent standard normal variables.

### Gamma 系列

- If $U \sim \Gamma(r, \lambda)$ and $V \sim \Gamma(s, \lambda)$ are independent, then

  $$X = \frac{U}{U + V} \sim Beta(r, s)$$

- If $U, V \sim Unif(0, 1)$ are independent, then

  $$X = \left\lfloor 1 + \frac{\log(V)}{\log(1 - (1 - \theta)^U)} \right\rfloor \sim Logarithmic(\theta)$$

### Examples

- **Ex1**：生成 Beta distribution
- **Ex2**：生成 logarithmic distribution

---

## Convolution and Mixture

### Convolution $F_X^{*(n)}$

- $X_1, \dots, X_n$ is iid distributed
- $S = X_1 + \dots + X_n$ is called the **n-fold convolution** of $X$
- **Example：Chi-square distribution**
  1. 製作一個 $n \times v$ 的矩陣，矩陣內填滿了服從 $N(0, 1)$ 的 r.v.
  2. 對矩陣內所有 element 取平方，並計算每個 row 的加總
  3. 每個 row 的加總形成的 vector 即為我們想要的 chi-square random sample

### Mixture

- 一個分佈是很多個分佈的組合

  $$F_X(x) = \int_{-\infty}^{\infty} F_{X|Y=y}(x)f_Y(y) \, dy, \qquad \int_{-\infty}^{\infty}f_Y(y)\,dy = 1$$

- 如何在 R 生成一個 mixture distribution：
  1. 生成 $k \in \{1, 2\}$，$P(1) = P(2) = 0.5$
  2. 若 $k = 1$，生成一個 $x$ 服從 $f_{X_1}(x_1)$；反之則生成一個 $x$ 服從 $f_{X_2}(x_2)$
- **Examples**：
  - 生成兩個 Normal 的 mixture
  - 生成兩個 gamma 的 mixture
  - 生成多個 gamma 的 mixture
    - $F_X = \sum_{i=1}^5 \theta_i F_{X_i}$
    - $X_i \sim \Gamma(r = 3, \lambda_i = 1/i)$, $\theta_i = i/15$, $i = 1,2,\dots,5$

      | $i$ | $X_i$ | $\theta_i$ |
      | --- | --- | --- |
      | 1 | $\Gamma(r = 3, \lambda_1 = 1)$ | $1/15$ |
      | 2 | $\Gamma(r = 3, \lambda_2 = 1/2)$ | $2/15$ |
      | 3 | $\Gamma(r = 3, \lambda_3 = 1/3)$ | $3/15$ |
      | 4 | $\Gamma(r = 3, \lambda_4 = 1/4)$ | $4/15$ |
      | 5 | $\Gamma(r = 3, \lambda_5 = 1/5)$ | $5/15$ |

  - Poisson-Gamma mixture distribution

---

## Multivariate Normal Distribution

If $Z \sim N_d(0, I_d)$, then $CZ + b \sim N_d(b, CC^T)$.
Suppose $\Sigma$ can be factored so that $\Sigma = CC^T$ for some matrix $C$.
Then $CZ + \mu \sim N_d(\mu, \Sigma)$.

- **目標**：想辦法對 $\Sigma$ 做 factorization，進一步得到 $\Sigma^{1/2}$
- **方法**：
  1. Spectral decomposition → `eigen()`
  2. Choleski factorization → `chol()`
  3. Singular value decomposition (SVD) → `svd()`

### Algorithm for Generating Random Sample from $N_d(\mu, \Sigma)$

1. 生成一個 $n \times d$ 的矩陣 $Z$，矩陣中都是 $N(0, 1)$ 的隨機樣本
2. 對 $\Sigma$ 作 factorization，得到 $\Sigma = Q^TQ$
3. 對 $Z$ 進行轉換 $X = ZQ + 1_n\mu^T$
4. 轉換後 $X$ 的每個 row 即為服從 $N_d(\mu, \Sigma)$ 的一個樣本

### Spectral Decomposition Method

- $\Sigma = P \Lambda P^T$
  - $\Lambda$：對角矩陣，矩陣內的元素為 $\Sigma$ 的 eigenvalue
  - $P$：eigenvector matrix
- 利用 Spectral Decomposition Method 將 $\Sigma$ 拆成 $Q^TQ$

  $$\Sigma^{1/2} = P \Lambda^{1/2} P^T = Q \Rightarrow \Sigma = Q^TQ$$

### Singular Value Decomposition (SVD)

For any $m \times n$ matrix $X$:

- There exist orthogonal matrices $U$ of size $m \times m$ and $V$ of size $n \times n$ such that $X = UDV^T$
- $D$：a $m \times n$ matrix with nondiagonal entries all zero

$$\Rightarrow U = V = P \quad\text{and}\quad \Sigma^{1/2} = UD^{1/2}V^T$$

### Choleski Method

For a symmetric positive-definite matrix $X$, there exists an upper triangular matrix $Q$ such that $X = QQ^T$.

---

## Wishart Distribution

### General Method

先生成 $n$ 個服從 $N_d(\mu, \Sigma)$ 的 random sample，再利用 $M = X^TX$ 的轉換得到 Wishart 的 random sample。

1. $M = X^TX$, where $X$ is an $n \times d$ data matrix of a random sample from a $N_d(0, \Sigma)$ distribution
2. Then $M \sim W_d(\Sigma, n)$, a Wishart distribution with scale matrix $\Sigma$ and $n$ degrees of freedom
3. When $d = 1$, $X \sim N(0, \sigma^2)$, $W_1(\sigma^2, n) = \sigma^2\chi^2(n)$

### Bartlett's Decomposition

Let $A = (A_{ij})$ be the lower triangular $d \times d$ random matrix with independent entries satisfying

- $A_{ij} \stackrel{iid}{\sim} N(0, 1),\ i>j$
- $A_{ii} \sim \sqrt{\chi^2(n-i+1)}$, for $i = 1,\dots, d$

**作法**：

1. Obtain Choleski decomposition $\Sigma = LL^T$, where $L$ 是下三角矩陣且對角的值非負
2. Then $LAA^TL^T \sim W_d(\Sigma, n)$

<br>

# Lecture 2

## Simulating Survival Data

### Cox Proportional Hazards Model

Survival analysis 建立在 hazard function $h_i(t)$ 以及 survival function $S(t|X_i)$ 上。

- **Hazard function**：$h_i(t) = h_0(t)\exp(X_i\beta)$
  - 概念：在 $t$ 這一時間，瞬間發生事件的風險

    $$h(t) = \lim_{\Delta t \to \infty}\frac{P(t \le T \le t+\Delta t\,|\, T\ge t)}{\Delta t}$$

  - $h_0(t)$：baseline hazard
  - $X_i$：baseline covariate
  - $\beta$：vector of log hazard ratio
- **Survival function**：$S(t|X_i) = \exp(-H_0(t)\exp(X_i\beta))$
  - 概念：個體存活超過 $t$ 的機率 $S(t) = P(T > t)$
  - $H_0(t)$：cumulative baseline hazard function
  - Proportional Hazard：同一個風險函數下，survival rate 只受 $X_i$ 影響，和 $t$ 無關

### Simulate the Survival Data

Simulate survival times $T_i$ under a Cox Proportional Hazards model:

1. Set $S(t|X_i) = \exp(-H_0(t)\exp(X_i\beta)) = U_i$, where $U_i \sim Unif(0, 1)$
2. $H_0(t) = \frac{-\ln(U_i)}{\exp(X_i\beta)} \Rightarrow T_i = H_0^{-1}\left(\frac{-\ln(U_i)}{\exp(X_i\beta)}\right)$

Hazard function $h(t)$ 是自己訂的，$T_i$ 對於幾個比較常見的 hazard function 具有 closed form：

- **Exponential distribution**：constant baseline hazard

  $$h_0(t) = \lambda \Rightarrow H_0(t) = \lambda t$$

  $$T_i = \frac{-\ln(U_i)}{\lambda \exp(X_i \beta)}$$

- **Weibull distribution**：monotonic baseline hazard

  $$h_0(t) = \lambda \nu t^{\nu-1}\Rightarrow H_0(t) = \lambda t^\nu$$

  $$T_i = \left(\frac{-\ln(U_i)}{\lambda \exp(X_i \beta)}\right)^{1/\nu}$$

- **Gompertz distribution**：exponentially changing baseline hazard

  $$h_0(t) = a\exp(bt) \Rightarrow H_0(t) = \frac{a}{b}(\exp(bt) - 1)$$

  $$T_i = \frac{1}{b} \ln\left(1 - \frac{b\ln(U_i)}{\lambda \exp(X_i \beta)}\right)$$

若是遇到 hazard function 太複雜、沒有 closed form 的狀況
$\Rightarrow$ 利用 `simsurv` to solve the root-finding problem numerically

## Censoring Data

受試者在實驗進行到一半時退出的狀況。

For each individual $i$:

- $T_i$：true event time
- $C_i$：censoring time
- $d_i$：observed status indicator $d_i = I(T_i \le C_i)$
  - $d_i = 1$：Uncensored
  - $d_i = 0$：Censored
- Censor rate：

  $$P(T_i > C_i) = \frac{\sum_{i=1}^{n}I(d_i = 0)}{n} = 1-\bar{d}$$

### Simulation of Right-Censoring Mechanisms

**1. Independent Censoring Time Generation**

- 對每個資料 $i$ 個別模擬一個 censoring time $C_i$
- $C_i \sim Uniform(0, \tau)$
- 透過調整 upper bound $\tau$，可以控制資料中 censored data 的比例

**2. Constructing the Observed Variables**

- 每個 survival dataset 都會有兩個 observed variables $(Y_i, d_i)$
- Observed survival time $Y_i = \min(T_i, C_i)$
- Event indicator $d_i = I(T_i \le C_i)$
  - 0：Event was actually observed
  - 1：Individual was censored

**如何控制 censor rate**

- 透過調整 $\tau$ 來控制，$\tau$ 越大則 more censoring，反之亦然
- 找到最佳 $\tau$ 的方式：
  - Mathematical expectation
  - Iterative simulation / Root-finding

### Models for Survival Data

#### Cox Proportional Hazards Model (Cox PH)

**Model**

$$h(t|X_i) = h_0(t)\exp(X_i\beta)$$

- Covariate 的效果是**乘在 hazard 上**：不同個體的 hazard 比值固定，不隨時間 $t$ 改變
- **Hazard Ratio (HR)**：$X$ 增加一單位時，hazard 變為原本的 $\exp(\beta)$ 倍

  $$HR = \frac{h(t|X_i = 1)}{h(t|X_i = 0)} = \exp(\beta)$$

  - $\beta < 0$（$HR < 1$）：風險降低，例如治療有效
  - $\beta > 0$（$HR > 1$）：風險升高
- **Semi-parametric**：不需要指定 $h_0(t)$ 的形式，利用 partial likelihood 估計 $\beta$

  $$L(\beta) = \prod_{i:\, d_i = 1}\frac{\exp(X_i\beta)}{\sum_{j \in R(Y_i)}\exp(X_j\beta)}$$

  - $R(Y_i) = \{j : Y_j \ge Y_i\}$：在時間 $Y_i$ 仍在 risk set 中的個體
  - 只有發生事件（$d_i = 1$）的個體會貢獻分子，censored data 只出現在 risk set 中
- **對 survival curve 的影響**：$S(t|X_i) = S_0(t)^{\exp(X_i\beta)}$，covariate 讓曲線在 **y 軸方向**變形

**R 實作**

```r
library(survival)
cox_model <- coxph(Surv(time, status) ~ treatment, data = sim_df)
summary(cox_model)   # coef = beta, exp(coef) = HR
```

**Monte Carlo Power Analysis**

課堂設定：Weibull baseline $\lambda = 0.04,\ \nu = 1.4$，$\beta = -0.7$（$HR \approx 0.5$），censor rate 30%（$\tau = 22.07$）

1. 對每個樣本數 $N$，重複模擬 $n_{sim}$ 次：
   1. 生成 treatment $X_i \sim Bernoulli(0.5)$
   2. 用 inverse transform 生成 $T_i$，並生成 $C_i \sim Unif(0, \tau)$
   3. 得到 $(Y_i, d_i)$ 後 fit `coxph()`，記錄 $\beta$ 的 p-value
2. $\text{Empirical Power} = \frac{1}{n_{sim}}\sum I(\text{p-value} < 0.05)$
3. 結果：樣本數越大 power 越高；censor rate 越高，要達到同樣 power 所需的樣本越多

#### Accelerated Failure Time (AFT)

**Model**

Covariate 的效果是**直接作用在時間上**，讓事件發生的時間加速或減速：

$$T_i = \exp(-X_i\beta)\,T_{0i} \quad\Longleftrightarrow\quad \log T_i = -X_i\beta + \log T_{0i}$$

- $T_{0i}$：baseline（$X_i = 0$）下的存活時間
- **Survival function**：$S(t|X_i) = S_0\big(t\exp(X_i\beta)\big)$
- **Acceleration factor**：$\exp(X_i\beta)$
  - $\exp(X_i\beta) > 1$：時間被加速，事件較早發生
  - $\exp(X_i\beta) < 1$：時間被減速，存活時間拉長
- **對 survival curve 的影響**：covariate 讓曲線沿著 **x 軸方向**拉伸或壓縮

**Simulate AFT Data（Weibull baseline）**

$$S(t|X_i) = \exp\left(-\lambda\big(t\exp(X_i\beta)\big)^\nu\right) = U_i, \quad U_i \sim Unif(0, 1)$$

$$\Rightarrow T_i = \exp(-X_i\beta)\left(\frac{-\ln(U_i)}{\lambda}\right)^{1/\nu}$$

```r
T_true <- ((-log(U_t) / lambda)^(1 / nu)) / exp(X * beta)
```

**Weibull 同時是 PH 與 AFT**

上面的 Weibull AFT model 的 hazard 為

$$h(t|X_i) = \lambda\nu t^{\nu - 1}\exp(\nu X_i\beta)$$

仍然是 $h_0(t) \times$ 常數的形式，所以也滿足 PH 假設，對應的 $HR = \exp(\nu\beta)$。

**R 實作**

```r
aft_model <- survreg(Surv(time, status) ~ treatment, data = sim_df,
                     dist = "weibull")
summary(aft_model)
```

- `survreg()` 的模型為 $\log T = \mu + \gamma X + \sigma W$，所以 $\gamma = -\beta$
- `scale` $\sigma = 1/\nu$

**Cox PH vs. AFT**

| | Cox PH | AFT |
| --- | --- | --- |
| Covariate 作用在 | hazard | 時間 |
| 參數解釋 | Hazard Ratio $\exp(\beta)$ | Acceleration factor $\exp(X_i\beta)$ |
| Survival curve 變化 | y 軸方向變形 | x 軸方向拉伸 / 壓縮 |
| Baseline 分佈 | 不需指定（semi-parametric） | 通常需指定（parametric） |
| R 函數 | `coxph()` | `survreg()` |
