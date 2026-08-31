# Structure Function

## 1. Time Series and Light Curves

### 1.1 Time series

A *time series* describes the evolution of a quantity with time. In a probabilistic description, the value of the quantity at a given time is treated as a random variable. For a finite set of time points, a time series can be represented as
$$
\mathcal{X} = \{(t_i,X_i)\}_{i=1}^{n},
$$
where $t_i$ is the time **(treated as a fixed quantity)** associated with the $i$-th observation and $X_i$ is the random variable describing the value of the time series at $t_i$. More generally, the sequence can be regarded as a discrete sampling of a stochastic process $X(t)$,
$$
X_i \equiv X(t_i).
$$

For a particular realization of the time series, the random variables $X_i$ take numerical values, which are denoted by lowercase symbols, $x_i \sim X_i$. Thus, $X_i$ and $x_i$ represent different mathematical objects: $X_i$ is a random variable, whereas $x_i$ is one realization of that random variable.

### 1.2 Astronomical light curves

An astronomical light curve is a time series describing the temporal variation of an observed astronomical quantity, such as flux, magnitude, or count rate. An idealized light curve can therefore be regarded as a realization of a stochastic process. If the process is sampled at times $t_1,\ldots,t_n$, the corresponding realization is
$$
\{(t_i,x_i)\}_{i=1}^{n}.
$$

In real observations, **the value recorded by an instrument is affected by measurement uncertainty**. It is therefore useful to distinguish the underlying signal from the measurement process. A simple measurement model is
$$
Y_i = X_i+\epsilon_i,
$$
where
- $X_i$ is the underlying random variable;
- $\epsilon_i$ is the measurement error;
- $Y_i$ is the random variable describing the measured quantity.

For a particular observation, $Y_i$ takes a realized value $y_i$. The actual dataset is therefore
$$
\mathcal{D} = \{(t_i,y_i,\sigma_i)\}_{i=1}^{n},
$$
where $\sigma_i$ characterizes the uncertainty associated with the $i$-th measurement.

Under the commonly adopted zero-mean error model,
$$
E[\epsilon_i]=0,
\qquad
\operatorname{Var}(\epsilon_i)=\sigma_i^2.
$$
The observed value $y_i$ should therefore not be identified with the underlying random variable $X_i$. Rather, it is a realization of the measurement random variable $Y_i$. This distinction can be summarized schematically as
$$
X_i
\longrightarrow
Y_i
\longrightarrow
y_i,
$$
where $X_i$ describes the underlying signal, $Y_i$ describes the noisy measurement, and **$y_i$ is the numerical value actually recorded**.

Light curves are generally irregularly sampled. In particular, $t_{i+1}-t_i$ need not be constant. Consequently, consecutive observations do not necessarily probe the same physical time scale. For example, two consecutive observations may be separated by days in one part of a light curve and by months in another. A variability statistic based only on consecutive observations would therefore mix different time scales. The *structure function (SF)* avoids this problem by considering pairs of observations and organizing them according to their time separation (or lag).

## 2. Observation Pairs and Time Lags

### 2.1 Observation pairs

Consider two random variables in a time series,
$$
X_i=X(t_i),\qquad X_j=X(t_j).
$$
Their difference is $\Delta X_{ij} = X_j - X_i$; the time separation of the pair is $\Delta t_{ij} = |t_j-t_i|$.
Each pair therefore provides a time separation and a corresponding change in the time-series value: $(\Delta t_{ij},\Delta x_{ij})$.

For a light curve containing $n$ observations, the number of distinct unordered pairs is
$$
N_{\rm pair} = \binom{n}{2}
= \frac{n(n-1)}{2}.
$$

> ‼️ The number of pairs should not be confused with the number of statistically independent measurements. Different pairs may share one or both observations, and therefore their pairwise statistics are **generally correlated**.

### 2.2 Pair is the basic unit

The purpose of an SF is to characterize variability as a function of time scale.
For two observations separated by a lag $\Delta t_{ij}$, the difference $\Delta x_{ij}$ provides one realization of the change occurring over that time scale. A single pair is not sufficient to characterize the statistical variability at that lag, but many pairs with similar lags can collectively provide an estimate.

Thus, the fundamental computational unit of an SF is an **observation pair**, rather than an individual observation or a consecutive pair of observations.
This leads to the basic structure
$$
\text{light curve}
\rightarrow
\text{observation pairs}
\rightarrow
\text{time lags}
\rightarrow
\text{pairwise statistics}.
$$
SF is then obtained by statistically summarizing those pairwise quantities as a function of lag.

### 2.3 Lag bins

A finite light curve provides only a finite set of pairwise lags. In general, there are no pairs at every possible value of $\tau$. Pairs are therefore grouped into **lag bins**. Let the $k$-th lag bin be
$$
B_k = [\tau_{k,\min},\tau_{k,\max}).
$$
An observation pair $(i,j)$ belongs to this bin when $\Delta t_{ij}\in B_k$. Define the set of pairs in the bin as
$$
P_k
=
\left\{
(i,j):
i<j,\;
\Delta t_{ij}\in B_k
\right\},
$$
and let $N_k=|P_k|$ be the number of pairs in the bin.

The bin is associated with a representative lag $\tau_k$. The exact definition of $\tau_k$ depends on the adopted binning scheme. Common choices include the arithmetic mean or geometric mean of the pairwise lags.

The use of a lag bin implicitly assumes that SF does not change substantially across the width of the bin. In other words, the pairs within a bin are treated as probing approximately the same time scale.

## 3. Structure Function

### 3.1 Population-level definition

The second-order SF describes the expected squared change of a stochastic process as a function of time lag. For a stationary process $X(t)$, it is defined as
$$
\boxed{
SF^2(\tau) = E\left[
\left( X(t+\tau)-X(t) \right)^2
\right], \qquad SF(\tau) = \sqrt{SF^2(\tau)}.
}
$$
Here, $SF^2(\tau)$ is a statistical quantity defined through an expectation over the underlying stochastic process. The second-order SF specifically characterizes the second moment of the pairwise difference.

The SF therefore answers the question: *How large is the typical change in the underlying signal when two points are separated by a time lag $\tau$?*

### 3.2 Stationarity

The definition above can be applied to a stationary process, for which the statistical properties of the process are invariant under a shift in time. For a stationary process, $\Delta X(t, \tau) \equiv X(t+\tau)-X(t)$ **depends on $\tau$ but not on the absolute time $t$**. Therefore, we can write $\Delta X(t, \tau)$ as $\Delta X(\tau)$. This allows the pairwise differences from different locations in the light curve to be treated as samples of the same statistical quantity at a given lag.

**Stationarity is therefore an important conceptual assumption behind the usual SF interpretation.** In practice, an astronomical light curve may only be approximately stationary over the time interval being analyzed.

### 3.3 Relation to the autocovariance function

For a second-order stationary process with *autocovariance function*
$$C(\tau)=\operatorname{Cov}[X(t),X(t+\tau)],$$
the second-order SF satisfies
$$
SF^2(\tau)=2\left[C(0)-C(\tau)\right].
$$
If the process has variance $\sigma_X^2=C(0)$, then
$$
\boxed{
SF^2(\tau) = 2\sigma_X^2 \left[1-\rho(\tau)\right],
}
$$
where $\rho(\tau)=\frac{C(\tau)}{C(0)}$ is the *autocorrelation function*.

Thus, SF and autocorrelation function contain closely related information about the temporal dependence of a stationary process.

### 3.4 From the statistical quantity to an estimator

The definition of $SF(\tau)$ contains an expectation and therefore refers to the underlying stochastic process. A real light curve, however, provides only a finite realization of that process.

The expectation must consequently be replaced by an estimator constructed from the observed data. This distinction is fundamental:
- $X_i$ is a random variable;
- $x_i$ is a realization of $X_i$;
- $\Delta X_{ij}$ is a random variable describing the difference between two process values;
- $\Delta x_{ij}$ is the corresponding realized difference;
- $SF(\tau)$ is a population-level statistical quantity;
- $\widehat{SF}(\tau)$ is an estimator calculated from finite data.

SF is therefore not itself "the average of the observed differences." Rather, an average of suitable observed pairwise statistics is one possible estimator of the underlying SF.

## 4. Estimating the SF from a Finite Light Curve

### 4.1 Pairwise statistics

For each observation pair, define a pairwise statistic

$$
X_p
=
h(x_i,x_j,\sigma_i,\sigma_j),
$$

where $p$ indexes the observation pair and $h$ specifies how that pair contributes to the SF estimator.

For a given lag bin, the collection of pairwise statistics is then summarized by an aggregation functional $g$,

$$
\widehat{SF}(\tau_k)
=
g\left(
\{X_p:p\in P_k\}
\right).
$$

This separates two concepts:

1. the **pairwise statistic**, which determines what information is extracted from each pair;
2. the **aggregation rule**, which determines how the pairwise values are combined within a lag bin.

This framework allows different SF estimators to be expressed within the same computational structure.

### 4.2 Second-order SF

The standard second-order estimator uses the squared pairwise difference,

$$
X_p
=
(\Delta x_p)^2.
$$

The estimator of the second-order structure function is then

$$
\boxed{
\widehat{SF^2}(\tau_k)
=
\frac{1}{N_k}
\sum_{p\in P_k}
(\Delta x_p)^2
}
$$

and

$$
\boxed{
\widehat{SF}(\tau_k)
=
\sqrt{
\widehat{SF^2}(\tau_k)
}.
}
$$

The corresponding computational sequence is

$$
\{(t_i,x_i)\}
\rightarrow
\{\Delta t_p,\Delta x_p\}
\rightarrow
\{(\Delta x_p)^2\}_{p\in P_k}
\rightarrow
\widehat{SF^2}(\tau_k).
$$

This is the default second-order SF estimator.

### 4.3 Noise-corrected second-order SF

The observed pairwise differences contain both intrinsic variability and measurement noise.

Under the measurement model

$$
Y_i=X_i+\epsilon_i,
$$

the difference between two measured random variables is

$$
Y_i-Y_j
=
(X_i-X_j)
+
(\epsilon_i-\epsilon_j).
$$

Suppose that the measurement errors satisfy

$$
E[\epsilon_i]=0,
$$

and are independent of the underlying signal and of each other, with

$$
\operatorname{Var}(\epsilon_i)=\sigma_i^2.
$$

Then

$$
E\left[
(Y_i-Y_j)^2
\right]
=
E\left[
(X_i-X_j)^2
\right]
+
\sigma_i^2+\sigma_j^2.
$$

The observed second moment therefore contains an additive contribution from measurement noise.

This motivates the noise-corrected pairwise statistic

$$
\boxed{
X_p
=
(\Delta y_p)^2
-
(\sigma_i^2+\sigma_j^2).
}
$$

A corresponding estimator is

$$
\boxed{
\widehat{SF^2}(\tau_k)
=
\frac{1}{N_k}
\sum_{p\in P_k}
\left[
(\Delta y_p)^2
-
(\sigma_i^2+\sigma_j^2)
\right].
}
$$

Under the assumptions above, this estimator is unbiased for the intrinsic second moment.

The assumptions are important. In particular, the correction relies on:

1. zero-mean measurement errors;
2. measurement errors independent of the underlying signal;
3. independent errors between observations;
4. correctly specified measurement variances $\sigma_i^2$.

If measurement errors are correlated, the noise contribution is not simply $\sigma_i^2+\sigma_j^2$. For example, if

$$
\operatorname{Cov}(\epsilon_i,\epsilon_j)\neq0,
$$

then

$$
\operatorname{Var}(\epsilon_i-\epsilon_j)
=
\sigma_i^2+\sigma_j^2
-
2\operatorname{Cov}(\epsilon_i,\epsilon_j).
$$

Therefore, the simple noise correction can be biased when correlated measurement errors are present.

### 4.4 Negative noise-corrected estimates

Because the noise-corrected estimator subtracts an estimate of the measurement-noise contribution, its finite-sample value can be negative:

$$
\widehat{SF^2}(\tau_k)<0.
$$

This does not imply that the population quantity $SF^2(\tau)$ is negative. Rather, it indicates that the estimated intrinsic second-order variability in that bin is consistent with being dominated by measurement noise, or that the finite-sample estimator has fluctuated below zero.

Since

$$
SF(\tau)=\sqrt{SF^2(\tau)},
$$

a negative estimate of $SF^2$ does not have a real-valued square root. AGNSF therefore represents the corresponding SF value as undefined (NaN) rather than taking an absolute value or otherwise forcing the result to be positive.

## 5. Mean-Absolute-Difference SF

The second-order SF is not the only possible way to summarize pairwise differences.

An alternative is to use the mean absolute difference,

$$
X_p
=
|\Delta x_p|.
$$

The mean absolute difference within a lag bin is

$$
\left\langle
|\Delta x|
\right\rangle
=
\frac{1}{N_k}
\sum_{p\in P_k}
|\Delta x_p|.
$$

AGNSF adopts the Gaussian-normalized convention

$$
\boxed{
SF^2
=
\frac{\pi}{2}
\left\langle
|\Delta x|
\right\rangle^2.
}
$$

The factor $\pi/2$ follows from the relation

$$
E|Z|
=
\sqrt{\frac{2}{\pi}}\,\sigma_Z
$$

for a zero-mean Gaussian random variable

$$
Z\sim N(0,\sigma_Z^2).
$$

Thus, under a Gaussian assumption for the pairwise differences,

$$
\sigma_Z^2
=
\frac{\pi}{2}
\left(E|Z|\right)^2.
$$

The mean-absolute-difference estimator is consequently a Gaussian-equivalent estimate of the second-order variance scale.

A noise-corrected version can be constructed by subtracting the average measurement-noise variance contribution,

$$
\widehat{SF^2}
=
\frac{\pi}{2}
\left\langle
|\Delta y|
\right\rangle^2
-
\left\langle
\sigma_i^2+\sigma_j^2
\right\rangle,
$$

where the averages are taken over the pairs in the lag bin.

The mean-absolute-difference estimator should not be confused with a median-based estimator. In particular, replacing the mean with a median changes the statistical estimator and its uncertainty properties.