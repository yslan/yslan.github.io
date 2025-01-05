---
title: Gauss-Lobatto-Legendre Points
date: 2025-01-05 13:00:00 -0600
categories: [Numerics, Numerical integral]
tags: [numerical integral,sem]     # TAG names should always be lowercase
toc: true
math: true
---

This is a brief note on the choice of Gauss-Lobatto-Legendre (GLL) points and
a summary of degree counting. GLL points are widely used in spectral element methods
(SEMs), but initially, I only relied on their definition and did not fully connect
them to the Gauss nodes, which are standard textbook material.
It wasn't until I TAed CS450 that I ~~was forced to~~ finally had the chance to
sit down and think it through.


## Legendre Polynomials

With the unweighted inner product defined as $(f,g)=\int_{-1}^1 f(x)g(x)dx$,
the Legendre polynomials, denoted as $L_N(x)$ with degree $N$, can be constructed
via the orthogonalization of the monomials basis $$\{1,x,x^2,\cdots,x^N\}$$.
Therefore, we have:

$$\int_{-1}^1L_N(x)\ x^j\ dx = 0, \quad\quad\forall j=0,\cdots,N-1$$

It can be shown (homework 9, CS450, FA24) that all Legendre polynomials have simple
root within the open interval $(-1,1)$. Using these roots, denoted as
$$\{\xi^{GL}_i\}_{i=1,\cdots,N}$$, as the nodes, we can define the Lagrange
polynomials $\ell^{GL}_i(x)$ and construct the Gauss-Legendre (GL) quadrature:

$$\int_{-1}^1 f(x)\ dx \approx \sum_{i=1}^N w^{GL}_i f(\xi^{GL}_i)\quad
{\rm where\ } w^{GL}_i:= \int_{-1}^1\ell^{GL}_i(x)\ dx,$$

Also, it can be shown (Part B, hw 9) that the quadrature rule is exact for any
polynomial $f$ of degree up to $2N-1$, using the polynomial remainder theorem
and the orthogonality relation. 


## Gauss-Lobatto

To include the endpoints, which is necessary to represent the boundary conditions
for SEM, the nodes must contain $\pm 1$. Therefore, the $n$ collocation points
are the roots of $$p_n(x) = (x^2-1)q_{n-2}(x)$$ where $$q_{n-2}\in\mathbb{P}_{n-2}$$ is some
polynomial of degree $n-2$. Here we examine two possible choices for $q$:
$q_{n-2}=L_{n-2}$ or $q_{n-2}=L'_{n-1}$.

- Case $q_{n-2}=L_{n-2}$

  Due to the extra order introduced by $x^2-1$, the orthogonality relations

  $$\int_{-1}^1 p_n(x)\ x^j\ dx = \int_{-1}^1 L_{n-2}(x)\ (x^{j+2}-x^j)\ dx = 0$$

  hold only for $j+2\leq n-3$, or $j=0,\cdots,n-5$.
  

- Case $q_{n-2}=L'_{n-1}$

  By applying the integration by part, we obtain

  $$\begin{align*}
  0 &= \int_{-1}^1 p_n(x)\ x^j\ dx = \int_{-1}^1 (x^2-1) x^j L'_{n-1}(x)\ dx
  = \int_{-1}^1 (x^2-1) x^j \ d L_{n-1}(x)\\
  &= \cancelto{0}{\left[(x^2-1) x^j L_{n-1}(x)\right]\bigg|_{-1}^1}
  - \int_{-1}^1 L_{n-1}(x) \left((j+2)x^{j+1}-jx^{j-1}\right)\ dx.
  \end{align*}$$

  This relaxes two extra orders and the equation holds when $j+1\leq n-2$ or
  $j=0,\cdots, n-3$.

The second choice leverages more orthogonality provided by the Legendre polynomials.
The Gauss-Lobatto-Legendre (GLL) nodes $$\{\xi^{GLL}_i\}_{i=1,\cdots,n}$$ are then
defined as the roots of $$(x^2-1)L'_{n-1}(x)$$ and the quadrature rule is exact
for polynomials up to degree $2n-3$.


## Comparison

In summary, for degree $N$ Legendre polynomials $L_N(x)$ and $n$ collocation points,
the following table summarizes the quadrature rule of Legendre family.
Here, the Gauss-Radau-Legendre (GRL) points are for the semi-open interval that
includes $-1$, which is convenient to avoid singularity of the radial axis under
the polar coordinate.

| Method | Roots of | Range | Orthogonality | Quad. Exact Deg. |
|:---:|:---:|:---:|:---:|:---:|
| Gauss          | $L_N(x)$           | $(-1,1)$ | $n-1$ | $2n-1$        |
| Gauss-Radau    | $L_N(x)+L_{N-1}(x)$| $[-1,1)$ | $n-2$ | $2n-2$        |
| Gauss-Lobatto  | $(x^2-1)L'_N(x)$   | $[-1,1]$ | $n-3$ | $2n-3 = 2N-1$ |


