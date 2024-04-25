# HeisenTesting
Python program to numerically check testing domains in the Heisenberg group following arguments in [BBPT].

$${\red Hi Lior - this is Rami editting :)}$$


Consists of 3 Python files: Aux_0 ,  Aux_Check, main

_**Needed Python libraries:**_
(1) Numpy
(2) itertools
(3) math


**Theoretical description of code:**

Implements algorithm in upcoming paper to check whether a set is a testing domain for substitution data of the form $\big( \mathcal{A}, \lambda_0, S_0 \big)$ for dilation data $\Big( H_3(\mathbb{R}), d_H, (D_\lambda)_{\lambda>0}, H_3(2\mathbb{Z}), [-1,1)^3  \Big)$. We recall the meaning of the previous expressions also given in [BHP]. The action in $H_3(\mathbb{R})$ is given by

$$(x,y,z)\cdot (a,b,c):= \big( x+a, y+b, z+c +\frac{1}{2}(xb-ay) \big). $$

The metric $d_H$ is left invariant metric on $H_3(\mathbb{R})$ satisfying $d_H \big( (x,y,z), (a,b,c)  \big)= \Vert (x,y,z)^{-1} \cdot (a,b,c)  \Vert_{CK} $, with $\Vert \cdot \Vert_{CK}$ being the Cygan-Koranyi norm on $H_3(\mathbb{R})$ given by 

$$ \Vert (x,y,z) \Vert_{CK}:= \sqrt[4]{(x^2+y^2)^2+16z^2}. $$

For every $\lambda>0$, $D_\lambda: H_3(\mathbb{R}) \to H_3(\mathbb{R})$ is a dilation given by $D_\lambda(x,y,z)=\big( \lambda x, \lambda y, \lambda^2 z \big)$. $\lambda_0$ is some underlying stretch factor, for which $D_{\lambda_0}$ preserves the lattice $\Gamma:= H_3(2\mathbb{Z})$. $V=[-1,1)^3$ is a fundemantal testing domain for the lattice $\Gamma$. i.e., $H_3(\mathbb{R})= \underset{\gamma \in \Gamma}{\sqcup} \gamma V$.

Using these notions, the sets $V(n,M)$, for $n\in \mathbb{N}$ and $M\Subset \Gamma$, are defined recursively in [BHP]. $V(1,M):= D[M]\cdot D[V]$, and $V(1):=V(1,\{e\})$ is simply $D[V]$. $V(n+1,M)$ is given by $D[V(n,M)\cap \Gamma]\cdot V(1)$. Using these notations, the program checks for inputs lamb for underlying stretch factor $\lambda_0\in \mathbb{N}$, N for iteration number $N\in \mathbb{N}$, and two nonempty finite sets $T_1,T_2\subset \Gamma$ whether every $x\in D^N[V]\cap \Gamma$ there exists a $\gamma_x \in \Gamma$ such that

$$ xT_1 \subseteq D^N(\gamma_x) V(N,T_2). $$

This is implemented in the function InclCheck(lamb, N, T_2, T_1) which returns True if the condition holds for all $x\in D^N[V]\cap \Gamma$ and False else.

If the last condition holds when  $T_1$ is a testing domain, then we deduce that $T_2$ is also a testing domain. It follows, that if $T_1,T_2,...,T_k \subseteq \Gamma$ are finite nonempty subsets such that InclCheck(lamb, N_j, T_j, T_(j+1)) returns True for all $j$ when $T_1$ is known to be a testing domain, then $T_k$ is also a testing domain. For that reason, the functions

Since $D^N[V]\cap \Gamma$ is of the order of $\lambda_0^{4N}$, we can see that the program has at least exponential runtime with respect to $N$. For this reason, we prefer to work with a sequence of checks of the program, rather than a direct application of the program InclCheck(lamb, N, T_1, T_k) for very large $N$.

To minimize checks and since all sets for which we check set inclusions are a unions of intervals with respect to the 'z' coordinate, these sets are saved as  two connsecutive triples corresponding to the edges of the interval. They are considered as the 'z'-faces of the sets. The set inclusions are checked by whether by inequalities of the 'z'-intervals.

