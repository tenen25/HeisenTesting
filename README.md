# HeisenTesting
Python program to numerically check testing domains in the Heisenberg group described in [BBPT], _Approximations of symbolic substitution systems_.

_**Needed Python libraries:**_
(1) Numpy
(2) itertools
(3) math


_**Parameters needed from user:**_

Parameters needed from the user are fed in the main function. They are:

(1) Underlying stretch factor - 'lamb'.
(2) Iteration number - 'N'.
(3) Set which is known to be a testing domain - 'known_set'.
(4) Set which is suspected to be a testing domain - 'init_set'.

The sets are input as lists of triples describing 'z'-intervals. The even or odd entries in the list describe the bottom or upper part of the 'z'-interval accordingly. For example the set $\\{-2,0\\}^2\times \\{-6,-4,-2,0,2\\}$ can be saved as

// 
[ (-2,-2,-6), (-2,-2,2), (-2,0,-6), (-2,0,2), (0,-2,-6), (0,-2,2), (0,0,-6), (0,0,2)  ]
 //

_**Preliminary computations:**_

The following sets are saved to the memory prior to the final check:

(1) $K = V(1)\cap H_3(2\mathbb{Z})$, 
(2) $K3 = V(2)\cap H_3(2\mathbb{Z})$ for the case when $\lambda_0=3$, 
(3) $T_1 = \\{ -2,0\\}^2\times \\{-6,-4,...,4,6\\}$, 
(4) $T_2 =\\{ 0 , -2 \\}^2 \times \\{ -4,-2,0,2 \\}$.

**Theoretical description of code:**

The program implements Algorithm 1 in [BBPT] to check whether a set is a testing domain for substitution data of the form $\big( \mathcal{A}, \lambda_0, S_0 \big)$ for dilation data $\Big( H_3(\mathbb{R}), d_H, (D_\lambda)_{\lambda>0}, H_3(2\mathbb{Z}), [-1,1)^3  \Big)$. **We note that the user can not change the input dilation datum, as we are only dealing with the same dilation datum in the Heisenberg group.** 
We attach a picture of Algorithm 1 here for the sake of convenience:

![Algorithm 1](https://github.com/tenen25/HeisenTesting/assets/75997072/4ffe97ad-bd8c-41ed-a627-912194068b05)

We recall the meaning of the previous expressions. The action in $H_3(\mathbb{R})$ is given by


$$(x,y,z)\cdot (a,b,c):= \big( x+a, y+b, z+c +\frac{1}{2}(xb-ay) \big). $$

The metric $d_H$ is left invariant metric on $H_3(\mathbb{R})$ satisfying $d_H \big( (x,y,z), (a,b,c)  \big)= \Vert (x,y,z)^{-1} \cdot (a,b,c)  \Vert_{CK} $, with $\Vert \cdot \Vert_{CK}$ being the Cygan-Koranyi norm on $H_3(\mathbb{R})$ given by 

$$ \Vert (x,y,z) \Vert_{CK}:= \sqrt[4]{(x^2+y^2)^2+z^2}. $$

For every $\lambda>0$, $D_\lambda: H_3(\mathbb{R}) \to H_3(\mathbb{R})$ is a dilation given by $D_\lambda(x,y,z)=\big( \lambda x, \lambda y, \lambda^2 z \big)$. $\lambda_0$ is some underlying stretch factor, for which $D_{\lambda_0}$ preserves the lattice $\Gamma:= H_3(2\mathbb{Z})$. $V=[-1,1)^3$ is a fundemantal testing domain for the lattice $\Gamma$. i.e., $H_3(\mathbb{R})= \underset{\gamma \in \Gamma}{\sqcup} \gamma V$. We denote $D:=D_{\lambda_0}$ for brevity. 

Using these notions, the sets $V(n,M)$, for $n\in \mathbb{N}$ and $M\Subset \Gamma$, are defined recursively in [BHP], _Symbolic substitution systems beyond abelian groups_. $V(1,M):= D[M]\cdot D[V]$, and $V(1):=V(1,\{e\})$ is simply $D[V]$. $V(n+1,M)$ is given by $D[V(n,M)\cap \Gamma]\cdot V(1)$. Using these notations, the program checks for inputs lamb for underlying stretch factor $\lambda_0\in \mathbb{N}$, N for iteration number $N\in \mathbb{N}$, and two nonempty finite sets $T_1,T_2\subset \Gamma$ whether every $x\in D^N[V]\cap \Gamma$ there exists a $\gamma_x \in \Gamma$ such that

$$ xT_1 \subseteq D^N(\gamma_x) V(N,T_2). $$

To this end we define three logical formulas:

(1) The formula $I(x, \gamma, \lambda_0, N, T_2, T_1)$ is defined as

$$ xT_1 \subseteq D^N(\gamma) V(N,T_2), $$

which has $x$ and $\gamma$ as free variables. A version of this is implemented in line 25 of Aux_Check.py .

(2) The formula $\phi(x, \lambda_0, N, T_2, T_1)$ is defined as

$$ \exists \gamma\in \Gamma \enspace I(x, \gamma, \lambda_0, N, T_2, T_1) , $$

which has $x$ as a free variable. A version of this formula is implemented by the functions Check_x(x, lamb, N, shifted, fixed) .

(3) The formula $\psi(x, \lambda_0, T_2, T_1 )$ is defined by

$$ \forall x\in D^N[V]\cap \Gamma \enspace \phi(x, \lambda_0, N, T_2, T_1)=  \forall x\in D^N[V]\cap \Gamma \enspace \exists \gamma\in \Gamma \enspace  I(x, \gamma, \lambda_0, N, T_2, T_1). $$

This is implemented in the function InclCheck(lamb, N, T_2, T_1) in Aux_Check.py.

If the last condition holds when  $T_1$ is a testing domain, then we deduce that $T_2$ is also a testing domain. It follows, that if $T_1,T_2,...,T_k \subseteq \Gamma$ are finite nonempty subsets such that InclCheck(lamb, N_j, T_j, T_(j+1)) returns True for all $j$ when $T_1$ is known to be a testing domain, then $T_k$ is also a testing domain. For that reason, the functions

Since $D^N[V]\cap \Gamma$ is of the order of $\lambda_0^{4N}$, being the size of $D^N[V]\cap \Gamma$, we can see that the program has at least exponential runtime with respect to $N$. For this reason, we prefer to work with a sequence of checks of the program, rather than a direct application of the program InclCheck(lamb, N, T_1, T_k) for very large $N$.

To minimize checks and since all sets for which we check set inclusions are a unions of intervals with respect to the 'z' coordinate, these sets are saved as  two connsecutive triples corresponding to the edges of the interval. They are considered as the 'z'-faces of the sets. The set inclusions are checked by whether by inequalities of the 'z'-intervals.


_**Current results from the code:**_

Recall the sets we defined earlier: 

(1) $K = V(1)\cap H_3(2\mathbb{Z})$, 
(2) $K3 = V(2)\cap H_3(2\mathbb{Z})$ for the case when $\lambda_0=3$, 
(3) $T_1 = \\{ -2,0\\}^2\times \\{-6,-4,...,4,6\\}$, 
(4) $T_2 =\\{ 0 , -2 \\}^2 \times \\{ -4,-2,0,2 \\}$.

For $\lambda_0=4$, the tuple $\big( T_2 ,1\big)$ is a testing tuple. This follows since InclCheck(4, 1, K, T_1) and InclCheck(4, 2, T_1, T_2) both return True.

For $\lambda_0=3$, the tuple $\big( T_1 ,1\big)$ is a testing tuple. This follows since InclCheck(3, 2, K3, K) and InclCheck(3, 2, K, T_1) both return True.

The corresponing gammas to x's are given in the files ~~FILL~~.

_**References:**_

[BBPT] R. Band, S. Beckus, F. Pogorzelski, L. Tenenbaum. Approximations of symbolic substitution systems. 

[BHP] S. Beckus, T. Hartnick, and F. Pogorzelski. Symbolic substitution systems beyond abelian groups. https://arxiv.org/abs/2109.15210



