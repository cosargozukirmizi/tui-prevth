# tui-prevth
tui for probabilistic evolution theory solution of ordinary differential equations with multinomial right-hand-sides


This program solves the initial value problem of ordinary differential equations with multinomial right-hand-side functions. Four different initial value problems are put in the source code. The code can be manipulated to accomodate other initial value problems with multinomial right-hand-sides. 

![tui-prevth](Screenshot_2023-08-22_12-21-22.png?raw=true "tui-prevth")

The calculations are performed using rational arithmetic. In this way, error accumulation because of finite nature of double arithmetic is avoided. 

The command line interface has five sets of input. 

ODE preset is for choosing the initial value problem. The ODEs are van der Pol, quartic anharmonic oscillator, Henon-Heiles and Rabinovich-Fabrikant. All initial values for all functions are taken as 1/2 (This can be changed easily within the source code, but currently not from the tui).

The ODEs are given below.

## vdPol

The van der Pol ODE is
```math
\begin{eqnarray}
 \dot{x} &=& \mu x - \frac{1}{3}\mu x^{3} - \mu y \\
 \dot{y} &=& \frac{1}{\mu} x 
\end{eqnarray}
```
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2})}(t) \equiv x(t)^{\ell_{1}}y(t)^{\ell_{2}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& \mu \, u^{(0,0)} u^{(1,0)} 
  - \frac{1}{3}\mu \, u^{(1,0)} u^{(2,0)} 
  - \mu \, u^{(0,0)} u^{(0,1)}  \\
  \dot{u}^{(0,1)} &=& \frac{1}{\mu}\, u^{(0,0)} u^{(1,0)} 
\end{eqnarray}
```
using the heuristic H1 used in beam search for pure quadratization. We will take $\mu$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(1,0)} 
  - 1/3\, u^{(1,0)} u^{(2,0)} 
  - 1\, u^{(0,0)} u^{(0,1)}
  \\
  \dot{u}^{(0,1)} &=& 1\, u^{(0,0)} u^{(1,0)}\\
  \dot{u}^{(0,0)} &=& 0\\
  \dot{u}^{(2,0)} &=& 2\, u^{(1,0)} u^{(1,0)} 
  - 2/3\, u^{(2,0)} u^{(2,0)} 
  - 2\, u^{(0,1)} u^{(1,0)}
\end{eqnarray}
```

## QuartAnhOsc

The Quartic anharmonic oscillator ODE is
```math
\begin{eqnarray}
 \dot{q} &=& \frac{1}{\mu} p \\
 \dot{p} &=& - k_{1} q - k_{2} q^{3}\\
\end{eqnarray}
``` 
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2})}(t) \equiv q(t)^{\ell_{1}}p(t)^{\ell_{2}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& \frac{1}{\mu} \, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& -k_{1}\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& k_{2}\, u^{(1,0)} u^{(2,0)}\\
\end{eqnarray}
```
using the heuristic H1 used in beam search for pure quadratization. We will take $k_{1}$, $k_{2}$ and $\mu$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& -1\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& 1\, u^{(1,0)} u^{(2,0)}\\
  \dot{u}^{(0,0)} &=& 0\\
  \dot{u}^{(2,0)} &=& 2\, u^{(0,1)} u^{(1,0)}
\end{eqnarray}
```

## HenHeil

The Henon-Heiles ODE is
```math
\begin{eqnarray}
  \dot{x} &=& p_{x}\\
  \dot{p}_{x} &=& -x \nonumber \\
 &-& 2\, \lambda xy\\
  \dot{y} &=& p_{y}\\
  \dot{p}_{y} &=& -y \nonumber \\
 &-& \lambda x^{2} \nonumber \\
 &+& \lambda y^{2} 
\end{eqnarray}
```
where there are four sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2},\ell_{3},\ell_{4})}(t) 
 \equiv x(t)^{\ell_{1}}p_{x}(t)^{\ell_{2}}y(t)^{\ell_{3}}
p_{y}(t)^{\ell_{4}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,1,0,0)}\\
  \dot{u}^{(0,1,0,0)} &=& -1\, u^{(0,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &-& 2\lambda \, u^{(0,0,1,0)} u^{(1,0,0,0)}\\
  \dot{u}^{(0,0,1,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,0,0,1)}\\
  \dot{u}^{(0,0,0,1)} &=& -1\, u^{(0,0,0,0)} u^{(0,0,1,0)} \nonumber \\
 &-& \lambda \, u^{(1,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &+& \lambda \, u^{(0,0,1,0)} u^{(0,0,1,0)}
\end{eqnarray}
``` 
using the heuristic H1 used in beam search for pure quadratization. We will take $\lambda$ as 1. After pure quadratization, the ODE set is 
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,1,0,0)}\\
  \dot{u}^{(0,1,0,0)} &=& -1\, u^{(0,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &-& 2\, u^{(0,0,1,0)} u^{(1,0,0,0)}\\
  \dot{u}^{(0,0,1,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,0,0,1)}\\
  \dot{u}^{(0,0,0,1)} &=& -1\, u^{(0,0,0,0)} u^{(0,0,1,0)} \nonumber \\
 &-& 1\, u^{(1,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &+& 1\, u^{(0,0,1,0)} u^{(0,0,1,0)}\\
  \dot{u}^{(0,0,0,0)} &=& 0
\end{eqnarray}
```

## RabFab

The Rabinovich-Fabrikant ODE is 
```math
\begin{eqnarray}
 \dot{x} &=& yz - y + x^{2}y + \gamma x \\
 \dot{y} &=& 3xz + x - x^{3} + \gamma y \\
 \dot{z} &=& -2\alpha z - 2 xyz 
\end{eqnarray}
```
where there are two sought functions. 
The definition
```math
\begin{equation}
 u^{(\ell_{1},\ell_{2},\ell_{3})}(t) 
 \equiv x(t)^{\ell_{1}}y(t)^{\ell_{2}}z(t)^{\ell_{3}} 
\end{equation}
```
can be used to rewrite in the form
```math
\begin{eqnarray}
\dot{u}^{(1,0,0)} &=& 1\, u^{(0,0,1)} u^{(0,1,0)} 
- 1\, u^{(0,0,0)} u^{(0,1,0)} \nonumber\\
&+& 1\, u^{(1,0,0)} u^{(1,1,0)} 
+ \gamma \, u^{(0,0,0)} u^{(1,0,0)}\\
\dot{u}^{(0,1,0)} &=& 3\, u^{(0,0,1)} u^{(1,0,0)} 
+ 1\, u^{(0,0,0)} u^{(1,0,0)} \nonumber\\
&-& 1\, u^{(1,0,0)} u^{(2,0,0)} 
+ \gamma \, u^{(0,0,0)} u^{(0,1,0)} \\
\dot{u}^{(0,0,1)} &=& -2\alpha \, u^{(0,0,0)} u^{(0,0,1)} 
- 2\, u^{(0,1,0)} u^{(1,0,1)} 
\end{eqnarray}
```
using the heuristic H1 used in beam search for pure quadratization. We will take $\alpha$ and $\gamma$ as 1. After pure quadratization, the ODE set is 
```math
\begin{align}
  \dot{u}^{(1,0,0)} &=  u^{(0,0,1)} u^{(0,1,0)} -  u^{(0,0,0)} u^{(0,1,0)} \\
  &+  u^{(1,0,0)} u^{(1,1,0)} +  u^{(0,0,0)} u^{(1,0,0)} \\
  \dot{u}^{(0,1,0)} &= 3 u^{(0,0,1)} u^{(1,0,0)} +  u^{(0,0,0)} u^{(1,0,0)} \\
  &-  u^{(1,0,0)} u^{(2,0,0)} +  u^{(0,0,0)} u^{(0,1,0)} \\
  \dot{u}^{(0,0,1)} &= -2 u^{(0,0,0)} u^{(0,0,1)} - 2 u^{(0,1,0)} u^{(1,0,1)} \\
  \dot{u}^{(0,0,0)} &= 0 \\
  \dot{u}^{(1,1,0)} &=  u^{(0,1,0)} u^{(0,1,1)} -  u^{(0,1,0)} u^{(0,1,0)}  \\
  &+  u^{(1,1,0)} u^{(1,1,0)} +  u^{(0,1,0)} u^{(1,0,0)} \\ 
  &+ 3 u^{(1,0,0)} u^{(1,0,1)} +  u^{(1,0,0)} u^{(1,0,0)}  \\
  &-  u^{(2,0,0)} u^{(2,0,0)} +  u^{(0,1,0)} u^{(1,0,0)} \\
  \dot{u}^{(2,0,0)} &= 2 u^{(0,1,0)} u^{(1,0,1)} - 2 u^{(0,1,0)} u^{(1,0,0)} \\
  &+ 2 u^{(1,0,0)} u^{(2,1,0)} + 2 u^{(1,0,0)} u^{(1,0,0)} \\
  \dot{u}^{(1,0,1)} &=  u^{(0,0,1)} u^{(0,1,1)} -  u^{(0,0,1)} u^{(0,1,0)} \\
  &+  u^{(1,0,1)} u^{(1,1,0)} +  u^{(0,0,1)} u^{(1,0,0)} \\
  &- 2 u^{(0,0,1)} u^{(1,0,0)} - 2 u^{(1,0,1)} u^{(1,1,0)} \\
  \dot{u}^{(0,1,1)} &= 3 u^{(0,0,1)} u^{(1,0,1)} +  u^{(0,0,1)} u^{(1,0,0)} \\
  &-  u^{(1,0,0)} u^{(2,0,1)} +  u^{(0,0,1)} u^{(0,1,0)} \\ 
  &- 2 u^{(0,0,1)} u^{(0,1,0)} - 2 u^{(0,1,1)} u^{(1,1,0)}
\end{align}
```
```math
\begin{align}
  \dot{u}^{(2,1,0)} &= 2 u^{(0,1,1)} u^{(1,1,0)} - 2 u^{(0,1,0)} u^{(1,1,0)} \\
  &+ 2 u^{(1,1,0)} u^{(2,1,0)} + 2 u^{(1,0,0)} u^{(1,1,0)} \\
  &+ 3 u^{(1,0,0)} u^{(2,0,1)} +  u^{(1,0,0)} u^{(2,0,0)} \\
  &-  u^{(2,0,0)} u^{(3,0,0)} +  u^{(1,0,0)} u^{(1,1,0)} \\
  \dot{u}^{(2,0,1)} &= 2 u^{(0,1,1)} u^{(1,0,1)} - 2 u^{(0,1,0)} u^{(1,0,1)} \\
  &+ 2 u^{(1,0,1)} u^{(2,1,0)} + 2 u^{(1,0,0)} u^{(1,0,1)} \\
  &- 2 u^{(1,0,0)} u^{(1,0,1)} - 2 u^{(1,0,1)} u^{(2,1,0)} \\
  \dot{u}^{(3,0,0)} &= 3 u^{(1,0,1)} u^{(1,1,0)} - 3 u^{(1,0,0)} u^{(1,1,0)} \\
  &+ 3 u^{(2,0,0)} u^{(2,1,0)} + 3 u^{(1,0,0)} u^{(2,0,0)}
\end{align}
```

## Further details

t_max is for choosing the maximum time. The system gives the solution between 0 and t_max in 10 equal intervals.

num_iter is the number of iterations. It determines the number of terms to take from the Taylor expansion of the solution. For convergent series, more terms mean better approximation.

f_num determines the function for which the solution will be given. There are only two choices here. Even if the ODE has more than two unknowns, it is possible to print the solution for only two of them for now.

p_prec is the printing precision. All calculations are performed in exact arithmetic. The rational number is converted to floating point number at the last step using very high number of bits for representation (using facilities of multiprecision library). The printing precision determines how many digits to print after the dot.

After each change of the radioboxes, the output is updated. The calculation is performed from scratch for each change to give a realistic observation of the computation time of prevth. Therefore, even the pure quadratization is performed again for each change of the radioboxes. 

# Install instructions specific to Ubuntu and Ubuntu-like distributions

Get the source of tui-prevth by
> wget https<span/>://github.com/cosargozukirmizi/tui-prevth/archive/refs/heads/main.zip

Unzip it by
> unzip main.zip

Change into directory by
> cd tui-prevth

Install GNU Multiprecision by
> sudo apt install libgmp-dev

Build by 
> cmake -B build .

Compile by
> cd build && make

Run by
> ./a.out
