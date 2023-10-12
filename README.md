# tui-prevth
tui for probabilistic evolution theory solution of ordinary differential equations with multinomial right-hand-sides


This program solves the initial value problem of ordinary differential equations with multinomial right-hand-side functions. Four different initial value problems are put in the source code. The code can be manipulated to accomodate other initial value problems with multinomial right-hand-sides. 

![tui-prevth](Screenshot_2023-08-22_12-21-22.png?raw=true "tui-prevth")

The calculations are performed using rational arithmetic. In this way, error accumulation because of finite nature of double arithmetic is avoided. 

The command line interface has five sets of input. 

ODE preset is for choosing the initial value problem. The ODEs are van der Pol, quartic anharmonic oscillator, Henon-Heiles and Rabinovich-Fabrikant. All initial values for all functions are taken as 1/2 (This can be changed easily within the source code, but currently not from the tui).

The ODEs are given below in LaTeX.

vdPol
-----
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& 1/3\, u^{(1,0)} u^{(2,0)} \nonumber \\
 &-& 1\, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& 1\, u^{(0,0)} u^{(1,0)}
\end{eqnarray}
```
QuartAnhOsc
-----------
```math
\begin{eqnarray}
  \dot{u}^{(1,0)} &=& 1\, u^{(0,0)} u^{(0,1)}\\
  \dot{u}^{(0,1)} &=& -1\, u^{(0,0)} u^{(1,0)} \nonumber \\
 &-& 1\, u^{(1,0)} u^{(2,0)}
\end{eqnarray}
```
HenHeil
-------
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,1,0,0)}\\
  \dot{u}^{(0,1,0,0)} &=& -1\, u^{(0,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &-& 2\, u^{(0,0,1,0)} u^{(1,0,0,0)}\\
  \dot{u}^{(0,0,1,0)} &=& 1\, u^{(0,0,0,0)} u^{(0,0,0,1)}\\
  \dot{u}^{(0,0,0,1)} &=& -1\, u^{(0,0,0,0)} u^{(0,0,1,0)} \nonumber \\
 &-& 1\, u^{(1,0,0,0)} u^{(1,0,0,0)} \nonumber \\
 &+& 1\, u^{(0,0,1,0)} u^{(0,0,1,0)}
\end{eqnarray}
```
Rabinovich-Fabrikant
--------------------
```math
\begin{eqnarray}
  \dot{u}^{(1,0,0)} &=& 1\, u^{(0,0,1)} u^{(0,1,0)} \nonumber \\
 &-& 1\, u^{(0,0,0)} u^{(0,1,0)} \nonumber \\
 &+& 1\, u^{(1,0,0)} u^{(1,1,0)} \nonumber \\
 &+& 1\, u^{(0,0,0)} u^{(1,0,0)}\\
  \dot{u}^{(0,1,0)} &=& 3\, u^{(0,0,1)} u^{(1,0,0)} \nonumber \\
 &+& 1\, u^{(0,0,0)} u^{(1,0,0)} \nonumber \\
 &-& 1\, u^{(1,0,0)} u^{(2,0,0)} \nonumber \\
 &+& 1\, u^{(0,0,0)} u^{(0,1,0)}\\
  \dot{u}^{(0,0,1)} &=& -2\, u^{(0,0,0)} u^{(0,0,1)} \nonumber \\
 &-& 2\, u^{(0,1,0)} u^{(1,0,1)}
\end{eqnarray}
```

Pure quadratization is performed by beam search with $\omega$ as 1. Then solution is obtained by probabilistic evolution theory.

t_max is for choosing the maximum time. The system gives the solution between 0 and t_max in 10 equal intervals.

num_iter is the number of iterations. It determines the number of terms to take from the Taylor expansion of the solution. For convergent series, more terms mean better approximation.

f_num determines the function for which the solution will be given. There are only two choices here. Even if the ODE has more than two unknowns, it is possible to print the solution for only two of them for now.

p_prec is the printing precision. All calculations are performed in exact arithmetic. The rational number is converted to floating point number at the last step using very high number of bits for representation (using facilities of multiprecision library). The printing precision determines how many digits to print after the dot.

After each change of the radioboxes, the output is updated. The calculation is performed from scratch for each change to give a realistic observation of the computation time of prevth. Therefore, even the pure quadratization is performed again for each change of the radioboxes. 

# Install instructions specific to Ubuntu and Ubuntu-like distributions

Get the source of ftxui-starter by
> wget https<span/>://github.com/ArthurSonzogni/ftxui-starter/archive/refs/heads/master.zip

Unzip it by
> unzip master.zip

Go to the directory to put the source of tui-prevth by 
> cd ftxui-starter-master/src/

Remove the current main.cpp by
> rm main.cpp

Get the main.cpp for tui-prevth by
> wget https<span/>://raw.githubusercontent.com/cosargozukirmizi/tui-prevth/main/main.cpp

Install GNU Multiprecision by
> sudo apt install libgmp-dev 

Go to the directory above by
> cd ..

Add a line to CMakeLists.txt by
> sed -i '2i set(CMAKE_CXX_STANDARD 20)' CMakeLists.txt

Add a target link library to CMakeLists.txt by
> sed -i '33i gmpxx' CMakeLists.txt

Add another target link library to CMakeLists.txt by
> sed -i '33i gmp' CMakeLists.txt

Make a build directory by
> mkdir build

Change into the build directory by
> cd build

Then build by
> cmake ..

Then compile by
> make -j

Then run by
> ./ftxui-starter
