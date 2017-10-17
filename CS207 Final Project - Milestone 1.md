# Chemical Kinetics Module

## 1. Introduction
`chemkin` stands for chemical kinetics and is an objected-oriented library for modeling kinetics of chemical reactions.

### 1.1 Key Chemical Kinetics Concepts

Chemical kinetics is the study of rates of chemical processes such as reaction rates, reaction mechanisms, etc. as well as the construction of mathematical models that can describe the characteristics of a chemical reaction ([wikipedia](https://en.wikipedia.org/wiki/Chemical_kinetics)). Typically, amounts of molecular species reacted (consumed)/formed and the rates of their consumption/formation are of interest.

Chemical reactions can be categorized as elementary or non-elementary and reversible or irreversible. Each class of chemical reactions can be modeled by ODEs (ordinary differential equations) or PDEs (partial differntial equations) using different strategies.

For instance, for a system  consisting of $N$ species undergoing $M$ **irreversible**, **elementary** reactions of the form:

$$\sum_{i=1}^{N}{\nu_{ij}^{\prime}{S}_{i}} \longrightarrow 
  \sum_{i=1}^{N}{\nu_{ij}^{\prime\prime}{S}_{i}}, \qquad \text{for } j = 1, \ldots, M$$

where\
$S_{i}$ = Chemical symbol of specie $i$ \
$\nu_{ij}^{\prime}$ = Stoichiometric coefficients of reactants \
$\nu_{ij}^{\prime\prime}$ = Stoichiometric coefficients of products

The rate of change of species $i$ (i.e. the reaction rate of species $i$) can be written as
  
$$\begin{aligned}
  f_{i} = \frac{d[i]}{dt} = \sum_{j=1}^{M}{\nu_{ij}\omega_{j}}, \qquad \text{for } i = 1, \ldots, N
\end{aligned}$$

where\
$\nu_{ij}$ = $\nu_{ij}^{\prime\prime}$ - $\nu_{ij}^{\prime}$ \
$\omega_{j}$ = the progress rate for each reaction


The progress rate $\omega_{j}$ is given by 

$$\begin{aligned}
  \omega_{j} = k_{j}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}}, \qquad j = 1, \ldots, M
\end{aligned}
$$

where \
$k_{j}$ = the forward reaction rate coefficient \
$x_{i}$ = the Concentration of specie $i$

Note: A complete table of notation
| Symbol | Meaning |
|:--------:|:-------:|
| $\mathcal{S}_{i}$ | Chemical symbol of specie $i$ |
| $\nu_{ij}^{\prime}$ | Stoichiometric coefficients of reactants |
| $\nu_{ij}^{\prime\prime}$ | Stoichiometric coefficients of products |
| $N$                       | Number of species in system |
| $M$                       | Number of elementary reactions |
| $f_{i}$                   | Rate of consumption or formation of specie $i$ (reaction rate) |
| $\omega_{j}$              | Progress rate of reaction $j$ |
| $x_{i}$                   | Concentration of specie $i$ |
| $k_{j}$                   | Reaction rate coefficient for reaction $j$ |

### 1.2 The `chemkin` Module

The high level functionality of the **chemkin** module is to take an XML file with reaction data as input and outputs the RHS (right-hand-side) of the ODE describing the rate of change of all molecular species involved in the chemical reaction(s) of interest (i.e. $\sum_{j=1}^{M}{\nu_{ij}\omega_{j}}, \text{ for } i = 1, \ldots, N$).

Features of the **chemkin** module include:

- Parsing XML file with reaction data

- Handling the calculation of 3+ classes of reaction rate coefficients (e.g. constant, Arrhenius and modified Arrhenius) given the appropriate parameters

- Handling the calculation of progress rates ($\omega_{j}$) and reaction rates ($f_{i}$)  for a system  consisting of $N$ species undergoing $M$ **irreversible**, **elementary** reactions of the above mentioned form

- (Future feature) Handling more classes of reaction rate coefficients and more types of reactions (e.g. reversible and/or non-elementary reactions)

## 2. Installation
TBD

## 3. Basic Usage and Examples

### 3.1 Basic Usage

### 3.2 Examples
