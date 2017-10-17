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

- Parsing XML file with chemical reaction data to get relevant paramters

- Handling the calculation of 3+ classes of reaction rate coefficients (e.g. constant, Arrhenius and modified Arrhenius) given the appropriate parameters

- Handling the calculation of progress rates ($\omega_{j}$) and reaction rates ($f_{i}$)  for a system  consisting of $N$ species undergoing $M$ **irreversible**, **elementary** reactions of the above mentioned form

- (Future feature) Handling more classes of reaction rate coefficients and more types of reactions (e.g. reversible and/or non-elementary reactions)

## 2. Installation

TBD

## 3. Basic Usage

### 3.1 Structure of the input file

_**DESCRBIBE THE INPUT FILE**_

Chemical reaction data should be stored in XML format with the following specifications:

1. a \<phase> element with a \<speciesArray> child element which lists the molecular species invovled in the reaction

2. a \<reactionData> element that stores revelant parameters of the reaction:
    - \<reaction reversible> tag indicates whether a reaction is reversible or irreversible, possible values = ["yes", "no"]
    - \<type> tag indicates whether a reaction is elementary or non-elementary, possible values = ["Elementary", "Non-Elementary]
    - \<equation> tag specifies the chemical reaction
    - \<rateCoeff> tag stores parameters relevant to calculate the rate coefficient. It has a child tag indicating the rate coefficient class, possible values = ["\<Constant>", "\<Arrhenius>", "\<modifiedArrhenius>"], and the associated parameters

**Example 3.1.** XML file for the following chemical reaction:

$$
\begin{aligned}
2H_2 + O_2 \quad \longrightarrow \quad 2OH + H_2
\end{aligned}
$$
~~~~
<?xml version="1.0"?>
<ctml>

    <phase>
        <speciesArray> H2 O2 OH HO2 H2O </speciesArray>
    </phase>

    <reactionData id="test_mechanism">

        <!-- reaction 01  -->
        <reaction reversible="no" type="Elementary" id="reaction01">
            <equation>2H2 + O2 [=] 2OH + H2</equation>
            <rateCoeff>
                <modifiedArrhenius>
                    <A units="m3/mol/s">1e+8</A>
                    <b>0.5</b>
                    <E units="J/mol">5e+04</E>
                </modifiedArrhenius>
            </rateCoeff>
            <reactants>H2:2 O2:1</reactants>
            <products>OH:2 H2:1</products>
        </reaction>
    </reactionData>

</ctml>
~~~~

### 3.2 Reading the input file

- `XmlParser(xml_file)` produces a list of reaction data from XML file "xml_file" contents.

- `

### 3.3. Calculating the reaction rate coefficient



### 3.4 Calculating the progress rate

### 3.5 Calculating the reaction rate

## 4. Examples

The neccessary code can be found at and downloaded from [here](https://github.com/G12-cs207-FinalProject/cs207-FinalProject).

Once you open the directory with the downloaded files, please note the following items:

- The main module file name is called **chemkin.&#8203;py**
- The test suite is called **test_chemkin.py**.
- Additional files include execution code stored in **run_chemkin.py** and a set of demo *.xml* reaction files.

You can run the test suite on your local machine by typing the following command in your command line when in the directory of the downloaded files.
```sh
$ pytest
```

### 4.1. A system of irreversible, elementary chemical reactions

- The system of chemical reactions of interest:

$$
\begin{aligned}
2H_2 + O_2 &\longrightarrow 2OH + H_2 \\
OH + HO_2 &\longrightarrow H_2O + O_2 \\
H_2O + O_2 &\longrightarrow HO_2 + OH
\end{aligned}
$$

- Input XML file:
~~~
<?xml version="1.0"?>

<ctml>
    <phase>
        <speciesArray> H2 O2 OH HO2 H2O </speciesArray>
    </phase>

    <reactionData id="test_mechanism">
        <!-- reaction 01  -->
        <reaction reversible="no" type="Elementary" id="reaction01">
            <equation>2H2 + O2 [=] 2OH + H2</equation>
            <rateCoeff>
                <modifiedArrhenius>
                    <A units="m3/mol/s">1e+8</A>
                    <b>0.5</b>
                    <E units="J/mol">5e+04</E>
                </modifiedArrhenius>
            </rateCoeff>
            <reactants>H2:2 O2:1</reactants>
            <products>OH:2 H2:1</products>
        </reaction>

        <!-- reaction 02 -->
        <reaction reversible="no" type="Elementary" id="reaction02">
            <equation>OH + HO2 [=] H2O + O2</equation>
            <rateCoeff>
                <Constant>
                    <k>1e+4</k>
                </Constant>
            </rateCoeff>
            <reactants>OH:1 HO2:1</reactants>
            <products>H2O:1 O2:1</products>
        </reaction>

        <!-- reaction 03 -->
        <reaction reversible="no" type="Elementary" id="reaction03">
            <equation>H2O + O2 [=] HO2 + OH</equation>
            <rateCoeff>
                <Arrhenius>
                    <A units="m3/mol/s">1e+7</A>
                    <E units="J/mol">1e+04</E>
                </Arrhenius>
            </rateCoeff>
            <reactants>H2O:1 O2:1</reactants>
            <products>HO2:1 OH:1</products>
        </reaction>
    </reactionData>
</ctml>
~~~

- Given the the species concentration $xi = [2.0, 1.0, 0.5, 1.0, 1.0]$ in units of $m^{3}/(mol\cdot s)$, calculate the reaction rate of each species at $Ti = [750, 1500, 2500]$ in units of K.

~~~
from chemkin import *

xml_file = './xml-files/rxns_hw5.xml'
xml_parser = XmlParser(xml_file)

species, rxn_data_list = xml_parser.load() # load data
n_species = len(species) # number of species

xi = [2.0, 1.0, 0.5, 1.0, 1.0] # specie concentration

species_idx_dict = {} # build the dictionary of key = species_name, value = species_index
for i, s in enumerate(species):
	species_idx_dict[s] = i

Ti = [750, 1500, 2500]

for T in Ti:
	sys_vi_p = [] # list of reactant Stoichiometric coefficients in each rxn
	sys_vi_dp = [] # list of product Stoichiometric coefficients in each rxn
	ki = [] # list of reation rate coefficients in each rxn

	for rxn_data in rxn_data_list: # 1 rxn per rxn_data
		if rxn_data.type != RxnType.Elementary: # skip non-elementary reactions for now
			continue
		if rxn_data.reversible != False: # skip reversible reactions for now
			continue
		
		rxn_id = rxn_data.rxn_id # save id

		rxn_vi_p = np.zeros((n_species,)) # save the Stoichiometric coefficients of the reactants in this rxn
		for s, vi in rxn_data.reactants.items():
			idx = species_idx_dict[s] # get index of the specii
			rxn_vi_p[idx] = vi
		sys_vi_p.append(list(rxn_vi_p))

		rxn_vi_dp = np.zeros((n_species,)) # save the Stoichiometric coefficients of the products in this rxn
		for s, vi in rxn_data.products.items():
			idx = species_idx_dict[s] # get index of the specii
			rxn_vi_dp[idx] = vi
		sys_vi_dp.append(list(rxn_vi_dp))
		
		coef_params = rxn_data.rate_coeff
		if isinstance(coef_params, list):
			if len(coef_params) == 3: # modified arrhenius coef
				A = coef_params[0]
				b = coef_params[1]
				E = coef_params[2]
				ki.append(ModArrheniusCoef(A, b, E, T).get_coef())
			else: # arrhenius coef
				A = coef_params[0]
				E = coef_params[1]
				ki.append(ArrheniusCoef(A, E, T).get_coef())
		else: # const coef
			ki.append(ConstCoef(coef_params).get_coef())

	rxn_rates = IrrevElemRxn(ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()
	
	print('------At Temperature', T, 'K------')
	for s, rate in zip(species, rxn_rates):
		print('    ', s, ':', rate)
	print('--------------------------------')
~~~

- The expected result:
~~~
-------At Temperature 750 K------
     H2 : -3607077.8728
     O2 : -5613545.18362
     OH : 9220623.05642
     HO2 : 2006467.31082
     H2O : -2006467.31082
--------------------------------
------At Temperature 1500 K------
     H2 : -281117620.765
     O2 : -285597559.238
     OH : 566715180.003
     HO2 : 4479938.47318
     H2O : -4479938.47318
--------------------------------
------At Temperature 2500 K------
     H2 : -1804261425.96
     O2 : -1810437356.94
     OH : 3614698782.9
     HO2 : 6175930.97566
     H2O : -6175930.97566
--------------------------------
~~~
