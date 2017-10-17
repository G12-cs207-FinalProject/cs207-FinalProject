# Chemical Kinetics Module

# Introduction
The following module allows users to calculate reaction rates for a set of species in a set of reactions provided by the user.

Currently, the module is functional for a system consisting of $$N$$ species undergoing $$M$$ **irreversible**, **elementary** reactions of the form 

$$\sum_{i=1}^{N}{\nu_{ij}^{\prime}{S}_{i}} \longrightarrow 
  \sum_{i=1}^{N}{\nu_{ij}^{\prime\prime}{S}_{i}}, \qquad j = 1, \ldots, M$$
  
 The rate of change of specie $$i$$ (the reaction rate) can be written as:
  
$$f_{i} = \sum_{j=1}^{M}{\nu_{ij}\omega_{j}}, \qquad i = 1, \ldots, N$$

Where the progress rate for each reaction is given by :

$$\omega_{j} = k_{j}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}}, \qquad j = 1, \ldots, M$$

and $$k_{j}$$ is the forward reaction rate coefficient.

# Installation

The neccessary code can be found at and downloaded from [here](https://github.com/G12-cs207-FinalProject/cs207-FinalProject).

Once you open the directory with the downloaded files, please note the following items:
- The main module file name is called **chemkin.&#8203;py** 
- The test suite is called **test_chemkin.py**. 
- Additional files include execution code stored in **run_chemkin.py** and a set of demo *.xml* reaction files.

You can run the test suite on your local machine by typing the following command in your command line when in the directory of the downloaded files.
```sh
$ pytest
```

# Basic Examples and Usage

When given a *.xml* file containing reaction data, the execution module will parse and extract reaction data from the file and calculate the reaction rates. Current functionality covers elementary, irreversible reactions. The user is expected to modify **run_chemkin.py** to include his own defined temperatures and species concentrations as well as path to the *.xml* containing the reaction data.

- Example $$1$$: Given a set of $$3$$ elementary, irreversible reactions, calculate the reaction rate for each species given the species concentration $$xi = [2.0, 1.0, 0.5, 1.0, 1.0]$$ in $$m^{3}/mol/s$$ at $$3$$ different temperatures $$Ti = [750, 1500, 2500]$$
-- Make sure you have **./xml-files/rxns_hw5.xml** in your current directory, then open and execute **run_chemkin.py** module in any Python  IDE environment supporting .py file type.
- Output are the reaction rates as follows:

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

 - Example $$2$$: Given a set of $$3$$ elementary, irreversible reactions, calculate the reaction rate for each species given the species concentration $$xi = [2.0, 1.0, 0.5, 1.0, 1.0]$$ in $$m^{3}/mol/s$$ at $$3$$ different temperatures $$Ti = [750, 1500, 2500]$$
-- Make sure you have **./xml-files/rxns_ideal.xml** in your current directory, then open and execute **run_chemkin.py** module in any Python  IDE environment supporting .py file type.
  -- 

------At Temperature 750 K------
     H : -749299.304081
     O : 749299.304081
     OH : 749299.305572
     H2 : -0.000745400452028
     O2 : -749299.304826
--------------------------------
------At Temperature 1500 K------
     H : -229675142.445
     O : 229675142.445
     OH : 229675142.457
     H2 : -0.00614143817624
     O2 : -229675142.451
--------------------------------
------At Temperature 2500 K------
     H : -2268284161.54
     O : 2268284161.54
     OH : 2268284161.57
     H2 : -0.0142765166881
     O2 : -2268284161.55
--------------------------------

```sh
$ npm install --production
$ NODE_ENV=production node app
```

### Plugins

Dillinger is currently extended with the following plugins. Instructions on how to use them in your own application are linked below.

| Plugin | README |
| ------ | ------ |
| Dropbox | [plugins/dropbox/README.md] [PlDb] |
| Github | [plugins/github/README.md] [PlGh] |
| Google Drive | [plugins/googledrive/README.md] [PlGd] |
| OneDrive | [plugins/onedrive/README.md] [PlOd] |
| Medium | [plugins/medium/README.md] [PlMe] |
| Google Analytics | [plugins/googleanalytics/README.md] [PlGa] |





