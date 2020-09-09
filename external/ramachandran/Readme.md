# Ramachandran-Angle-Probability.py
-----------------------------------

A simple python class to generate random angles based on the angle density extracted experimental data.
Two classes are provided. One class reads from a data file in (Py)Rosetta and the other one reads from
Dorn's Angle probability list.

## Usage
--------

A simple example for (Py)Rosetta that uses matplotlib to generate plots for each aminoacid and a Ramachandran plot:

``` Python
    import matplotlib.pyplot as plt

    rama = Ramachandran()
    rama.load(path="/home/h3nnn4n/PyRosetta4.Release.python35.linux.release-147/"
                   "setup/pyrosetta/database/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb")
    rama.process()

    for k, v in rama.parsed_data.items():
        plt.imshow(v.transpose(), cmap='viridis', interpolation='nearest', origin="lower")
        plt.savefig(k + ".png")
```

To get a random angle for an amino acid one can use `phi, psi = rama.get_random_angle('T')`. Note that either the 1 or 3 letter name works.

This program relies on data released with (Py)Rosetta. Rosetta can be found here and it is free for academic use: https://www.rosettacommons.org/

For dorn's a simple example is show bellow:

``` Python
    import APL

    APL = apl()
    apl.load()

    import matplotlib.pyplot as plt
    for k, v in apl.data.items():
        print(k, v)
        plt.imshow(v.transpose(), cmap='viridis', interpolation='nearest', origin="lower")
        plt.savefig("%s_%s.png" % (k))
```

Random angles can be generated using `phi, psi = apl.get_random_angle('GLU', 'H')`. Note that wither the 1 or 3 letter name works.
The second argument is the secondary structure.

More info and APL can be found on the links below:
https://www.ncbi.nlm.nih.gov/pubmed/26495908
http://sbcb.inf.ufrgs.br/apl


LICENSE
-------

Only the files in this repository are released under the MIT liecense. Check the LICENSE file for more details.
The files belonging to (Py)Rosetta have their own licenses.
