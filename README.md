# RANS
RANS Modeling of Airfoil. An overview of our project is given in the research project pdf and the code is based off of Gigi's (Luigi Martinelli) thesis. 

## How to Run
The main file is the flo103.py. This can be run directly without any additional inputs. There are 2 sample airfoils in the main repository and they can be changed directly by changing the name in flo103.py. 

## How to Make Changes to Model
### Input
The current format of the input file is not reproducible. The easiest way will be to change the constructor from the input class to read a new geometry and new parameters. Input.py contains the attributes needed for an input object that a new implementation would also need.

### RANS Model
We are using the Multigrid method to solve a Navier-Stokes problem. To swap out for a different problem, create a new Model subclass with the necessary methods. 

### Integrator
The current integrator used is ImplicitEuler.py. If you implement the constructor and the step() function as a child of the Integrator class, any integration scheme should work. 

## Final Paper
The final paper (latex document) can be found [here](https://www.overleaf.com/project/61bf84dbe37b215a352b7f58)


## UML Diagram
The UML Diagram can be found [here](https://github.com/andybroth/RANS/blob/e60dc63318b2e0e1277bf4da68918732ce84af3d/Git%20Stokesed%20UML%20Diagram.PNG)


## Documentation
To generate documentation, `cd` into the `docs/` directory where there should be a `conf.py` file. (Install Sphinx if it is not already installed.) From there run the command `make html` to build html documentation or `make latexpdf` for latex files (you need Latex installed). 

If you need to build the documentation from absolute scratch (i.e. there is no `docs/` directory), run the command `sphinx-quickstart docs` in the `RANS/` directory. Then `cd docs/` and open `conf.py`. About line 33, change `extensions = []` to `extensions = ['autoapi.extension', 'sphinx.ext.napoleon']` and add a line below that with `autoapi_dirs = ['../bin', '../tests']`. Save your changes and then you can `make html` and `make latexpdf`.

If necessary, documentation can be found in the `ourdocs/` directory. The `html` file is in `ourdocs/_build/html/index.html` and the `latex` file is in `ourdocs/_build/latex/ransmodelingofairfoils.tex`.
