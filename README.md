# RANS
RANS Modeling of Airfoil

## How to Run
The main file is the flo103.py. This can be run directly without any additional inputs. 

## How to Make Changes to Model
### Input
The current format of the input file is not reproducible. The easiest way will be to change the constructor from the input class to read a new geometry and new parameters. Input.py contains the attributes needed for an input object that a new implementation would also need.

### RANS Model
We are using the Multigrid method to solve a Navier-Stokes problem. To swap out for a different problem, create a new Model subclass with the necessary methods. 

### Integrator
The current integrator used is ImplicitEuler.py. If you implement the constructor and the step() function as a child of the Integrator class, any integration scheme should work. 


## Class Diagram
Link to online source [here](https://lucid.app/lucidchart/3740e36c-b01f-494a-b61c-b08bc9aa8092/edit?invitationId=inv_f79f431d-b5df-4b4d-bf69-b01c7a08117e)

