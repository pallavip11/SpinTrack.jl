# Installation

Please install the official Julia Binary from [here](https://julialang.org/downloads/). PLease download the Long-term support (LTS) release: v1.6.7 (July 19, 2022).

Now, open Julia REPL and input th following commands

``` 
import Pkg;
Pkg.add(url="https://github.com/pallavip11/SpinTrack.jl")
Pkg.add("Plots")
```
To run Julia in tandem with Jupyter notebook, install `IJulia`:

```
using Pkg;
Pkg.add("IJulia")
```
For your own Python/Jupyter installation, just run jupyter notebook yourself in a terminal. Another simple method to run the notebok is:

```
using IJulia
notebook()
```
to launch the IJulia notebook in browser.

On the first time running the notebook, there will be a prompt for whether to install Jupyter. Hit enter to have it use the Conda.jl package to install a minimal Python+Jupyter distribution (via Miniconda) that is private to Julia (not in your PATH).
For more details refer to this [link](https://julialang.github.io/IJulia.jl/stable/manual/running/).

# Running Cases

## Basic Case

To run the most basic example:

``` 
using SpinTrack;
using Plots;

m = FNAL_muon_ring()
sol = get_solution(u1_long(1e-6), m);
plot(sol)
```
NOTE: The semi-colon is used to supress output of the script. Use it according to your convenience.

If the plotting backend `Plots` does not work, try using a different plotting [backend](https://docs.juliaplots.org/latest/backends/). 
## Variations

- To run the simulation for a different time interval, change the value of 'turns' argument in ring parameters as such (default=500):
  ``` 
  m = FNAL_muon_ring()
  m.turns = 200;
  sol = get_solution(u1_long(1e-6), m);
  plot(sol)
  ```  
  This is because ending time is calculated as: no. of turns * length of the ring.
  
- To change other parameters in the code, clone the repository into your own github or make a copy of the code locally. To compile the local copy, use the   following command:
  ``` 
  Pkg.add(path="User/local_path/SpinTrack.jl")
  ```  
- To experiment with using different particles, explore this [file](src/particles.jl).
- To experiment with different ring designs, explore the different designs in this [folder](src/ring_designs).

