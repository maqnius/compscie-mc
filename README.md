[![CircleCI](https://circleci.com/gh/maqnius/compscie-mc/tree/master.svg?style=svg)](https://circleci.com/gh/maqnius/compscie-mc/tree/master)
[![codecov](https://codecov.io/gh/maqnius/compscie-mc/branch/master/graph/badge.svg)](https://codecov.io/gh/maqnius/compscie-mc/)


# Particlesim
Simple tool for particle markov simulations and simulated annealing using Ewald Summation.

## Setup

To install the package, simply go to the package directory and 
run `python setup.py install` with in terminal.

## Initialization
There are basically two ways to initialize a simulation:

### Using a csv input file
+ The basic parameters are set in a config file that looks like this

```
[general]
k_cutoff = 1.0
sigma_ewald = 1.0
r_cutoff = 3.0
box-size = 20.0

[manual]
csv_path = /path/to/csvfile/example_init_positions.
[ewald_summation]
use_ewald = yes
sigma = 1.0

[lennard_jones]
use_lennard_jones = yes

```

+ The initial particle configuration is defined by a csv file at the location given in the config and should look similar to this example

```
label,x,y,z,Charge,LJ-Epsilon,LJ-Sigma
Na,1,0,0,1,0.5,1
Cl,2,2,1,-1,0.5,3
Cl,1,0,1,-1,1,1
Na,1,1,1,1,0.5,1

```

+ Then the following commands are necessary to initialize the setup

```python
# Get a ProblemCreator file from the config file
creator = ProblemCreator(test_config_path)

# Initialize a SystemConfiguration Object from the config file
system_config = creator.generate_problem()
```

### Using a config input file
+ The basic parameters are set in a config file similar to the previous case but without a `[manual]` section.
But instead, it contains sections for each particle type beginning with `particle_class_*`. Each class will be collected in the following.

```
[general]
k_cutoff = 1.0
Sigma_ewald = 1.0
r_cutoff = 3.0
box-size = 20.0

[particle_class_1]
label = Natrium
type = NA
charge = 1
distribution = uniform
number = 10

[particle_class_2]
label = Chlor
type = CL
charge = -1
distribution = uniform
number = 20

```
The Lennard-Jones parameters for the atom types are defined in lib/particle_types.py. This file can be edited to fit your
purpose. In the following, a particle setup is created according to your statements for the attributes `distribution` and `number`.
(At this point, only the creation of `uniform` particle distributions is supported.)

+ The initialization process is equal to the previous method

```python
# Get a ProblemCreator file from the config file
creator = ProblemCreator(test_config_path)

# Initialize a SystemConfiguration Object from the config file
system_config = creator.generate_problem()
```

## Run

To run the simulation simply create a sampler with
```
sampler = particlesim.api.Sampler(system_configuration=system_conf)
```

And run one simple metropolis simulation with
```
xyz_trajectory, potential = sampler.metropolis(iteration_number=1000,beta=100)

```

Or metropolis with simulated annealing
```
xyz_trajectory, potential = sampler.metropolis(iteration_number=1000, step=0.1, beta=100)
```

`beta` can be a single value - in this case it gets interpreted as the maximum value -,
a tuple that defines the minimum and maximum value or a list with a length of `iteration_number` to define
a beta value for each simulation step.

## Evaluation
At this point a trajectory export to xyz-Format and csv output is supported and can be accessed with the `particlesim.utils.xyz`
and `particlesim.config_parser` module.
