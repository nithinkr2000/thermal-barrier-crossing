# thermal-barrier-crossing
A package to simulate the crossing of a potential barrier via random thermal fluctuations. Intended to prove that the initial conditions can be optimized to achieve faster results. 

## Installation / Compilation

To compile the code, go to folder `scripts`, and execute :

```bash
g++ -fopenmp main_func.cpp -o 1dmcmc
```
This will create an executable called `1dmcmc`.

## Example execution

After successful compilation, you can run the executable that you created during compilation. (If you followed the instructions above, it will be called `1dmcmc`.) To run the example files, go to subfolder `scripts`.

Using the example files in this repo, the output will be written to a subfolder, called `example_run`, as specified in the parameter file `pg0`. Therefore, first create this directory:

```bash
mkdir example_run
```

Afterwards, run the program:
```bash
./1dmcmc
```

This program needs three input arguments. You will be prompted to provide those at execution.


1. First, you will be prompted for indicating the filename. To use the example files, type `pg` and _Enter_. 

2. Next choose the exchange function to be used: Either TREX `t` or HREX `h`, confirm with _Enter_. Afterwards, the program will be running. 

3. Lastly, indicate the number of replicas. In the provided example, you should prompt `6`. By default, you will be using as many OpenMP threads, as there are replicas in your setup.

The default values in the example will lead to a very short simulation, which finishes within 1 s. An appropriate choice for the number of steps and the number of executions in the input file would be:
```bash
n_steps=5000
n_ex=1000
```

