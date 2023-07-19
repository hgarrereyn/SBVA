# SBVA (Structured Bounded Variable Addition)

SBVA is a tool for reducing SAT formulas using _structured bounded variable addition_.

Read our SAT'23 paper: [Effective Auxiliary Variables via Structured Reencoding](https://arxiv.org/pdf/2307.01904.pdf) (preprint)

🏆 SBVA-CaDiCaL (CaDiCaL w/ SBVA as preprocessor) was the overall main-track winner of the [2023 SAT competition](https://satcompetition.github.io/2023/)! (1st overall, 1st satisfiable, 2nd unsatisfiable)

## Installation

SBVA requires Eigen as a dependency. To compile, run the following:

```
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar xf eigen-3.4.0.tar.gz
make
```

## Usage

### As a preprocessor

```
./sbva [-i input] [-o output] [-p proof] [-t timeout] [-s max_replacements] [-v]
```

Options:
* `-i`: specify input file (default is stdin)
* `-o`: specify output file (default is stdout)
* `-p`: if specified, save a DRAT proof of the transformation to this file
* `-t`: if specified, the BVA algorithm will exit after this many seconds and the semi-reduced formula will be returned
* `-s`: if specified, limits the number of replacements / new auxiliary variables
* `-v`: enable verbose logging

Examples:

```sh
# Reduce a formula
./sbva -i problem.cnf -o out.cnf

# Reduce a formula and generate a proof
./sbva -i problem.cnf -o out.cnf -p proof.drat
```

### With a solver

A wrapper script is provided to run SBVA along with a SAT solver and automatically fixup the resulting model (if SAT) or DRAT proof (if UNSAT).

On certain types of problems, SBVA itself can take a long time to run even if the original formula would solve quickly in a SAT solver, so the wrapper script also suports running SBVA with a timeout and falling back to running the original formula.

For sane defaults, use the `sbva_wrapped` script:

```
./sbva_wrapped [solver] [input.cnf] [output.proof]
```

The solver will be invoked like:
```
<solver> <input.cnf> <output.drat> --no-binary
```

(CaDiCaL and Kissat use this format)

The result of the solver will be printed to `stdout`, e.g.

```
s SATISFIABLE
<model>
```

or

```
s UNSATISFIABLE
```

If the problem is unsatisfiable, a DRAT proof will be saved to `[output.proof]`.

Examples:
```sh
# Run SBVA with CaDiCaL
./sbva_wrapped cadical problem.cnf output.proof
```

You can also invoke the `wrapper.py` script directly, e.g.:

```
python3 wrapper.py \
    --input problem.bva \
    --output proof.out \
    --bva ./sbva \
    --t1 200 \
    --t2 400 \
    --solver cadical
```

The options `t1` and `t2` control the inner and outer timeouts (in seconds) respectively. SBVA will attempt to finish gracefully after `t1` seconds, but if it gets stuck in a busy section of the algorithm, an outer timeout will kill it after `t2` seconds.

### Example Formulas

You can find example CNF formulas (on which SBVA is effective) in the `examples/` directory. These are from the [Global Benchmark Database](https://benchmark-database.de/) and the [packing k-coloring problem](https://arxiv.org/abs/2301.09757).


## Authors

SBVA was developed by Andrew Haberlandt and Harrison Green with advice from Marijn Heule.
