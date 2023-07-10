
import argparse
import subprocess
from tempfile import NamedTemporaryFile
import sys

def find_stats(s):
    for line in s.split(b"\n"):
        line = line.decode()
        if line[0] == 'p':
            nv, nc = map(int, line.split()[2:])
            return nv, nc

    return (-1,-1)


def run_solver_bare(args, cnf_str):
    argv = [
        args.solver,
        "-",
    ]
    if args.output:
        argv += [args.output]
    argv += ["--no-binary"]

    subprocess.run(argv, input=cnf_str)

def run_solver_reduced(args, reduced_cnf, reduced_drat, orig_cnf_str):
    with NamedTemporaryFile() as solver_drat:
        p = subprocess.run([
            args.solver,
            reduced_cnf,
            solver_drat.name,
            '--no-binary'
        ], stdout=subprocess.PIPE)

        rcode = p.returncode

        if rcode not in [0,10,20]:
            # If solver failed, run on original input
            run_solver_bare(args)
        else:
            res = p.stdout.decode('latin-1')

            if rcode == 20:
                # UNSAT: Fixup the DRAT proof
                proof = ''
                with open(reduced_drat, 'r') as f:
                    proof += f.read()
                with open(solver_drat.name, 'r') as f:
                    proof += f.read()

                if args.output:
                    with open(args.output, 'w') as f:
                        f.write(proof)

                print('s UNSATISFIABLE')
            else:
                # SAT: Remove auxiliary variables from solution
                nv, _ = find_stats(orig_cnf_str)

                lits = []
                for line in res.split('\n'):
                    if line.startswith('v'):
                        line_lits = line.split(' ')[1:]
                        lits += [int(l) for l in line_lits if abs(int(l)) <= nv]
                out = "s SATISFIABLE\nv " + " ".join(map(str, lits)) + "\n"

                print(out)


def run(args):
    if args.input:
        with open(args.input, "rb") as f:
            input_cnf_str = f.read()
    else:
        input_cnf_str = sys.stdin.buffer.read()

    stats = find_stats(input_cnf_str)
    print(f"c before: {stats}")

    extra_opts = []
    if args.no_heuristic:
        extra_opts.append('-n')

    if args.limit > 0:
        extra_opts += ["-s", str(args.limit)]

    with NamedTemporaryFile() as bva_cnf, NamedTemporaryFile() as bva_drat:
        # Run BVA
        p = subprocess.run([
            'timeout', str(args.t2), args.bva,
            '-o', bva_cnf.name,
            '-p', bva_drat.name,
            '-t', str(args.t1)
        ] + extra_opts, input=input_cnf_str)

        with open(bva_cnf.name, "rb") as f:
            stats = find_stats(f.read())
            print(f"c after: {stats}")

        if p.returncode == 0:
            # BVA ran successfully
            run_solver_reduced(args, bva_cnf.name, bva_drat.name, input_cnf_str)
        else:
            # Run original solver on input
            run_solver_bare(args)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-o', '--output', type=str)
    parser.add_argument('--bva', type=str, required=True)
    parser.add_argument('--t1', type=int, help='Inner timeout', required=True)
    parser.add_argument('--t2', type=int, help='Outer timeout', required=True)
    parser.add_argument('--solver', type=str, required=True)
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--no-heuristic', action="store_true")
    args = parser.parse_args()
    run(args)
