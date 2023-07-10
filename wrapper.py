
import argparse
import subprocess
from tempfile import NamedTemporaryFile

def find_stats(fname):
    with open(fname, 'r') as f:
        for line in f:
            if line[0] == 'p':
                nv, nc = map(int, line.split()[2:])
                return nv, nc

    return (-1,-1)


def run_solver_bare(args):
    subprocess.run([
        args.solver,
        args.input,
        args.output,
        '--no-binary'
    ])


def run_solver_reduced(args, reduced_cnf, reduced_drat):
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

                with open(args.output, 'w') as f:
                    f.write(proof)

                print('s UNSATISFIABLE')
            else:
                # SAT: Remove auxiliary variables from solution
                nv, _ = find_stats(args.input)

                lits = []
                for line in res.split('\n'):
                    if line.startswith('v'):
                        line_lits = line.split(' ')[1:]
                        lits += [int(l) for l in line_lits if abs(int(l)) <= nv]
                out = "s SATISFIABLE\nv " + " ".join(map(str, lits)) + "\n"

                print(out)


def run(args):
    with NamedTemporaryFile() as bva_cnf, NamedTemporaryFile() as bva_drat:
        # Run BVA
        p = subprocess.run([
            'timeout', str(args.t2), args.bva,
            '-i', args.input,
            '-o', bva_cnf.name,
            '-p', bva_drat.name,
            '-t', str(args.t1),
        ])

        if p.returncode == 0:
            # BVA ran successfully
            run_solver_reduced(args, bva_cnf.name, bva_drat.name)
        else:
            # Run original solver on input
            run_solver_bare(args)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('--bva', type=str, required=True)
    parser.add_argument('--t1', type=int, help='Inner timeout', required=True)
    parser.add_argument('--t2', type=int, help='Outer timeout', required=True)
    parser.add_argument('--solver', type=str, required=True)
    args = parser.parse_args()
    run(args)
