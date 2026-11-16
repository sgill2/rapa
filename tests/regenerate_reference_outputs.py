"""
Regenerate examples/solutions_to_compare/ with the current code.

Run this only after a deliberate change to RAPA's results, then review the diff before committing:

    python tests/regenerate_reference_outputs.py
    git diff --stat examples/solutions_to_compare

Each example is run in a temporary folder (with -l and --pymol, as the reference runs were) and absolute paths in
the text outputs are replaced by '.', so the references contain no machine-specific paths.
"""

import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

import rapa  # noqa: E402

# reference name -> input PDB in examples/
EXAMPLES = {'1bcd': '1bcd_noH.pdb', '1xl2': '1xl2_noH.pdb'}
REFERENCE_DIR = os.path.join(REPO, 'examples', 'solutions_to_compare')
TEXT_SUFFIXES = ('.info', '.log', '.pml', '.json', '.csv')


def run_example(name, pdb_name, workdir):
    """Run one example in workdir and return its output folder."""
    shutil.copy(os.path.join(REPO, 'examples', pdb_name), workdir)
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        rapa.run(pdb_name, out_prefix=name, output_dir=f'{name}_outputs', log=True, pymol=True)
    finally:
        os.chdir(cwd)
    return os.path.join(workdir, f'{name}_outputs')


def strip_paths(folder, workdir):
    prefixes = sorted({workdir, os.path.realpath(workdir)}, key=len, reverse=True)
    for file_name in os.listdir(folder):
        if file_name.endswith(TEXT_SUFFIXES):
            path = os.path.join(folder, file_name)
            with open(path) as f:
                text = f.read()
            for prefix in prefixes:
                text = text.replace(prefix, '.')
            with open(path, 'w') as f:
                f.write(text)


def main():
    for name, pdb_name in EXAMPLES.items():
        with tempfile.TemporaryDirectory() as workdir:
            out = run_example(name, pdb_name, workdir)
            strip_paths(out, workdir)
            target = os.path.join(REFERENCE_DIR, f'{name}_outputs')
            shutil.rmtree(target, ignore_errors=True)
            shutil.copytree(out, target)
            print(f"{name}: wrote {sorted(os.listdir(target))}")


if __name__ == '__main__':
    main()
