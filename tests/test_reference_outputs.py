"""
Regression test: running the examples must reproduce examples/solutions_to_compare/.

Compares, for 1bcd and 1xl2:
  - the set of output PDB files
  - every atom record: same atoms in the same order, same residue names, coordinates within COORD_TOL
  - the per-residue summary: same states, flips and decision sources, energies within ENERGY_TOL

Tolerances allow for floating-point differences between platforms and numpy/Biopython versions.
If a change to the results is intended, regenerate the references with tests/regenerate_reference_outputs.py
and review the diff.

Run with:   python tests/test_reference_outputs.py     (or: pytest tests/)
Takes about 15 seconds.
"""

import csv
import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))  # so the tests also run without installing

import rapa  # noqa: E402

REFERENCE_DIR = os.path.join(REPO, 'examples', 'solutions_to_compare')
EXAMPLES = {'1bcd': '1bcd_noH.pdb', '1xl2': '1xl2_noH.pdb'}
COORD_TOL = 0.01    # Angstrom
ENERGY_TOL = 0.02   # kcal/mol
CATEGORICAL = ['output_file', 'chain', 'resnum', 'icode', 'input_name', 'output_name', 'flipped', 'changed', 'source']
NUMERIC = ['energy', 'gap_to_next_best']


def read_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atoms.append(((line[12:16].strip(), line[17:20], line[21], line[22:27].strip()),
                              tuple(float(line[i:i + 8]) for i in (30, 38, 46))))
    return atoms


def read_summary(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def parse_state_energies(text):
    return {k: float(v) for k, v in (item.split(':') for item in text.split(';'))} if text else {}


class TestReferenceOutputs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.results = {name: rapa.run(os.path.join(REPO, 'examples', pdb), out_prefix=name,
                                      output_dir=os.path.join(cls.tmp.name, f'{name}_outputs'))
                       for name, pdb in EXAMPLES.items()}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def reference(self, name, file_name):
        return os.path.join(REFERENCE_DIR, f'{name}_outputs', file_name)

    def test_same_output_files(self):
        for name, result in self.results.items():
            with self.subTest(example=name):
                expected = sorted(f for f in os.listdir(os.path.join(REFERENCE_DIR, f'{name}_outputs'))
                                  if f.endswith('.pdb'))
                self.assertEqual(sorted(os.path.basename(f) for f in result.output_files), expected)

    def test_atoms_and_coordinates_match(self):
        for name, result in self.results.items():
            for path in result.output_files:
                file_name = os.path.basename(path)
                with self.subTest(example=name, output=file_name):
                    new, ref = read_atoms(path), read_atoms(self.reference(name, file_name))
                    self.assertEqual([a[0] for a in new], [a[0] for a in ref],
                                     'atom list or residue names differ from the reference')
                    worst = max((max(abs(x - y) for x, y in zip(a[1], b[1])), a[0]) for a, b in zip(new, ref))
                    self.assertLess(worst[0], COORD_TOL, f'largest coordinate difference at {worst[1]}')

    def test_residue_states_match(self):
        for name in self.results:
            with self.subTest(example=name):
                new = read_summary(self.results[name].summary_csv)
                ref = read_summary(self.reference(name, f'{name}_summary.csv'))
                self.assertEqual([[r[c] for c in CATEGORICAL] for r in new],
                                 [[r[c] for c in CATEGORICAL] for r in ref])
                for a, b in zip(new, ref):
                    label = f"{a['output_file']} {a['chain']}:{a['resnum']}"
                    for column in NUMERIC:
                        if a[column] == '' or b[column] == '':
                            self.assertEqual(a[column], b[column], f'{label} {column}')
                        else:
                            self.assertAlmostEqual(float(a[column]), float(b[column]), delta=ENERGY_TOL,
                                                   msg=f'{label} {column}')
                    ea, eb = parse_state_energies(a['state_energies']), parse_state_energies(b['state_energies'])
                    self.assertEqual(sorted(ea), sorted(eb), f'{label} candidate states')
                    for state in ea:
                        self.assertAlmostEqual(ea[state], eb[state], delta=ENERGY_TOL, msg=f'{label} {state}')


if __name__ == '__main__':
    unittest.main(verbosity=2)
