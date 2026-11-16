"""
Regression test for --ligand on streptavidin-biotin (examples/1stp_btnH.pdb).

Biotin's ureido O3 accepts an H-bond from Asn23 ND2. With the ligand included, RAPA should:
  - type biotin's atoms correctly (ureido NH donors, carbonyl and carboxylate O acceptors)
  - keep Asn23 in its crystal orientation (ND2 toward O3)
  - score Asn23's kept orientation clearly lower, and its flipped orientation clearly higher, than without the ligand
  - give the same Asn23 result with --around 8
  - not list residues outside the --around region (only given defaults) as changes in the .info report

Run with:   python tests/test_ligand_1stp.py     (or: pytest tests/)
Takes a few seconds.
"""

import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))  # so the tests also run without installing

import rapa  # noqa: E402
from rapa import ligand as lig  # noqa: E402
from rapa.errors import RapaInputError  # noqa: E402
from Bio.PDB import PDBParser  # noqa: E402

PDB = os.path.join(REPO, 'examples', '1stp_btnH.pdb')


def asn23(result):
    rows = [r for r in result.residues if r['chain'] == 'A' and r['resnum'] == 23]
    assert len(rows) == 1, rows
    return rows[0]


class TestLigand1STP(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        out = cls.tmp.name
        cls.no_ligand = rapa.run(PDB, output_dir=os.path.join(out, 'nolig'))
        cls.with_ligand = rapa.run(PDB, ligands=['BTN'], output_dir=os.path.join(out, 'lig'))
        cls.around = rapa.run(PDB, ligands=['BTN'], around=8, output_dir=os.path.join(out, 'around'))
        # His87 and His127 are ~5.9 A from biotin: outside a 5 A region, so they only get the HIE default
        cls.around5 = rapa.run(PDB, ligands=['BTN'], around=5, output_dir=os.path.join(out, 'around5'))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_biotin_atom_typing(self):
        structure = PDBParser(QUIET=True).get_structure('x', PDB)
        btn = [r for r in structure.get_residues() if r.resname == 'BTN'][0]
        summary = lig.type_ligand_residue(btn)
        self.assertTrue(summary['has_H'])
        self.assertEqual(sorted(summary['donor']), ['N1', 'N2'])
        self.assertEqual(sorted(summary['acceptor']), ['O11', 'O12', 'O3'])
        self.assertEqual(summary['both'], [])

    def test_asn23_keeps_crystal_orientation(self):
        for result in (self.no_ligand, self.with_ligand, self.around):
            row = asn23(result)
            self.assertEqual(row['output_name'], 'ASN')
            self.assertFalse(row['flipped'])
            self.assertEqual(row['source'], 'assigned')

        # in the output, ND2 (the donor) still faces biotin O3
        structure = PDBParser(QUIET=True).get_structure('x', self.with_ligand.output_files[0])
        chain = structure[0]['A']
        btn = [r for r in chain if r.resname == 'BTN'][0]
        self.assertLess(chain[23]['ND2'] - btn['O3'], 3.2)
        self.assertGreater(chain[23]['OD1'] - btn['O3'], 3.8)

    def test_ligand_hbond_changes_asn23_energies(self):
        before = asn23(self.no_ligand)['state_energies']
        after = asn23(self.with_ligand)['state_energies']
        # the ND2-H...O3 H-bond lowers the kept orientation
        self.assertLess(after['ASN'], before['ASN'] - 2.0)
        # flipping would put OD1 against O3 (acceptor-acceptor repulsion)
        self.assertGreater(after['ASNR'], before['ASNR'] + 2.0)

    def test_around_matches_full_run(self):
        self.assertEqual(asn23(self.around)['state_energies'], asn23(self.with_ligand)['state_energies'])
        self.assertLess(sum(1 for r in self.around.residues if r['source'] != 'outside_region'),
                        sum(1 for r in self.with_ligand.residues))

    def test_change_report_skips_outside_region_defaults(self):
        defaulted = sorted(r['resnum'] for r in self.around5.residues
                           if r['source'] == 'outside_region' and r['changed'])
        self.assertEqual(defaulted, [87, 127])
        self.assertEqual(asn23(self.around5)['source'], 'assigned')
        with open(self.around5.info_file) as f:
            info = f.read()
        report = info[info.index('Changes relative to the input'):]
        self.assertNotIn('(outside_region)', report)
        self.assertIn(f"Not listed: {len(defaulted)} residue(s) outside the --around region", report)

    def test_unknown_ligand_name_is_an_input_error(self):
        with tempfile.TemporaryDirectory() as out:
            with self.assertRaises(RapaInputError):
                rapa.run(PDB, ligands=['XYZ'], output_dir=out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
