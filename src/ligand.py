"""
Ligand support for RAPA.

Ligand (HETATM) residues selected with --ligand take part in the hydrogen-bond
energy evaluation as fixed ("known") partners. Ligand atoms are typed as
donor / acceptor / both from their element and the explicit hydrogens present
in the input PDB:

    O bonded to H  -> both (e.g. hydroxyl)
    O without H    -> acceptor
    N bonded to H  -> donor
    N without H    -> acceptor
    anything else  -> no hydrogen-bonding behavior

Ligand hydrogens are kept from the input and used directly when a ligand atom
acts as a donor. Ligands carry no lone-pair atoms, so when a ligand atom acts
as an acceptor its lone pair is assumed to point straight at the partner atom
(the best-case orientation, gamma = 0).

If a ligand has no hydrogens at all, donors cannot be identified: its O atoms
are treated as acceptors, its N atoms are ignored, and a warning is reported.

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Released under the MIT License (see LICENSE).
"""

import numpy as np

import global_constants as gc


def is_ligand_residue(residue):
    """True if this residue was selected as a ligand with --ligand."""
    return getattr(residue, 'is_ligand', False)


def is_ligand_atom(atom):
    """True if this atom belongs to a residue selected as a ligand."""
    return is_ligand_residue(atom.get_parent())


def _is_hydrogen(atom):
    element = (atom.element or '').strip().upper()
    if element:
        return element == 'H'
    return atom.get_id().lstrip('0123456789').upper().startswith('H')


def get_ligand_hydrogens_bonded_to(atom):
    """Hydrogen atoms of the same ligand residue that are bonded to the given heavy atom."""
    return [h for h in atom.get_parent() if _is_hydrogen(h) and (h - atom) < gc.ligand_hbond_H_dist]


def pseudo_lonepair_toward(acceptor_atom, partner_atom):
    """
    Lone-pair coordinate for a ligand acceptor: 1 Angstrom from the acceptor, pointing straight at the partner.
    """
    direction = partner_atom.coord - acceptor_atom.coord
    norm = np.linalg.norm(direction)
    if norm < 1e-8:
        return acceptor_atom.coord.copy()
    return acceptor_atom.coord + direction / norm


def type_ligand_residue(residue):
    """
    Mark a residue as a ligand and assign a hydrogen-bonding behavior to each of its atoms.
    Returns a dict summarizing the typing: {'donor': [...], 'acceptor': [...], 'both': [...], 'has_H': bool}
    """
    residue.is_ligand = True
    hydrogens = [a for a in residue if _is_hydrogen(a)]
    has_H = len(hydrogens) > 0
    summary = {'donor': [], 'acceptor': [], 'both': [], 'has_H': has_H}

    for atom in residue:
        element = (atom.element or '').strip().upper()
        behavior = gc.doesNothing
        if element in ('N', 'O') and not _is_hydrogen(atom):
            bonded_H = [h for h in hydrogens if (h - atom) < gc.ligand_hbond_H_dist]
            if element == 'O':
                behavior = gc.both if bonded_H else gc.acceptor
            elif bonded_H:
                behavior = gc.donor
            elif has_H:
                # nitrogen with no hydrogen in a ligand that does have explicit hydrogens
                behavior = gc.acceptor
            # else: ligand has no hydrogens at all, N could be donor or acceptor, so ignore it

        atom.ligand_behavior = behavior
        if behavior is gc.donor:
            summary['donor'].append(atom.get_id())
        elif behavior is gc.acceptor:
            summary['acceptor'].append(atom.get_id())
        elif behavior is gc.both:
            summary['both'].append(atom.get_id())

    return summary


def find_ligand_residues(structure, ligand_names):
    """All residues in the structure whose residue name is in ligand_names."""
    return [res for res in structure.get_residues() if res.resname.strip() in ligand_names]


def setup_ligands(structure, ligand_names):
    """
    Find and type all ligand residues. Returns a list of report lines describing what was found.
    Exits with an error if a requested ligand name is not present in the structure.
    """
    ligands = find_ligand_residues(structure, ligand_names)
    found_names = {res.resname.strip() for res in ligands}
    missing = sorted(set(ligand_names) - found_names)
    if missing:
        raise SystemExit(f"ERROR: --ligand residue name(s) not found in the input PDB: {', '.join(missing)}")

    lines = []
    for res in ligands:
        summary = type_ligand_residue(res)
        label = f"{res.resname.strip()}-{res.id[1]} of chain {res.get_parent().id}"
        lines.append(f" {label}: donors {summary['donor'] or '-'}, acceptors {summary['acceptor'] or '-'}, "
                     f"donor/acceptor {summary['both'] or '-'}\n")
        if not summary['has_H']:
            lines.append(f"  WARNING: {label} has no hydrogens, so its donors cannot be identified. "
                         f"O atoms are treated as acceptors and N atoms are ignored. "
                         f"Provide the ligand with explicit hydrogens for a better result.\n")
    return lines
