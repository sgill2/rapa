"""
User control over which residues RAPA assigns.

--fix CHAIN:RESNUM:STATE   pin a residue to a given state instead of letting RAPA assign it
--around DIST              only assign residues within DIST Angstrom of the --ligand residues;
                           everything else keeps the state it has in the input PDB

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Released under the MIT License (see LICENSE).
"""

import Bio.PDB

import state_assignment as sa
import ligand as lig
import global_constants as gc

# states a residue can be pinned to with --fix, keyed by the residue's generic name.
# an "R" suffix means the side chain is flipped 180 degrees relative to the input
FIXABLE_STATES = {
    'HIS': ['HID', 'HIE', 'HIP', 'HIDR', 'HIER', 'HIPR'],
    'ASN': ['ASN', 'ASNR'],
    'GLN': ['GLN', 'GLNR'],
    'ASP': ['ASP'],
    'GLU': ['GLU'],
}
GENERIC_NAME = {'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HIS': 'HIS', 'ASN': 'ASN', 'GLN': 'GLN',
                'ASP': 'ASP', 'ASH': 'ASP', 'GLU': 'GLU', 'GLH': 'GLU'}


def residue_key(residue):
    """(chain, residue number, insertion code) - stable across re-parsing the structure."""
    return (residue.get_parent().id, residue.id[1], residue.id[2].strip())


def format_key(key):
    chain, num, icode = key
    return f"{chain}:{num}{icode}"


def parse_fix_args(fix_args):
    """
    Parse --fix values. Each value may hold several comma-separated entries like A:94:HIP or A:52A:ASNR.
    Returns {(chain, resnum, icode): STATE}.
    """
    fixes = {}
    for value in fix_args or []:
        for entry in value.split(','):
            entry = entry.strip()
            if not entry:
                continue
            parts = entry.split(':')
            if len(parts) != 3:
                raise SystemExit(f"ERROR: --fix entry '{entry}' should look like CHAIN:RESNUM:STATE, e.g. A:94:HIP")
            chain, num_str, state = parts[0].strip(), parts[1].strip(), parts[2].strip().upper()
            icode = ''
            if num_str and num_str[-1].isalpha():
                num_str, icode = num_str[:-1], num_str[-1]
            try:
                num = int(num_str)
            except ValueError:
                raise SystemExit(f"ERROR: --fix entry '{entry}' has an invalid residue number '{parts[1]}'")
            all_states = {s for states in FIXABLE_STATES.values() for s in states}
            if state not in all_states:
                raise SystemExit(f"ERROR: --fix entry '{entry}': unsupported state '{state}'. "
                                 f"Supported states: {', '.join(sorted(all_states))} (R = flipped side chain). "
                                 f"ASH/GLH cannot be pinned.")
            fixes[(chain, num, icode)] = state
    return fixes


def validate_fixes(structure, fixes):
    """Check that each fixed residue exists and that the state matches the residue type."""
    residues = {residue_key(res): res for res in structure.get_residues() if res.id[0] == ' '}
    for key, state in fixes.items():
        if key not in residues:
            raise SystemExit(f"ERROR: --fix residue {format_key(key)} not found in the input PDB")
        generic = GENERIC_NAME.get(residues[key].resname)
        if generic is None or state not in FIXABLE_STATES[generic]:
            raise SystemExit(f"ERROR: --fix {format_key(key)}:{state} does not match residue type "
                             f"{residues[key].resname}")


def get_residues_within(structure, center_residues, distance):
    """Keys of amino-acid residues with any atom within `distance` Angstrom of any atom of center_residues."""
    protein_atoms = [a for a in structure.get_atoms() if a.get_parent().id[0] == ' ']
    ns = Bio.PDB.NeighborSearch(protein_atoms)
    inside = set()
    for res in center_residues:
        for atom in res:
            for close in ns.search(atom.coord, distance):
                inside.add(residue_key(close.get_parent()))
    return inside


def get_out_of_region_keys(structure, ligand_names, distance):
    """
    Keys of the HIS/ASN/GLN/ASP/GLU-type residues that are outside --around DIST of the ligands.
    These keep their input state.
    """
    ligands = lig.find_ligand_residues(structure, ligand_names)
    inside = get_residues_within(structure, ligands, distance)
    return {residue_key(res) for res in structure.get_residues()
            if res.id[0] == ' ' and res.resname in GENERIC_NAME and residue_key(res) not in inside}


def pin_residue(structure, residue, state):
    """Set a residue to `state` (one of FIXABLE_STATES) and mark it known so RAPA will not change it."""
    model_id = residue.get_parent().get_parent().id
    chain_id = residue.get_parent().id
    res = structure[model_id][chain_id][residue.id]
    if res.resname == 'HIS':
        sa.setup_HIS(res, structure, state[:3])
        res = structure[model_id][chain_id][residue.id]
        if state.endswith('R'):
            sa.setup_rotamer_in_structure(res, structure)
    elif res.resname in ('ASN', 'GLN') and state.endswith('R'):
        sa.setup_rotamer_in_structure(res, structure)
    res.isKnown = 1


def release_orphaned_ASP_GLU(structure):
    """
    ASP/GLU are resolved in pairs. If one member of a pair has been pinned (known), the other cannot be
    resolved as a pair, so it is kept deprotonated as well. Returns the residues that were marked known.
    """
    released = []
    for res in list(structure.get_residues()):
        if res.resname in ('ASP', 'GLU') and getattr(res, 'isKnown', 1) == 0:
            partner = sa.get_ASP_GLU_pair(res, structure, gc.deltaD)
            if partner is None or partner.isKnown == 1:
                res.isKnown = 1
                released.append(res)
    return released


def apply_pins(structure, fixes, out_of_region_keys):
    """
    Pin user-fixed residues and residues outside the --around region.
    Returns report lines for the .info file.
    """
    lines = []
    residues = {residue_key(res): res for res in structure.get_residues() if res.id[0] == ' '}

    for key, state in fixes.items():
        pin_residue(structure, residues[key], state)
        lines.append(f" {format_key(key)} fixed to {state} (--fix)\n")

    kept = []
    for key in sorted(out_of_region_keys - set(fixes)):
        res = residues.get(key)
        if res is None:
            continue
        if res.resname == 'HIS':
            # HID/HIE/HIP in the input keep their name; only an unspecified HIS gets here
            pin_residue(structure, res, 'HIE')
            kept.append(f"{format_key(key)} HIS->HIE")
        elif res.resname in ('ASN', 'GLN', 'ASP', 'GLU'):
            pin_residue(structure, res, res.resname)
            kept.append(f"{format_key(key)} {res.resname}")
    if out_of_region_keys:
        lines.append(f" {len(out_of_region_keys)} HIS/ASN/GLN/ASP/GLU-type residue(s) outside the --around "
                     f"region keep their input state\n")
        if kept:
            lines.append(f"  of these, set as known: {', '.join(kept)}\n")

    for res in release_orphaned_ASP_GLU(structure):
        lines.append(f" {res.resname}-{res.id[1]} of chain {res.get_parent().id} kept as {res.resname}: "
                     f"its ASP/GLU pair partner is fixed or outside the --around region\n")
    return lines
