"""
Machine-readable summary and change report for a RAPA run.

For every output PDB, one row per residue whose state RAPA can change (HIS/ASN/GLN/ASP/GLU families, plus
any residue RAPA touched, e.g. CYS->CYM next to a metal):

    output_file, chain, resnum, icode, input_name, output_name, flipped, changed, source,
    energy, gap_to_next_best, state_energies

source is one of:
    assigned            lowest-energy state, clearly better than the alternatives
    degenerate_branch   several states within the energy cutoff; this output is one of the branches
    no_polar_neighbors  no polar atoms nearby (HIS defaults to HIE, ASN/GLN kept as in the input)
    no_known_neighbors  only other unresolved residues nearby; kept as in the input
    metal               HIS/CYS coordinating a metal (HIM/CYM)
    fixed               pinned with --fix
    outside_region      outside the --around region; kept as in the input
    pair_partner_fixed  ASP/GLU whose pair partner was fixed or outside the --around region
    default             not investigated (e.g. an ASP/GLU with no ASP/GLU partner, or -k/--keep)

energy is the H-bond energy (kcal/mol) of the chosen state; gap_to_next_best is the energy of the best other
state minus the chosen one (larger = more confident). 0 means a tie that RAPA settled with its default rules;
negative means this state was not the lowest-energy one (a degenerate branch, or HIP, which is chosen when both
nitrogens are individually favorable). HIP/HIPR energies are listed in state_energies but not used for the gap.

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Released under the MIT License (see LICENSE).
"""

import csv
import json
import os

FAMILY = {'HIS': 'HIS', 'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HIM': 'HIS',
          'ASN': 'ASN', 'GLN': 'GLN',
          'ASP': 'ASP', 'ASH': 'ASP', 'GLU': 'GLU', 'GLH': 'GLU'}

COLUMNS = ['output_file', 'chain', 'resnum', 'icode', 'input_name', 'output_name', 'flipped', 'changed',
           'source', 'energy', 'gap_to_next_best', 'state_energies']

_input_names = {}   # (chain, resnum, icode) -> residue name in the input PDB
_outputs = []       # [(output file name, [row dicts])]


def reset():
    _input_names.clear()
    _outputs.clear()


def _key(residue):
    return (residue.get_parent().id, residue.id[1], residue.id[2].strip())


def set_input_names(structure):
    """Record residue names as they are in the input PDB (call before RAPA renames anything)."""
    for res in structure.get_residues():
        if res.id[0] == ' ':
            _input_names[_key(res)] = res.resname


def _round(value):
    return None if value is None else round(float(value), 3)


def _state_label(res):
    """Name of the residue's state as used in RAPA's energy tables, e.g. HIER or ASNR for a flipped side chain."""
    if res.resname in ('HID', 'HIE', 'HIP', 'ASN', 'GLN') and getattr(res, 'isRotamer', 0) == 1:
        return res.resname + 'R'
    return res.resname


def _row(file_name, res):
    key = _key(res)
    input_name = _input_names.get(key, res.resname)
    flipped = getattr(res, 'isRotamer', 0) == 1
    source = getattr(res, 'rapa_source', 'default')
    state_energies = getattr(res, 'rapa_state_energies', None) if source in ('assigned', 'degenerate_branch') else None
    energy = getattr(res, 'rapa_energy', None) if state_energies is not None else None

    gap = None
    if state_energies and energy is not None:
        label = getattr(res, 'rapa_state', None) or _state_label(res)
        # HIP is chosen by its own rule (both nitrogens favorable), not by total energy, so it is not a competitor
        others = [e for name, e in state_energies.items() if name != label and name not in ('HIP', 'HIPR')]
        if others:
            gap = min(others) - energy

    return {
        'output_file': file_name,
        'chain': key[0],
        'resnum': key[1],
        'icode': key[2],
        'input_name': input_name,
        'output_name': res.resname,
        'flipped': flipped,
        'changed': flipped or res.resname != input_name,
        'source': source,
        'energy': _round(energy),
        'gap_to_next_best': _round(gap),
        'state_energies': {k: _round(v) for k, v in sorted(state_energies.items())} if state_energies else None,
    }


def record_output(pdb_path, structure):
    """Collect the residue rows for one output structure (call before the structure is written)."""
    rows = []
    for res in structure.get_residues():
        if res.id[0] != ' ':
            continue
        input_name = _input_names.get(_key(res), res.resname)
        if res.resname in FAMILY or input_name in FAMILY or hasattr(res, 'rapa_source'):
            rows.append(_row(os.path.basename(pdb_path), res))
    _outputs.append((os.path.basename(pdb_path), rows))


def write_files(out_folder, prefix, settings):
    """Write {prefix}_summary.csv and {prefix}_summary.json. Returns the two paths."""
    csv_path = os.path.join(out_folder, f"{prefix}_summary.csv")
    json_path = os.path.join(out_folder, f"{prefix}_summary.json")

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        for _, rows in _outputs:
            for row in rows:
                flat = dict(row)
                if flat['state_energies']:
                    flat['state_energies'] = ';'.join(f"{k}:{v}" for k, v in flat['state_energies'].items())
                writer.writerow(flat)

    with open(json_path, 'w') as f:
        json.dump({'settings': settings,
                   'outputs': [{'file': name, 'residues': rows} for name, rows in _outputs]}, f, indent=2)

    return csv_path, json_path


def _describe(row):
    if row['output_name'] != row['input_name'] and row['flipped']:
        change = f"{row['input_name']} -> {row['output_name']}, flipped"
    elif row['flipped']:
        change = f"{row['output_name']} flipped"
    else:
        change = f"{row['input_name']} -> {row['output_name']}"
    detail = row['source']
    if row['energy'] is not None:
        detail += f", E={row['energy']:.2f}"
    if row['gap_to_next_best'] is not None:
        detail += f", next best +{row['gap_to_next_best']:.2f}" if row['gap_to_next_best'] >= 0 \
            else f", next best {row['gap_to_next_best']:.2f}"
    return f"{row['chain']}:{row['resnum']}{row['icode']} {change} ({detail})"


def change_report():
    """
    Human-readable list of residues whose state differs from the input: changes common to every output first,
    then the changes specific to each output.
    """
    if not _outputs:
        return ''

    def signature(row):
        return (row['chain'], row['resnum'], row['icode'], row['output_name'], row['flipped'])

    changed_per_file = [{signature(r): r for r in rows if r['changed']} for _, rows in _outputs]
    common = set.intersection(*(set(c) for c in changed_per_file))

    lines = ['#' * 100 + '\n', 'Changes relative to the input PDB (flipped = side chain rotated 180 degrees):\n']
    if len(_outputs) > 1:
        lines.append(f"In all {len(_outputs)} outputs:\n")
    common_rows = sorted((changed_per_file[0][s] for s in common), key=lambda r: (r['chain'], r['resnum'], r['icode']))
    lines += [f" {_describe(r)}\n" for r in common_rows] or [" none\n"]

    if len(_outputs) > 1:
        for (name, _), changed in zip(_outputs, changed_per_file):
            specific = sorted((r for s, r in changed.items() if s not in common),
                              key=lambda r: (r['chain'], r['resnum'], r['icode']))
            # residues that differ between outputs but match the input in this one
            lines.append(f"Only in {name}:\n")
            lines += [f" {_describe(r)}\n" for r in specific] or [" none (other differences keep the input state)\n"]

    num_changed = len(common_rows)
    lines.append(f"{num_changed} residue(s) changed in every output; per-residue details for all outputs are in "
                 f"the _summary.csv/_summary.json files\n")
    return ''.join(lines)
