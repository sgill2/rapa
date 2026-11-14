#!/usr/bin/env python
"""
Main file to run rapa

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


"""

import argparse
import contextlib
import sys, os
from dataclasses import dataclass
import timeit
from datetime import datetime

import setup_protein as stp
import hydrogen_placement_sp2 as hsp2
import state_assignment as sa
import global_constants as gc
import ligand as lig
import residue_selection as rsel
import summary
from errors import RapaInputError
from misc import prettify_time, Logger, generate_pymol_script, clean_pdb_atom_id, put_pdb_record_back

def get_unknown_list_for_pml(structure):
    all_unknown_residues = [] # this will be [('A', 253, 'HIE'), ('B', 101, 'SER')...]]
    for unknown_res in stp.get_unknown_residue_list(structure):
        chainID = unknown_res.parent.id
        res_num = unknown_res.id[1]
        all_unknown_residues.append((chainID, res_num, unknown_res.resname))
    return all_unknown_residues


def is_hydrogen_atom_line(line):
    """
    Determine whether an ATOM/HETATM PDB line represents a hydrogen atom.

    Primary check uses the element symbol column (columns 77-78, 0-indexed
    76:78) per the PDB format spec. If that column is blank/absent (common in
    older or loosely-formatted PDBs), falls back to a heuristic on the atom
    name field (columns 13-16, 0-indexed 12:16).
    """
    element = line[76:78].strip() if len(line) >= 78 else ''
    if element:
        return element.upper() == 'H'

    # Fallback: infer from atom name when element column is missing.
    atom_name = line[12:16].strip() if len(line) >= 16 else ''
    if not atom_name:
        return False
    # Standard PDB atom-name convention: hydrogens either start with 'H'
    # (e.g. "H", "HA", "HB2") or, when the name is 4 characters and starts
    # with a digit due to right-justification rules, have 'H' as the next
    # character (e.g. "1HB2", "2HD1").
    if atom_name[0] == 'H':
        return True
    if atom_name[0].isdigit() and len(atom_name) > 1 and atom_name[1] == 'H':
        return True
    return False


def strip_hydrogens_from_pdb(in_path, out_folder, keep_resnames=()):
    """
    Read the PDB at in_path, remove all hydrogen ATOM/HETATM records, and
    write the result to {out_folder}/{basename}_stripH.pdb.
    Hydrogens of residues named in keep_resnames (the --ligand residues) are kept.

    This is an intermediate file: the caller is responsible for deleting it.
    Returns the full path of the written file.
    """
    if not os.path.exists(in_path):
        raise FileNotFoundError(f"Input PDB file not found: {in_path}")

    base = os.path.splitext(os.path.basename(in_path))[0]
    out_path = os.path.join(out_folder, f"{base}_stripH.pdb")

    num_removed = 0
    num_kept = 0
    with open(in_path, 'r') as fin, open(out_path, 'w') as fout:
        for line in fin:
            if (line.startswith(('ATOM', 'HETATM')) and is_hydrogen_atom_line(line)
                    and line[17:20].strip() not in keep_resnames):
                num_removed += 1
                continue
            if line.startswith(('ATOM', 'HETATM')):
                num_kept += 1
            fout.write(line)

    print(f"Removed {num_removed} hydrogen atom record(s) from {in_path}; "
          f"{num_kept} non-hydrogen atom record(s) written to {out_path}")

    return out_path


def parse_arguments(argv):
    """
    process arguments and parse them to global when needed
    """
    parser = argparse.ArgumentParser(description='RAPA')
    parser.add_argument('-pID', '--protein_id',dest='protID', type=str, required=True, metavar='',
                        help='Input protein ID name. It is a required argument of data type string.')

    parser.add_argument('-o', '--out_name',dest='out_prefix', type=str, default='', metavar='',
                        help='Prefix of output file name. It is an optional argument of data type string. If not provided, it will use the default value as "[pID]_out". ')

    parser.add_argument('-s', '--single_pdb_out',
                        help='This flag requests the program to provide a single output pdb file. It is done by setting the cut off that defines degeneracy of unknown residues to zero.',
                        action="store_true")

    parser.add_argument('-l', '--loud', dest='log_file',
                        help='Print out more details and save it to a log file, default is off',
                        action='store_true')

    parser.add_argument('-d', '--debug',
                        help='Print out even more details for debugging',
                        action="store_true")

    parser.add_argument('-k_hlp', '--keep_hlp',
                        help = 'Keep the PDB created with the added hydrogen and lone pair coordinates ({prefix}_HLPsp2.pdb), Default: off (delete file after)',
                        action="store_true")

    parser.add_argument('-k', '--keep', dest='keep_asp_glu',
                        help='Keep all ASP/GLU residues in their deprotonated form; never convert them to ASH/GLH. Default: off',
                        action="store_true")

    parser.add_argument('-hlp', '--HLPsp2_known',
                        help='A flag to indicate if the additional PDB (having suffix "_HLPsp2.pdb"), with SP2 hydrogen and lone pair coordinates exist. This can occur if the user requested the program to generate this additional PDB by providing "-k_hlp" as a flag in the previous run for the same PDB. Note, it is assumed RAPA for the given PDB is needed to run twice. During the first run the user provided the "-k_hlp" flag and generated the additional PDB (with suffix "_HLPsp2.pdb"). This was then used in the second run by providing "-hlp" in the second run. This two step approach is recommended for debugging as it stores an intermediate calculation and saves computational time.',
                        action="store_true")

    parser.add_argument('--ligand', dest='ligand', type=str, default='', metavar='',
                        help='Comma-separated residue name(s) of ligands (HETATM) to include as hydrogen-bond '
                             'partners, e.g. --ligand LIG or --ligand LIG,NAG. Ligand hydrogens in the input are '
                             'kept and used to identify donors; without them only O atoms are used (as acceptors).')

    parser.add_argument('--fix', dest='fix', action='append', default=[], metavar='',
                        help='Pin residues to a state instead of letting RAPA assign them: CHAIN:RESNUM:STATE, '
                             'comma-separated or repeated, e.g. --fix A:94:HIP,A:120:ASNR. States: HID, HIE, HIP, '
                             'ASN, GLN, ASP, GLU, and HIDR/HIER/HIPR/ASNR/GLNR for a flipped side chain.')

    parser.add_argument('--around', dest='around', type=float, default=None, metavar='',
                        help='Only assign residues within this distance (Angstrom) of the --ligand residues. '
                             'All other residues keep their input state (an unspecified HIS becomes HIE). '
                             'Requires --ligand.')

    parser.add_argument('--max_outputs', dest='max_outputs', type=int, default=0, metavar='',
                        help='Write at most this many output PDBs. Degenerate states are explored lowest-energy '
                             'first, so the first N combinations are kept. Default: 0 (no limit)')

    parser.add_argument('--pymol',
                        help='Write a PyMOL script (inspect.pml) that loads all output PDBs, colors the investigated '
                             'residues by how their state was decided (changed, flipped, low confidence, degenerate, '
                             'fixed, metal), labels them with state and energy gap, and stores a scene per residue',
                        action='store_true')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1.1 ')

    parser.add_argument('--output_dir', dest='output_dir', type=str, default=None, metavar='',
                        help='Folder for the outputs. Default: [out_name]_outputs in the current directory')

    args = parser.parse_args(argv)
    try:
        check_options(args)
    except RapaInputError as e:
        print(f"ERROR: {e}. Exiting...")
        sys.exit(1)
    return args


def check_options(args):
    """
    Validate and normalize the options shared by the command line and run(). Raises RapaInputError.
    Adds: input_pdb (path), protID (file name without extension), ligand_names, fixes.
    """
    # the input can be given with or without .pdb, and may include a directory
    path = args.protID
    if not os.path.isfile(path) and os.path.isfile(path + '.pdb'):
        path = path + '.pdb'
    if not os.path.isfile(path):
        raise RapaInputError(f"the input pdb file {path} does not exist")
    args.input_pdb = path
    args.protID = os.path.splitext(os.path.basename(path))[0]

    if isinstance(args.ligand, str):
        args.ligand_names = {name.strip() for name in args.ligand.split(',') if name.strip()}
    else:
        args.ligand_names = {str(name).strip() for name in args.ligand}
    if args.around is not None:
        if not args.ligand_names:
            raise RapaInputError('--around requires --ligand')
        if args.around <= 0:
            raise RapaInputError('--around must be a positive distance')
    if args.max_outputs < 0:
        raise RapaInputError('--max_outputs must be 0 (no limit) or a positive number')
    args.fixes = rsel.parse_fix_args(args.fix)

    if not args.out_prefix:
        args.out_prefix = args.protID
    if not args.output_dir:
        args.output_dir = f"{args.out_prefix}_outputs"
    args.output_folder = args.output_dir


def configure_run(args):
    """Create the output folder and reset the module-level settings for a new run."""
    if os.path.isdir(args.output_folder):
        print(f"WARNING: specified output folder {args.output_folder} already exist, will overwrite")
    else:
        os.makedirs(args.output_folder)
    print(f"The calculation result will be saved to the folder: {args.output_folder}")

    args.log_name = os.path.join(args.output_folder, f"{args.out_prefix}.log")

    gc.debug = args.debug
    gc.log_file = args.log_file
    gc.out_folder = args.output_folder
    gc.out_info_file = os.path.join(args.output_folder, f"{args.out_prefix}.info")
    gc.max_outputs = args.max_outputs
    gc.outputs_truncated = False
    gc.ECutOff = 0 if args.single_pdb_out else 1
    summary.reset()


@dataclass
class RapaResult:
    """What a RAPA run produced. `residues` holds one summary row per residue per output (see summary.py)."""
    output_files: list
    output_folder: str
    info_file: str
    summary_csv: str
    summary_json: str
    residues: list
    degenerate_residues: list
    outputs_truncated: bool
    pymol_script: str = None
    log_file: str = None

    @property
    def changes(self):
        """Rows for residues whose state differs from the input."""
        return [row for row in self.residues if row['changed']]

    def residues_in(self, output_file):
        """Rows for one output file (full path or file name)."""
        name = os.path.basename(output_file)
        return [row for row in self.residues if row['output_file'] == name]


def _run(args, program_call):
    """Run RAPA for already-validated options. Returns a RapaResult."""
    starttime = timeit.default_timer()
    print(f'PROGRAM CALL: {program_call}')

    # first check that the input pdb file exists
    pdb_file = args.input_pdb
    if not os.path.isfile(pdb_file):
        raise RapaInputError(f"the input pdb file {pdb_file} does not exist")

    # Remove hydrogens from the input PDB before doing anything else.
    # RAPA requires hydrogen-free input; this makes that requirement automatic
    # instead of relying on the user to have pre-stripped hydrogens themselves.
    # The stripped copy is an intermediate file written to the output folder
    # (the user's input folder is left untouched) and deleted at the end.
    stripped_pdb_path = strip_hydrogens_from_pdb(pdb_file, args.output_folder, keep_resnames=args.ligand_names)

    # then clean up the (hydrogen-free) input pdb file atom ID
    cleaned_id_pdb_name, cleaned_id_pdb_full_path = clean_pdb_atom_id(stripped_pdb_path)
    os.remove(stripped_pdb_path)

    # set up structure, know the number of model and chain in the pdb
    cleaned_id_pdb_prefix = os.path.splitext(cleaned_id_pdb_name)[0]
    structure = stp.setup_structure(args.protID, outFolder=args.output_folder, fName=cleaned_id_pdb_prefix)
    chains = [item for sublist in structure for item in sublist]
    num_models = len(structure.child_list)
    num_chains = len(chains)

    if gc.log_file:
        print(f"The PDB ID being used: {args.protID} and start time is: "
              f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"\nThe structure has number of models {num_models}: {structure.child_list} and number of chains: {num_chains}: {chains}\n")

    with open(gc.out_info_file, "w") as fInfo:
        fInfo.write(f"PDB ID: {args.protID} \n")
        fInfo.flush()
        fInfo.write(f"Models: {num_models}: {structure.child_list}\n")
        fInfo.write(f"Chains: {num_chains}: {chains}\n\n")
        fInfo.flush()

    if gc.debug:
        print(f"PDB ID: {args.protID}")
        print(f"Models: {num_models}: {structure.child_list}")
        print(f"Chains: {num_chains}: {chains}\n")


    # check the user's residue selections against the input structure
    try:
        if args.ligand_names:
            lig.setup_ligands(structure, args.ligand_names)  # exits if a ligand name is not found
        rsel.validate_fixes(structure, args.fixes)
    except RapaInputError:
        os.remove(cleaned_id_pdb_full_path)
        raise
    out_of_region_keys = set()
    if args.around is not None:
        out_of_region_keys = rsel.get_out_of_region_keys(structure, args.ligand_names, args.around)

    # remember the input residue names for the change report
    summary.set_input_names(structure)

    # rename HID/HIE/HIP/ASH/GLH->HIS/HIS/HIS/ASP/GLU
    # residues outside the --around region (and not --fix'ed) keep their input names
    stp.change_HIDEP_ASH_GLH(structure, keep_keys=out_of_region_keys - set(args.fixes))

    # initialize the known/unknown residues(ASN/GLN/HIS, and mark if 2 ASPs/GLUs are in h-bonding distance to each other)
    stp.set_initial_known_residues_and_rotamers(structure)

    gc.xc_orig, gc.yc_orig, gc.zc_orig = stp.get_centroid(structure)
    # normalize so the new centroid=(0,0,0)
    stp.normalize_atom_coords(structure)

    ###########################Adding hydrogen and  connected to sp2 side chain atoms of unknown residue ##########################
    if not args.HLPsp2_known:
        lastSerial = list(structure.get_residues())[-1].child_list[-1].serial_number
        # This does not use neighboring residues
        lastSerial, _ = hsp2.place_lonepair_on_backbone(structure, lastSerial)
        lastSerial, _ = hsp2.add_sp2_sidechain_lonepairs(structure, lastSerial)
        lastSerial, _ = hsp2.placeHydrogens_backbone(structure, lastSerial)
        # Adding hydrogen to side chains:
        lastSerial, _ = hsp2.add_sp2_sidechain_hydrogens(structure, lastSerial)
        fName = os.path.join(args.output_folder, f'{args.protID}_HLPsp2.pdb')
        stp.write_to_PDB(structure, fName, removeHLP=False, removeHall=False, set_original_centroid=False)
        if gc.log_file:
            print("After adding the hydrogens....")
        protIDName = args.protID + "_HLPsp2"

    else:
        protIDName = args.protID + "_HLPsp2"

        if gc.log_file:
            print("Hydrogen and lone pairs were already present!")


    # setting up the new structure with hydrogen and lone pairs
    structure = stp.setup_structure(protIDName, outFolder=args.output_folder, fName=None)
    stp.set_initial_known_residues_and_rotamers(structure)

    # flagging the unknown ASP/GLUs (skipped with -k/--keep so ASP/GLU are never converted to ASH/GLH)
    over2ASPs = 0
    over2GLUs = 0
    if args.keep_asp_glu:
        if gc.log_file:
            print("Keeping all ASP/GLU residues as ASP/GLU (-k/--keep); no ASH/GLH will be generated.")
    else:
        unknownASP, unknownASP_atomInfo = stp.get_all_unknown_ASP_GLU(structure, resName='ASP')
        unknownGLU, unknownGLU_atomInfo = stp.get_all_unknown_ASP_GLU(structure, resName='GLU')
        if (unknownASP):
            structure, uniqueIDsASP, over2ASPs = stp.accomodate_for_ASPs_GLUs(structure, unknownASP, unknownASP_atomInfo)
        if (unknownGLU):
            structure, uniqueIDsGLU, over2GLUs = stp.accomodate_for_ASPs_GLUs(structure, unknownGLU, unknownGLU_atomInfo)

    # ligands take part in the energy evaluation as known hydrogen-bond partners
    selection_lines = []
    if args.ligand_names:
        selection_lines.append("Ligands used as hydrogen-bond partners:\n")
        selection_lines += lig.setup_ligands(structure, args.ligand_names)

    # pin --fix residues and residues outside the --around region
    if args.fixes or out_of_region_keys:
        selection_lines.append("Residues not assigned by RAPA:\n")
        selection_lines += rsel.apply_pins(structure, args.fixes, out_of_region_keys)

    if selection_lines:
        msg = '\n' + ''.join(selection_lines)
        with open(gc.out_info_file, "a") as fInfo:
            fInfo.write(msg)
        print(msg)

    # Flag the side chain residue hydrogen position for SER, THR, LYS, TYR is unknown.
    stp.set_initial_residue_side_chain_hydrogen_unknown(structure)

    # Get all unkown residues before start resolving structures
    all_unknown_res = get_unknown_list_for_pml(structure)

    # Print out unknown residues that RAPA is dealing with
    lines = ['\n'+"#"*100, f"\nTotal number of unknown residues: {len(all_unknown_res)}\n", "Unknown residues:\n"]
    for (chain, res_id, res_name) in all_unknown_res:
        lines.append(f" {res_name}-{res_id} of chain {chain}\n")
    msg = ''.join(lines)
    with open(gc.out_info_file, "a") as fInfo:
        fInfo.write(msg)
        fInfo.flush()
    print(msg)

    # first set the states of HIS/CYS around Zn
    structure = sa.set_states_around_metal(structure)
    generated_files, _, branched_residues = sa.resolve_residue_ambiguities_in_one_structure(structure, set_original_centroid=True,
                                                                         generated_files=None, pdb_file_num=None,
                                                                         outprefix=args.out_prefix, branched_residues=None)
    total_run_time = timeit.default_timer() - starttime
    # print out branching info
    if len(branched_residues) != 0:
        lines = ['#'*100+'\n', f'Residues w/ degenerate states: (energy of different states within '
                              f'{gc.ECutOff} kcal/mol)\n']
        for (chain, res_id, res_name) in list(set(branched_residues)):
            lines.append(f" {res_name}-{res_id} of chain {chain}\n")
        msg = ''.join(lines)
        with open(gc.out_info_file, "a") as fInfo:
            fInfo.write(msg)
            fInfo.flush()
        print(msg)

    # machine-readable per-residue summary and the human-readable change report
    settings = {'input_pdb': args.input_pdb, 'energy_cutoff_kcal_mol': gc.ECutOff,
                'ligands': sorted(args.ligand_names), 'around_angstrom': args.around,
                'fixed': {rsel.format_key(k): v for k, v in args.fixes.items()},
                'keep_asp_glu': args.keep_asp_glu, 'max_outputs': args.max_outputs or None,
                'outputs_truncated': gc.outputs_truncated}
    summary_csv, summary_json = summary.write_files(args.output_folder, args.out_prefix, settings)
    msg = summary.change_report() + f"Summary files: {summary_csv}, {summary_json}\n"
    with open(gc.out_info_file, "a") as fInfo:
        fInfo.write(msg)
    print(msg)

    if gc.outputs_truncated:
        msg = (f"WARNING: output capped at {gc.max_outputs} PDB file(s) by --max_outputs; "
               f"further degenerate combinations were not written.\n")
        with open(gc.out_info_file, "a") as fInfo:
            fInfo.write(msg)
        print(msg)

    pymol_script = None
    if args.pymol:
        print("Generating pymol script based on output PDB files")
        pymol_script = os.path.join(args.output_folder, 'inspect.pml')
        rows_by_file = dict(summary.get_outputs())
        generate_pymol_script([(f, rows_by_file.get(os.path.basename(f), [])) for f in generated_files],
                              pymol_script, ligand_names=args.ligand_names)

    generated_files_str = '\n'.join(generated_files)
    lines = ['#'*100+'\n',f'Done. Execution time is: {prettify_time(total_run_time)}\n',f"Number of files generated: {len(generated_files)} and files generated:\n{generated_files_str}\n"]
    msg = ''.join(lines)
    with open(gc.out_info_file, "a") as fInfo:
        fInfo.write(msg)
        fInfo.flush()
    print(msg)

    if (over2ASPs > 0 or over2GLUs > 0):
        if gc.log_file:
            print(f"WARNING: There are more than 2 ASPs/GLUs.\n\n over2ASPs: {over2ASPs}\n over2GLUs: {over2GLUs}. Human intervention is required. Look for WARNING in log_file file for details")

        with open(gc.out_info_file, "a") as fInfo:
            fInfo.write(
                f"WARNING: There are more than 2 ASPs/GLUs.\n\n over2ASPs: {over2ASPs}\n over2GLUs: {over2GLUs}. Human intervention is required. Look for WARNING in log_file file for details\n")
            fInfo.flush()

    if not args.keep_hlp:
        file_hlp = os.path.join(args.output_folder, f"{args.protID}_HLPsp2.pdb")
        if os.path.exists(file_hlp):
            os.remove(file_hlp)
        else:
            if gc.debug:
                print(f"File {file_hlp} does not exist.")
    put_pdb_record_back(cleaned_id_pdb_full_path, generated_files)
    os.remove(cleaned_id_pdb_full_path)

    print("************************************************************")
    print(f"RAPA exiting. Run for the given PDB: {args.protID} is completed")
    print("************************************************************")

    return RapaResult(output_files=generated_files, output_folder=os.path.abspath(args.output_folder),
                      info_file=gc.out_info_file, summary_csv=summary_csv, summary_json=summary_json,
                      residues=[row for _, rows in summary.get_outputs() for row in rows],
                      degenerate_residues=sorted(set(branched_residues)),
                      outputs_truncated=gc.outputs_truncated, pymol_script=pymol_script,
                      log_file=args.log_name if args.log_file else None)


def _run_with_logging(args, program_call, verbose=True):
    """Set up the run, send output to the log file (-l) and/or hide it (verbose=False), and always restore stdout."""
    original_stdout = sys.stdout
    logger = None
    try:
        with contextlib.ExitStack() as stack:
            if not verbose:
                devnull = stack.enter_context(open(os.devnull, 'w'))
                stack.enter_context(contextlib.redirect_stdout(devnull))
            configure_run(args)
            if gc.log_file:
                # If the log file option is on, print out info and direct stdout to file as well
                logger = Logger(args.log_name)
                sys.stdout = logger
            try:
                return _run(args, program_call)
            finally:
                if logger is not None:
                    sys.stdout = logger.terminal
                    logger.log.close()
    finally:
        sys.stdout = original_stdout


def run(pdb, out_prefix=None, output_dir=None, ligands=(), fix=(), around=None, max_outputs=0,
        single_output=False, keep_asp_glu=False, pymol=False, log=False, debug=False, keep_hlp=False,
        hlp_known=False, verbose=False):
    """
    Run RAPA from Python.

        result = rapa.run('1xl2.pdb', ligands=['189'], around=8, fix=['A:25:ASP'])
        result.output_files      # paths of the output PDBs
        result.changes           # residues whose state differs from the input
        result.residues          # every residue row (see summary.py for the columns)

    Arguments mirror the command-line options:
        pdb           input PDB path (".pdb" may be omitted)
        out_prefix    -o; output_dir: --output_dir
        ligands       --ligand, a list of residue names
        fix           --fix, a list of "CHAIN:RESNUM:STATE" strings
        around        --around (Angstrom); max_outputs: --max_outputs
        single_output -s; keep_asp_glu: -k; pymol: --pymol; log: -l; debug: -d
        keep_hlp      -k_hlp; hlp_known: -hlp
        verbose       print progress to stdout (default off)

    Raises RapaInputError for invalid input or options. Returns a RapaResult.
    """
    if isinstance(fix, str):
        fix = [fix]
    args = argparse.Namespace(protID=str(pdb), out_prefix=out_prefix or '', output_dir=output_dir,
                              single_pdb_out=single_output, log_file=log, debug=debug, keep_hlp=keep_hlp,
                              keep_asp_glu=keep_asp_glu, HLPsp2_known=hlp_known, ligand=list(ligands),
                              fix=list(fix), around=around, max_outputs=max_outputs, pymol=pymol)
    check_options(args)
    options = {k: v for k, v in dict(out_prefix=out_prefix, output_dir=output_dir, ligands=list(ligands),
                                       fix=list(fix), around=around, max_outputs=max_outputs,
                                       single_output=single_output, keep_asp_glu=keep_asp_glu).items() if v}
    return _run_with_logging(args, f"rapa.run({str(pdb)!r}, {options})", verbose=verbose)


def main(argv):
    args = parse_arguments(argv)
    try:
        _run_with_logging(args, 'python {}'.format(' '.join(sys.argv)))
    except RapaInputError as e:
        print(f"ERROR: {e}. Exiting...")
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
