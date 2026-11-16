# 1. Overview
RAPA (Rotamer And Protonation Assignment) is a command line tool that aids in protein preparation by assigning various protonation and rotamer states of unassigned amino-acids of the given structure and outputs all the models that are energetically accessible. Details of the algorithm can be found in the following manuscript:\
[DOI: 10.1021/acs.jcim.5c00859](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00859)

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.

# 2. Algorithm
The code follows the following flowchart (fig1 in the [paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00859)):
![Flowchart](flowchart.jpeg)

# 3. Getting started

Minimum of Python 3.10 is required.

## Installation
1. Download or git clone using: https://github.com/Ventus-Therapeutics/rapa.git
2. From the rapa folder, install RAPA and its dependencies (Biopython, NumPy):
   ```bash
   $pip install .            # or "pip install -e ." to edit the code in place
   ```
   This provides the `rapa` command and the `rapa` Python package.

   Without installing: install the dependencies (``` pip install -r requirement.txt ``` or ``` conda install --yes --file requirement.txt ```) and run `python -m rapa ...` with the `src` folder on `PYTHONPATH`:
   ```bash
   $PYTHONPATH=<path to rapa folder>/src python -m rapa -pID input.pdb
   ```

**Upgrading from 1.1:** the command is now `rapa` (or `python -m rapa`) instead of `<path to rapa folder>/src/rapa.py`; the options are unchanged.


## Usage
### Preparing Files for RAPA (not supported by this repository)
#### RAPA requires that the input PDB file meet the following criteria
1. **No gaps** – Protein structure must be contiguous. There must be no missing residues or missing atoms.  
2. **No steric clashes** – The structure must not contain steric clashes.  
3. **No hydrogen atoms** – The PDB file must not contain hydrogen atoms.  
4. **AMBER format** – All residue names must follow the AMBER format.  

---
### Running RAPA  

#### Basic Usage  
RAPA is executed with a PDB file as input. For instance, to run RAPA on the PDB file `1xl2_noH.pdb` ( **having no gaps, no steric clashes, no hydrogen atoms and has all residue names in AMBER format**) in the examples folder, the following steps can be executed:
``` bash
$cd <path to rapa folder>/examples/
$rapa -pID [input_pdb_name.pdb] -o [out_name]   
$rapa -pID 1xl2_noH.pdb -o 1xl2_out  
```
##### Notes
1. The provided PDBs in the example folder have **no gaps, no steric clashes, no hydrogen atoms and has all residue names in AMBER format**, hence meeting the conditions to run RAPA on it.

#### Command-Line Arguments  

##### `-pID` or `--protein_id` (Required)  
- **Type:** String  
- **Description:** Specifies the input PDB file. It can include a directory (e.g. `examples/1bcd_noH.pdb`) and the `.pdb` extension can be omitted. The input folder is never modified.
- **Example:**  
  ```bash 
  $cd <path to rapa folder>/examples/
  $rapa -p 1bcd_noH
  or
  $rapa -p 1bcd_noH.pdb
  ```

##### `-o` or `--out_name` (Optional)  
- **Type:** String  
- **Default:** `[pID]`  
- **Description:** Defines the prefix for the output PDB files and folders. If not provided, it will use pID of the input PDB. The output PDBs will be {prefix}_#.pdb in {prefix}_output folder
- **Example:**  
  ```bash  
  $cd <path to rapa folder>/examples/
  $rapa -pID 1bcd_noH.pdb -o test 
  ```  
  This will generate output files in the folder `test_outputs`: `test_0.pdb`, `test_1.pdb`, etc.  
  If omitted: 
  ```bash  
  $rapa -pID 1bcd_noH.pdb'  
  ```
  This will generate output files in the folder `1bcd_noH_outputs`: `1bcd_noH_0.pdb`, `1bcd_noH_1.pdb`, etc
  

##### `--ligand` (Optional)
- **Type:** String (comma-separated residue names)
- **Description:** Includes the named HETATM residues (ligands) as fixed hydrogen-bond partners when RAPA scores states. Ligand atoms are typed from their element and the hydrogens in the input PDB: O bonded to H is donor/acceptor, O without H is an acceptor, N bonded to H is a donor, N without H is an acceptor. Ligand hydrogens in the input are kept (in the output too). Ligand lone pairs are assumed to point straight at their partner (best-case geometry). **Provide the ligand with explicit hydrogens:** without them donors cannot be identified, so only the ligand's O atoms are used (as acceptors) and a warning is written to the `.info` file.
- **Example:**
  ```bash
  $rapa -pID 1xl2_noH.pdb --ligand 189
  ```

##### `--fix` (Optional)
- **Type:** String, `CHAIN:RESNUM:STATE` (comma-separated or repeated)
- **Description:** Pins residues to a state instead of letting RAPA assign them; RAPA then assigns everything else around them. States: `HID`, `HIE`, `HIP`, `ASN`, `GLN`, `ASP`, `GLU`, plus `HIDR`, `HIER`, `HIPR`, `ASNR`, `GLNR` for a side chain flipped 180° relative to the input. ASH/GLH cannot be pinned. When one ASP/GLU of an interacting pair is pinned, its partner is kept deprotonated too.
- **Example:**
  ```bash
  $rapa -pID 1bcd_noH.pdb --fix A:64:HIP,A:67:ASNR
  ```

##### `--around` (Optional, requires `--ligand`)
- **Type:** Float (Å)
- **Description:** Only assigns HIS/ASN/GLN/ASP/GLU residues with any atom within this distance of the `--ligand` residues. All other residues keep their input state (HID/HIE/HIP/ASH/GLH names are kept; an unspecified HIS becomes HIE). Much faster on large systems. Note that kept ASH/GLH residues outside the region do not contribute hydrogen-bond energy to residues inside it.
- **Example:**
  ```bash
  $rapa -pID 1xl2_noH.pdb --ligand 189 --around 8
  ```

##### `--max_outputs` (Optional)
- **Type:** Integer
- **Default:** `0` (no limit)
- **Description:** Writes at most this many output PDB files. Degenerate states are explored lowest-energy first, so the first N combinations are kept. If the cap is reached, a warning is written to the `.info` file.

##### `--output_dir` (Optional)
- **Type:** String
- **Default:** `[out_name]_outputs` in the current directory
- **Description:** Folder for all output files.

##### `--pymol` (Optional)
- **Description:** Writes `inspect.pml`, a PyMOL script that loads every output PDB and colors the investigated residues by how their state was decided: red = degenerate (differs between outputs), slate = fixed (`--fix`), grey = metal (HIM/CYM), yellow = low confidence (next best state within 1 kcal/mol, or a tie), magenta = flipped, orange = protonation changed, green lines = investigated but unchanged. Residues outside the `--around` region are not highlighted. Each residue is labelled with its state(s) and energy gap (e.g. `HID +0.5 | HID flip -0.5` for a degenerate residue), each category is a named selection (`rapa_low_confidence`, ...), each non-trivial residue has a scene, and `--ligand` residues are shown in cyan. Run it with `pymol inspect.pml`.

##### Listing Available Options  
Use the `-h` flag to view all available command-line options:  
```bash  
$rapa -h  
```

### Using RAPA from Python
After installing (or with `src` on the Python path), call `rapa.run`. Its arguments mirror the command-line options; progress output is off unless `verbose=True`. Invalid input raises `rapa.RapaInputError` instead of exiting.
```python
import rapa

result = rapa.run('examples/1xl2_noH.pdb', ligands=['189'], around=8,
                  fix=['A:25:ASP'], max_outputs=4, output_dir='1xl2_run', pymol=True)

result.output_files          # paths of the output PDBs
result.changes               # rows for residues whose state differs from the input
result.residues              # one row per residue per output (columns as in _summary.csv)
result.residues_in(result.output_files[0])
result.degenerate_residues   # [(chain, resnum, resname), ...] that were branched on
result.outputs_truncated     # True if max_outputs stopped further branching
result.info_file, result.summary_csv, result.summary_json, result.pymol_script
```
Other arguments: `out_prefix`, `single_output`, `keep_asp_glu`, `log`, `debug`, `keep_hlp`, `hlp_known`, `verbose`. Several runs can be made in one Python process. Runs are not thread-safe (RAPA keeps per-run settings at module level).

---
# 4. Files and Folders Generated by RAPA  
All output files will be in the `[input_pdb_prefix or specified_output_prefix]_outputs` folder (or `--output_dir`)
The following files can be found in that folder:

### 4.1 [input_pdb_prefix or specified_output_prefix].info
Information files with:  
- **Basic information:**  
  - PDB file name, ID, number of models, and chains.  
- **Residue and atom details:**
  - All unknown residues investigated by RAPA
  - Residues with degenerate states after calculation (if any)
  - Lists ASP/GLU residues investigated by RAPA.  
  - Displays unresolved ASP/GLU residues with their residue numbers, atom IDs, and distance details.  
      - If more than 2 ASP/GLU residues are within hydrogen bonding distance, RAPA does **not** resolve them. Manual intervention is required.  
- **Run summary:**  
  - Execution time.  
  - Number of output PDB files produced.  
  - File path of generated outputs.  

- **Changes relative to the input:** every residue whose state differs from the input PDB (renamed and/or side chain flipped), with how it was decided and its energy margin. With several outputs, changes common to all outputs are listed first, then those specific to each output. Residues outside the `--around` region that were only given a default (an unspecified HIS becomes HIE) are counted in one line instead of listed; they remain in the summary files with `source = outside_region`.

### 4.1.1 [prefix]_summary.csv and [prefix]_summary.json
One row per HIS/ASN/GLN/ASP/GLU-type residue (plus any other residue RAPA changed, e.g. CYS→CYM) per output PDB, for use in scripts and spreadsheets. The JSON file has the same rows grouped by output file, plus the run settings.

| Column | Meaning |
|---|---|
| `output_file`, `chain`, `resnum`, `icode` | which output and residue |
| `input_name`, `output_name` | residue name in the input and in this output |
| `flipped` | side chain rotated 180° relative to the input |
| `changed` | name changed or side chain flipped |
| `source` | how the state was decided: `assigned` (lowest energy), `degenerate_branch` (one of several states within the cutoff), `no_polar_neighbors`, `no_known_neighbors`, `metal`, `fixed` (`--fix`), `outside_region` (`--around`), `pair_partner_fixed`, `default` (not investigated) |
| `energy` | H-bond energy (kcal/mol) of the chosen state |
| `gap_to_next_best` | energy of the best alternative minus the chosen one; larger = more confident, 0 = tie settled by RAPA's default rules, negative = not the lowest-energy state (degenerate branch or HIP). HIP/HIPR are not counted as alternatives because RAPA chooses HIP with its own rule (both nitrogens favorable) |
| `state_energies` | energies of all candidate states (e.g. `HID:-2.05;HIE:3.68;...`) |

### 4.2 [input_pdb_prefix or specified_output_prefix]_#.pdb   
Final PDB files with states that RAPA picked based on the calculation.
### 4.3 [input_pdb_prefix or specified_output_prefix].log (Generated if `-l` flag is on)
Log file contains additional details of the calculation (more detailed than .info file). The messages printed out are saved here.
- looping and function call information
  -  Record when the recursive function is making a new function call
  -  Record how many iterations the loop has gone through
  -  Number of unknown states after each iteration
- Energy information
  - Details of energy during the evaluation of degenerate states

### 4.4 [input_pdb_prefix]_HLPsp2.pdb  (Generated if `--keep_hkp` flag is on)
pdb file with sp2 hydrogen and lone pair coordinates computed and added.  


---
# 5. Example Runs  
Example PDB files are available in the `examples` folder.
- `1bcd_noH.pdb`, `1xl2_noH.pdb`: proteins prepared for RAPA (reference outputs in `solutions_to_compare`).
- `1stp_btnH.pdb`: streptavidin with biotin, hydrogens added to biotin. Biotin's ureido O accepts an H-bond from Asn23, so this shows `--ligand` at work:
  ```bash
  $rapa -pID examples/1stp_btnH.pdb --ligand BTN --pymol
  ```
  With the ligand, Asn23's crystal orientation scores about 4 kcal/mol lower and its flipped orientation about 7 kcal/mol higher than without it (see `state_energies` in the summary files).

### Tests
```bash
$pip install ".[test]"
$pytest                     # about 15 seconds
```
- `tests/test_reference_outputs.py` runs the 1bcd and 1xl2 examples and checks that the output PDBs (atoms, residue names, coordinates within 0.01 Å) and every residue's state, flip, decision source and energies (within 0.02 kcal/mol) match `examples/solutions_to_compare/`.
- `tests/test_ligand_1stp.py` checks ligand typing, the Asn23 energies above, `--around`, and the change report.

Both also run without pytest, e.g. `python tests/test_reference_outputs.py`. GitHub Actions runs them on Python 3.10-3.13 for every push and pull request (`.github/workflows/tests.yml`).

If a change to RAPA's results is intended, regenerate the reference outputs and review the differences before committing:
```bash
$python tests/regenerate_reference_outputs.py
$git diff --stat examples/solutions_to_compare
```

---
# 6. Copyright and license
Code released under the [MIT License](https://github.com/Ventus-Therapeutics/rapa/blob/main/LICENSE)


MIT License
 

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


# 7. Authors:
This code is developed and released by Ventus Therapeutics U.S., Inc.
In conjunction with Tom Kurtzman's group at CUNY.
