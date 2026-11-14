# 1. Overview
RAPA (Rotamer And Protonation Assignment) is a command line tool that aids in protein preparation by assigning various protonation and rotamer states of unassigned amino-acids of the given structure and outputs all the models that are energetically accessible. Details of the algorithm can be found in the following manuscript:\
[DOI: 10.1021/acs.jcim.5c00859](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00859)

Copyright (c) 2025 Ventus Therapeutics U.S., Inc.

# 2. Algorithm
The code follows the following flowchart (fig1 in the [paper](https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00859)):
![Flowchart](flowchart.jpeg)

# 3. Getting started

Minimum of Python 3.9 is required.

## Installation
1. Download or git clone using: https://github.com/Ventus-Therapeutics/rapa.git
2. Run:  ``` pip install -r requirement.txt ``` or ``` conda install --yes --file requirement.txt ``` to install the dependencies.


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
$<path to rapa folder>/src/rapa.py -pID [input_pdb_name.pdb] -o [out_name]   
$<path to rapa folder>/src/rapa.py -pID 1xl2_noH.pdb -o 1xl2_out  
```
##### Notes
1. The provided PDBs in the example folder have **no gaps, no steric clashes, no hydrogen atoms and has all residue names in AMBER format**, hence meeting the conditions to run RAPA on it.

#### Command-Line Arguments  

##### `-pID` or `--protein_id` (Required)  
- **Type:** String  
- **Description:** Specifies the input PDB file name. You can specify the file with the file extension or just the file prefix
- **Example:**  
  ```bash 
  $cd <path to rapa folder>/examples/
  $<path to rapa folder>/src/rapa.py -p 1bcd_noH
  or
  $<path to rapa folder>/src/rapa.py -p 1bcd_noH.pdb
  ```

##### `-o` or `--out_name` (Optional)  
- **Type:** String  
- **Default:** `[pID]`  
- **Description:** Defines the prefix for the output PDB files and folders. If not provided, it will use pID of the input PDB. The output PDBs will be {prefix}_#.pdb in {prefix}_output folder
- **Example:**  
  ```bash  
  $cd <path to rapa folder>/examples/
  $<path to rapa folder>/src/rapa.py -pID 1bcd_noH.pdb -o test 
  ```  
  This will generate output files in the folder `test_outputs`: `test_0.pdb`, `test_1.pdb`, etc.  
  If omitted: 
  ```bash  
  $<path to rapa folder>/src/rapa.py -pID 1bcd_noH.pdb'  
  ```
  This will generate output files in the folder `1bcd_noH_outputs`: `1bcd_noH_0.pdb`, `1bcd_noH_1.pdb`, etc
  

##### `--ligand` (Optional)
- **Type:** String (comma-separated residue names)
- **Description:** Includes the named HETATM residues (ligands) as fixed hydrogen-bond partners when RAPA scores states. Ligand atoms are typed from their element and the hydrogens in the input PDB: O bonded to H is donor/acceptor, O without H is an acceptor, N bonded to H is a donor, N without H is an acceptor. Ligand hydrogens in the input are kept (in the output too). Ligand lone pairs are assumed to point straight at their partner (best-case geometry). **Provide the ligand with explicit hydrogens:** without them donors cannot be identified, so only the ligand's O atoms are used (as acceptors) and a warning is written to the `.info` file.
- **Example:**
  ```bash
  $<path to rapa folder>/src/rapa.py -pID 1xl2_noH.pdb --ligand 189
  ```

##### `--fix` (Optional)
- **Type:** String, `CHAIN:RESNUM:STATE` (comma-separated or repeated)
- **Description:** Pins residues to a state instead of letting RAPA assign them; RAPA then assigns everything else around them. States: `HID`, `HIE`, `HIP`, `ASN`, `GLN`, `ASP`, `GLU`, plus `HIDR`, `HIER`, `HIPR`, `ASNR`, `GLNR` for a side chain flipped 180° relative to the input. ASH/GLH cannot be pinned. When one ASP/GLU of an interacting pair is pinned, its partner is kept deprotonated too.
- **Example:**
  ```bash
  $<path to rapa folder>/src/rapa.py -pID 1bcd_noH.pdb --fix A:64:HIP,A:67:ASNR
  ```

##### `--around` (Optional, requires `--ligand`)
- **Type:** Float (Å)
- **Description:** Only assigns HIS/ASN/GLN/ASP/GLU residues with any atom within this distance of the `--ligand` residues. All other residues keep their input state (HID/HIE/HIP/ASH/GLH names are kept; an unspecified HIS becomes HIE). Much faster on large systems. Note that kept ASH/GLH residues outside the region do not contribute hydrogen-bond energy to residues inside it.
- **Example:**
  ```bash
  $<path to rapa folder>/src/rapa.py -pID 1xl2_noH.pdb --ligand 189 --around 8
  ```

##### `--max_outputs` (Optional)
- **Type:** Integer
- **Default:** `0` (no limit)
- **Description:** Writes at most this many output PDB files. Degenerate states are explored lowest-energy first, so the first N combinations are kept. If the cap is reached, a warning is written to the `.info` file.

##### Listing Available Options  
Use the `-h` flag to view all available command-line options:  
```bash  
$<path to rapa folder>/src/rapa.py -h  
```

---
# 4. Files and Folders Generated by RAPA  
All output files will be in `[input_pdb_prefix or specified_output_prefix]_outputs` folder
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
