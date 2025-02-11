# TopMPI
A computational tool specifically designed for the identification of multiplexed TD-MS/MS spectra
## System Requirements
- **Python** (version 3.7 or later)  
  - Download from [python.org](https://www.python.org/downloads/)
- **pip** (Python package manager)  
  - It comes with Python by default. If not, install it using:
    ```sh
    python -m ensurepip --default-pip
    ```
- **virtualenv** (Recommended for isolated environments)  
  - Install it using:
    ```sh
    pip install virtualenv
    ```
- **TopPIC Suite** (version 1.7.6 or later)
    - Install instructions can be found [here](https://github.com/toppic-suite/toppic-suite)

## Installation Steps

### 1. Clone the Repository

If you haven’t already, clone the repository from GitHub:

```sh
git clone https://github.com/todanielwang/TopMPI.git
cd TopMPI
```

### 2. Create and Activate a Virtual Environment (Optional but Strongly Recommended)

Create a virtual environment:

```sh
python3 -m venv venv
```

Activate it:

- **Windows (cmd/PowerShell):**
  ```sh
  venv\Scripts\activate
  ```
- **Mac/Linux:**
  ```sh
  source venv/bin/activate
  ```

### 3. Install Dependencies

Install all required Python packages using:

```sh
pip install -r requirements.txt
```

## Uninstallation

To remove the virtual environment and dependencies, simply delete the `venv` folder:

```sh
rm -rf venv  # Mac/Linux
rmdir /s /q venv  # Windows
```
## Manual
### 1. Input
* A protein database file in the FASTA format
* A mass spectrum data file in the msalign format
* A text file containing LC-MS feature information (optional)
* A text file of fixed PTMs (optional)
* A text file of variable PTMs (optional)
* A text file of PTMs for the characterization of unexpected mass shifts (optional)
### 2. Output
TopMPI outputs 4 msalign files, 2 sets of TopPIC outputs, and 6 tab separated value (TSV) files. All outputs will be placed into one folder under name TopMPI. 

Format of TopPIC inputs and outputs can be found [here](https://www.toppic.org/software/toppic/manual.html). 

For example, when the input data file to TopMPI is spectra_ms2.msalign, the output includes:

* **spectra_TopMPI/First_ms2.msalign**: a msalign file with all MS/MS spectra in spectra_ms2.msalign using the first precursor as their precursor. 
* **spectra_TopMPI/Second_ms2.msalign**: a msalign file with only multiplexed MS/MS spectra in spectra_ms2.msalign using the second precursor as their precursor. 
* **spectra_TopMPI/Primary_ms2.msalign**: a msalign file with all MS/MS spectra in spectra_ms2.msalign using the TopMPI selected precursor as their precursor. 
* **spectra_TopMPI/Secondary_ms2.msalign**: a msalign file with only multiplexed MS/MS spectra in spectra_ms2.msalign with their primary matched experimental peaks removed, using the TopMPI selected secondary precursor as their precursor. 
* **spectra_TopMPI/First_\***: The output of TopPIC searched using First_ms2.msalign while no E-value/FDR filter were applied to the output results. 
* **spectra_TopMPI/Second_\***: The output of TopPIC searched using Second_ms2.msalign while no E-value/FDR filter were applied to the output results. 
* **spectra_TopMPI/Primary_ms2_toppic_prsm_single.tsv**: a TSV file containing identified proteoform spectrum-matches (PrSMs) selected by TopMPI from output of First_ms2.msalign and Second_ms2.msalign with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
* **spectra_TopMPI/Primary_ms2_toppic_proteoform_single.tsv**: a TSV file containing identified proteoforms selected by TopMPI from output of First_ms2.msalign and Second_ms2.msalign with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
* **spectra_TopMPI/Secondary_ms2_toppic_prsm_single.tsv**: a TSV file containing identified proteoform spectrum-matches (PrSMs) selected by TopMPI from secondary spectra and their primary matched peaks removed with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
* **spectra_TopMPI/Secondary_ms2_toppic_proteoform_single.tsv**: a TSV file containing identified proteoforms selected by TopMPI from from secondary spectra and their primary matched peaks removed with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
* **spectra_TopMPI/Total_PrSMs.tsv**: a TSV file containing the combined list of primary PrSMs and secondary PrSMs filtered using an E-value or spectrum-level FDR cutoff.
* **spectra_TopMPI/Total_Proteoforms.tsv**: a TSV file containing the combined list of primary proteoforms and secondary proteoforms, refiltered using TopFD features and an E-value or proteoform-level FDR cutoff. 

To browse identified proteins, proteoforms, and PrSMs in e.g. First_ms2.msalign, use a chrome browser to open the file spectra_TopMPI/First_html/topmsv/index.html. **Google Chrome** is recommended (Firefox and Edge are not recommended).

* Note that the folder spectra_html/topfd **needs** to be copied to spectra_TopMPI/First_html/ in order to use this function. 

### 3. Command line usage
To run TopMPI, open a terminal window and run the following command. 
```
TopMPI.py toppic-executable database-file spectrum-file [TopPIC options] [TopMPI options]
```
For TopPIC options, please see [here](https://www.toppic.org/software/toppic/manual.html). 

#### TopMPI Options
-h [ --help ]

Print the help message.
```
--alpha <a number between 0 and 1>
```
The intensity ratio between the first and second precursor required for a spectrum to treated as multiplexed. TopMPI only reselect primary precursor for multiplexed spectra. Default value: 0.2
```
--beta <a number between 0 and 1>
```
The percentage of shared matched experimental peaks required for the identifications of the two precursors to be treated as inconsistent. Default value: 0.9
```
--delta <an integer>
```
The offset to calculate the number of normalized matched fragment masses (NMFMs) based on number of unknown mass shifts. Default value: 5
```
--gamma <an integer>
```
The number of normalized matched fragment masses (NMFMs) difference required to switch to the second precursor Default value: 4

#### Examples
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta using an TopPIC executable toppic with a feature file spectra_ms2.feassture (reported by TopFD). The user does not need to specify the feature file name. Like TopPIC, TopMPI will automatically obtain the feature file name from the spectrum file name spectra_ms2.msalign.
```
TopMPI.sh toppic proteins.fasta spectra_ms2.msalign
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta without feature files.

```
TopMPI.sh toppic proteins.fasta spectra_ms2.msalign -x
```

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file and the NMFMs offset set to 4
```
TopMPI.sh toppic proteins.fasta spectra_ms2.msalign --delta 4
```
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. In an identified proteoform, at most 2 mass shifts are allowed and the maximum allowed mass shift value is 10,000 Dalton. Furthermore, the NMFM difference required to switch precursor is set to 5. 
```
TopMPI.sh toppic proteins.fasta spectra_ms2.msalign -s 2 -M 10000 --gamma 5
```

<!-- Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. The error tolerance for precursor and fragment masses is 5 ppm.

toppic -e 5 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Use the target decoy approach to compute spectrum level and proteoform-level FDRs, filter identified proteoform spectrum-matches by a 5% spectrum level FDR, and filter identified proteoforms by a 5% proteoform-level FDR.

toppic -d -t FDR -v 0.05 -T FDR -V 0.05 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign with alternating CID, HCD, and ETD spectra against a protein database file proteins.fasta with a feature file. Combine alternating CID, HCD, and ETD spectra to increase proteoform coverage.

toppic -r 3 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. After proteoforms with unexpected mass shifts are identified, TopPIC matches the mass shifts to four common PTMs: acetylation, phosphorylation, oxidation and methylation, and uses an MIScore cutoff 0.1 to filter reported PTM sites. The modification file common_mods.txt can be found here.

toppic -B common_mods.txt -H 0.1 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file.Use 6 CPU threads to speed up the computation.

toppic -u 6 proteins.fasta spectra_ms2.msalign -->