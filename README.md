TopMPI
1. Input
An TopPIC executable, install instructions can be found [here](https://github.com/toppic-suite/toppic-suite)
A protein database file in the FASTA format
A mass spectrum data file in the msalign format
A text file containing LC-MS feature information (optional)
A text file of fixed PTMs (optional)
A text file of variable PTMs (optional)
A text file of PTMs for the characterization of unexpected mass shifts (optional)
2. Output
TopMPI outputs 4 msalign files, 3 sets of TopPICoutputs, and 3 tab separated value (TSV) files, all output into one folder under name TopMPI. Format of TopPIC inputs and outputs can be found [here](https://www.toppic.org/software/toppic/manual.html). For example, when the input data file is spectra_ms2.msalign, the output includes:

spectra_TopMPI/First_ms2.msalign: a msalign file with all MS/MS spectra in spectra_ms2.msalign using the first precursor as their precursor. 
spectra_TopMPI/Second_ms2.msalign: a msalign file with only multiplexed MS/MS spectra in spectra_ms2.msalign using the second precursor as their precursor. 
spectra_TopMPI/Primary_ms2.msalign: a msalign file with all MS/MS spectra in spectra_ms2.msalign using the TopMPI selected precursor as their precursor. 
spectra_TopMPI/Secondary_ms2.msalign: a msalign file with only multiplexed MS/MS spectra in spectra_ms2.msalign with their matched experimental peaks removed if thier corresponding PrSM is reported in Primary_toppic_prsm_single.tsv,_using the TopMPI selected secondary precursor as their precursor. 
spectra_TopMPI/First_*: The output of TopPIC searched using First_ms2.msalign while no E-value/FDR filter were applied to the output results. 
spectra_TopMPI/Second_*: The output of TopPIC searched using Second_ms2.msalign while no E-value/FDR filter were applied to the output results. 
spectra_TopMPI/Primary_toppic_prsm_single.tsv: a TSV file containing identified proteoform spectrum-matches (PrSMs) selected by TopMPI from output of First_ms2.msalign and Second_ms2.msalign with an E-value or spectrum-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
spectra_TopMPI/Primary_toppic_proteoform_single.tsv: a TSV file containing identified proteoforms selected by TopMPI from output of First_ms2.msalign and Second_ms2.msalign with an E-value or proteoform-level FDR cutoff. When an identified proteoform is shared by multiple proteins, only one protein is reported.
spectra_TopMPI/Secondary_*: The output of TopPIC searched using Secondary_ms2.msalign while a E-value/FDR filter was applied. 
To browse identified proteins, proteoforms, and PrSMs in e.g. First_ms2.msalign, use a chrome browser to open the file spectra_TopMPI/First_html/topmsv/index.html. Note that spectra_html/topfd needs to be copied to spectra_TopMPI/First_html/ in order to use this function. Google Chrome is recommended (Firefox and Edge are not recommended).
3 Command line usage
To run TopMPI, open a terminal window and run the following command. For TopPIC command flags, please see [here](https://www.toppic.org/software/toppic/manual.html). 

TopMPI.sh toppic-executable database-file-name spectrum-file-names [--TopPIC-flag=<flags>] [--TopMPI-flag=<flags>]

TopMPI Options
<!-- -h [ --help ]

Print the help message. -->

-a [ --alpha ] <a number between 0 and 1>

Set the fragmentation method(s) of MS/MS spectra. When "FILE" is selected, the fragmentation methods of spectra are given in the input spectrum data file. Default value: FILE.



Examples
Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file spectra_ms2.feature (reported by TopFD). The user does not need to specify the feature file name. TopPIC will automatically obtain the feature file name from the spectrum file name spectra_ms2.msalign.

toppic proteins.fasta spectra_ms2.msalign

Search two deconvoluted MS/MS spectrum files spectra1_ms2.msalign and spectra2_ms2.msalign against a protein database file proteins.fasta with feature files. In addition, all identifications are combined and reported using a file name "combined."

toppic -c combined proteins.fasta spectra1_ms2.msalign spectra2_ms2.msalign

Search all deconvoluted MS/MS spectrum files in the current folder against a protein database file proteins.fasta with feature files.

toppic proteins.fasta *_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta without feature files.

toppic -x proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file and a fixed modification: carbamidomethylation on cysteine.

toppic -f C57 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. In an identified proteoform, at most 2 mass shifts are allowed and the maximum allowed mass shift value is 10,000 Dalton.

toppic -s 2 -M 10000 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Two variable PTMs: oxidation on M and methylation on K are used. The modification file two_var_mods.txt can be found here.

toppic -b two_var_mods.txt proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. The error tolerance for precursor and fragment masses is 5 ppm.

toppic -e 5 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. Use the target decoy approach to compute spectrum level and proteoform-level FDRs, filter identified proteoform spectrum-matches by a 5% spectrum level FDR, and filter identified proteoforms by a 5% proteoform-level FDR.

toppic -d -t FDR -v 0.05 -T FDR -V 0.05 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign with alternating CID, HCD, and ETD spectra against a protein database file proteins.fasta with a feature file. Combine alternating CID, HCD, and ETD spectra to increase proteoform coverage.

toppic -r 3 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file. After proteoforms with unexpected mass shifts are identified, TopPIC matches the mass shifts to four common PTMs: acetylation, phosphorylation, oxidation and methylation, and uses an MIScore cutoff 0.1 to filter reported PTM sites. The modification file common_mods.txt can be found here.

toppic -B common_mods.txt -H 0.1 proteins.fasta spectra_ms2.msalign

Search a deconvoluted MS/MS spectrum file spectra_ms2.msalign against a protein database file proteins.fasta with a feature file.Use 6 CPU threads to speed up the computation.

toppic -u 6 proteins.fasta spectra_ms2.msalign