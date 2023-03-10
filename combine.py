import read_msalign
import sys
import pandas as pd

def main():
    args = sys.argv[1:]
    if (len(args) != 2):
        raise Exception(
            "Please pass in the ms.align file and the feature file")
    
    
    
    spec_list = read_msalign.read_spec_file(args[0])

    mzwindow = 4

    startRatio = 720

    feature_file = pd.read_csv(args[1], sep="\t")


    nomatchList = []
    for spectrum in spec_list:
        isolationWindow = (spectrum.header.spec_scan - 1) % 21
        lowerbound = startRatio + (isolationWindow - 1) * mzwindow 
        upperbound = startRatio + (isolationWindow) * mzwindow
        rt = spectrum.header.retention_time / 60

        query = feature_file[(feature_file['RepMz'] > lowerbound) & (feature_file['RepMz'] < upperbound) \
                             & (feature_file['MinElutionTime'] < rt) & (feature_file['MaxElutionTime'] > rt)]
        
        if (query.empty):
            #if (len(spectrum.peak_list) > 5):
            nomatchList.append(spectrum.header.spec_scan)
            continue
        
        maxMatch = query.loc[query['Abundance'].idxmax()]
        spectrum.header.ms_one_id = int(maxMatch['FeatureID'])
        spectrum.header.mono_mz = maxMatch['RepMz']
        spectrum.header.charge = int(maxMatch['MinCharge'])
        spectrum.header.mono_mass = maxMatch['MonoMass']
        spectrum.header.inte = maxMatch['Abundance']

    with open("modified_" + args[0], "w") as outputfile:
        for spectrum in spec_list:
            if spectrum.header.spec_scan in nomatchList:
                continue
            outputfile.write("BEGIN IONS\n")
            outputfile.write("ID=" + str(spectrum.header.spec_id) + "\n")
            outputfile.write("FRACTION_ID=0\nFILE_NAME=/home/daniel/Desktop/datafiles/RealData/yeast/study4/20220718_Yeast_DIA_720-800_75cmCapil_300nL_30kV_01.mzML\n")
            outputfile.write("SCANS=" + str(spectrum.header.spec_scan) + "\n")
            outputfile.write("RETENTION_TIME=" + str(spectrum.header.retention_time) + "\n")
            outputfile.write("LEVEL=2\n")
            outputfile.write("ACTIVATION=" + str(spectrum.header.activation) + "\n")
            outputfile.write("MS_ONE_ID=" + str(spectrum.header.ms_one_id) + "\n")
            outputfile.write("MS_ONE_SCAN=" + str(spectrum.header.ms_one_scan) + "\n")
            outputfile.write("PRECURSOR_MZ=" + str(spectrum.header.mono_mz) + "\n")
            outputfile.write("PRECURSOR_CHARGE=" + str(spectrum.header.charge) + "\n")
            outputfile.write("PRECURSOR_MASS=" + str(spectrum.header.mono_mass) + "\n")
            outputfile.write("PRECURSOR_INTENSITY=" + str(spectrum.header.inte) + "\n")
            for peak in spectrum.peak_list:
                outputfile.write(str(peak.mass) + "\t" + str(peak.intensity) + "\t" + str(peak.charge) + "\n")
            outputfile.write("END IONS\n\n")

    
        
        

        
        


        




if __name__ == "__main__":
    main()