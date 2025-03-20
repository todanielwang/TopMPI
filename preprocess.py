import pandas as pd
import util
import os
import read_msalign

def preprocess(basedir, commonprefix, outputdir, featurecutoff=0.5, isFeature=True):
    spectra = read_msalign.read_spec_file(os.path.join(basedir, f"{commonprefix}_ms2.msalign"))

    if isFeature:
        features = util.read_tsv(os.path.join(basedir, f"{commonprefix}_ms1.feature"))
        goodfeatures = features[features["EC_score"] > featurecutoff]["Feature_ID"].tolist()

        filtered_spectra = []

        for spec in spectra[:]:  
            valid_indices = [i for i, pre_id in enumerate(spec.header.pre_id_list) if int(pre_id) in goodfeatures]

            if valid_indices:
                spec.header.pre_mz_list = [spec.header.pre_mz_list[i] for i in valid_indices]
                spec.header.pre_charge_list = [spec.header.pre_charge_list[i] for i in valid_indices]
                spec.header.pre_mass_list = [spec.header.pre_mass_list[i] for i in valid_indices]
                spec.header.pre_inte_list = [spec.header.pre_inte_list[i] for i in valid_indices]
                spec.header.pre_id_list = [spec.header.pre_id_list[i] for i in valid_indices]

                filtered_spectra.append(spec)  # Keep valid objects
            else:
                del spec  # Explicitly delete invalid objects

        # Remove all references from spectra
        spectra.clear()  

        # Assign only the valid objects back
        spectra.extend(filtered_spectra)


    for spec in spectra:
        while len(spec.header.pre_mz_list) > 1:
            if not int(spec.header.pre_charge_list[0]) == int(spec.header.pre_charge_list[1]):
                del spec.header.pre_mz_list[1]
                del spec.header.pre_charge_list[1]
                del spec.header.pre_mass_list[1]
                del spec.header.pre_inte_list[1]
                del spec.header.pre_id_list[1]

    read_msalign.write_spec_file(os.path.join(outputdir, "First_ms2.msalign"))

    os.rename(os.path.join(outputdir, "First_ms2_modified.msalign"), os.path.join(outputdir, "First_ms2.msalign"))


