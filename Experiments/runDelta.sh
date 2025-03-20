for fraction in {1..6}; do
    for gamma in {4..9}; do
        python3 noisefilter.py /home/daniel/Desktop/datafiles/RealData/SW480Noise/F${fraction}/480_F${fraction}_01_TopMPI/ 100 ${gamma}
	mv /home/daniel/Desktop/datafiles/RealData/SW480Noise/F${fraction}/480_F${fraction}_01_TopMPI/FirstPrSM_full.tsv /home/daniel/Desktop/datafiles/RealData/SW480Noise/F${fraction}/480_F${fraction}_01_TopMPI/FirstPrSM_full${gamma}.tsv
    done
done
