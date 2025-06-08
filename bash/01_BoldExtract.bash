#!/bin/bash
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#####  B I P O L A R #####################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# Linux
cd /media/kaanka5312/HD-B1/BPB_proc

for subj in {1..2} {4..9} {10..13} {15..23};do

	echo ${subj}
	cd sub-${subj}.results

		cd Zscore_data_bandpass
		# Extracting the BOLD signal and save
		
		for ROI in {1..180} {1001..1180};do
		
		3dROIstats -quiet \
		-mask /home/kaanka5312/projects/MultGroup_WC/REPLICATION/MASK/glasser_${ROI}.nii \
		rest_Zscore+tlrc \
		> BOLD_${ROI}.1D
		done

cd /media/kaanka5312/HD-B1/BPB_proc	
done

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#####  C O N T R O L #####################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# Linux
cd /media/kaanka5312/HD-B1/BIDS/derivatives/afni/

for subj in {35..53} {55..57} {59..70};do

	echo ${subj}
	cd sub-${subj}.results

		cd Zscore_data_bandpass
		# Extracting the BOLD signal and save
		
		for ROI in {1..180} {1001..1180};do
		
		3dROIstats -quiet \
		-mask /home/kaanka5312/projects/MultGroup_WC/REPLICATION/MASK/glasser_${ROI}.nii \
		rest_Zscore+tlrc \
		> BOLD_${ROI}.1D
		done

cd /media/kaanka5312/HD-B1/BIDS/derivatives/afni/	
done
