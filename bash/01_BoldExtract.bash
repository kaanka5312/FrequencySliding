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

##### SESSION -1 #####

cd /Volumes/HD-B1/Thesis/SES-1_BIDS/derivatives/afni_proc

for subj in {24..34} {36..42} ;do

echo ${subj}
cd sub-${subj}.results

	cd Zscore_data_bandpass
		# Extracting the BOLD signal and save
		
	for ROI in {1..180} {1001..1180};do
		
	3dROIstats -quiet \
	-mask /Users/kaankeskin/projects/MultGroup_WC/REPLICATION/MASK/glasser_${ROI}.nii \
	rest_Zscore+tlrc \
	> BOLD_${ROI}.1D
	done

cd /Volumes/HD-B1/Thesis/SES-1_BIDS/derivatives/afni_proc
done

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#####  S E S S I O N - 2 #################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

cd /mnt/e/Thesis/SES-2_BIDS/derivatives/afni_proc

for subj in 1 13 19 20 23 36 4 5 6 8 ;do

subj=${subj}_ses2
echo ${subj}
cd sub-${subj}.results

	cd Zscore_data_bandpass
		# Extracting the BOLD signal and save
		
	for ROI in {1..180} {1001..1180};do
	echo ${ROI}
	3dROIstats -quiet \
	-mask /home/kaanka5312/MultGroup_WC/REPLICATION/MASK/glasser_${ROI}.nii \
	rest_Zscore+tlrc \
	> BOLD_${ROI}.1D
	done

cd /mnt/e/Thesis/SES-2_BIDS/derivatives/afni_proc
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
