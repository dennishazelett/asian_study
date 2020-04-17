cat GSM1663008_iEEC16_FAIRE_rep1.hg19.narrowPeak GSM1663009_iEEC16_FAIRE_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663008_009_iEEC16_FAIRE_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663008_009_iEEC16_FAIRE_combined_sorted_replicates.hg19.narrowPeak > GSM1663008_009_iEEC16_FAIRE_merged_replicates.hg19.bed
rm GSM1663008_009_iEEC16_FAIRE_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663010_iEEC16_H3K27ac_rep1.hg19.narrowPeak GSM1663011_iEEC16_H3K27ac_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663010_011_iEEC16_H3K27ac_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663010_011_iEEC16_H3K27ac_combined_sorted_replicates.hg19.narrowPeak > GSM1663010_011_iEEC16_H3K27ac_merged_replicates_hg19.bed
rm GSM1663010_011_iEEC16_H3K27ac_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663012_iEEC16_H3K4me1_rep1.hg19.narrowPeak GSM1663013_iEEC16_H3K4me1_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663012_013_iEEC16_H3K4me1_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663012_013_iEEC16_H3K4me1_combined_sorted_replicates.hg19.narrowPeak > GSM1663012_013_iEEC16_H3K4me1_merged_replicates.hg19.bed
rm GSM1663012_013_iEEC16_H3K4me1_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663016_iFTSEC246_FAIRE_rep1.hg19.narrowPeak GSM1663017_iFTSEC246_FAIRE_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663016_017_iFTSEC246_FAIRE_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663016_017_iFTSEC246_FAIRE_combined_sorted_replicates.hg19.narrowPeak > GSM1663016_017_iFTSEC246_FAIRE_merged_replicates.hg19.bed
rm GSM1663016_017_iFTSEC246_FAIRE_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663018_iFTSEC246_H3K27ac_rep1.hg19.narrowPeak GSM1663019_iFTSEC246_H3K27ac_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663018_019_iFTSEC246_H3K27ac_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663018_019_iFTSEC246_H3K27ac_combined_sorted_replicates.hg19.narrowPeak > GSM1663018_019_iFTSEC246_H3K27ac_merged_replicates_hg19.bed
rm GSM1663018_019_iFTSEC246_H3K27ac_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663020_iFTSEC246_H3K4me1_rep1.hg19.narrowPeak GSM1663021_iFTSEC246_H3K4me1_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663020_021_iFTSEC246_H3K4me1_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663020_021_iFTSEC246_H3K4me1_combined_sorted_replicates.hg19.narrowPeak > GSM1663020_021_iFTSEC246_H3K4me1_merged_replicates.hg19.bed
rm GSM1663020_021_iFTSEC246_H3K4me1_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663024_iFTSEC33_FAIRE_rep1.hg19.narrowPeak GSM1663025_iFTSEC33_FAIRE_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663024_025_iFTSEC33_FAIRE_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663024_025_iFTSEC33_FAIRE_combined_sorted_replicates.hg19.narrowPeak > GSM1663024_025_iFTSEC33_FAIRE_merged_replicates.hg19.bed
rm GSM1663024_025_iFTSEC33_FAIRE_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663026_iFTSEC33_H3K27ac_rep1.hg19.narrowPeak GSM1663027_iFTSEC33_H3K27ac_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663026_027_iFTSEC33_H3K27ac_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663026_027_iFTSEC33_H3K27ac_combined_sorted_replicates.hg19.narrowPeak > GSM1663026_027_iFTSEC33_H3K27ac_merged_replicates.hg19.bed
rm GSM1663026_027_iFTSEC33_H3K27ac_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663028_iFTSEC33_H3K4me1_rep1.hg19.narrowPeak GSM1663029_iFTSEC33_H3K4me1_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663028_029_iFTSEC33_H3K4me1_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663028_029_iFTSEC33_H3K4me1_combined_sorted_replicates.hg19.narrowPeak > GSM1663028_029_iFTSEC33_H3K4me1_merged_replicates.hg19.bed
rm GSM1663028_029_iFTSEC33_H3K4me1_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663033_iOSE11_FAIRE_rep1.hg19_rCRSchrm.narrowPeak GSM1663034_iOSE11_FAIRE_rep2.hg19_rCRSchrm.narrowPeak | sort -k 1,1 -k2,2n > GSM1663033_034_iOSE11_FAIRE_combined_sorted_replicates.hg19_rCRSchrm.narrowPeak
mergeBed -i GSM1663033_034_iOSE11_FAIRE_combined_sorted_replicates.hg19_rCRSchrm.narrowPeak > GSM1663033_034_iOSE11_FAIRE_merged_replicates.hg19_rCRSchrm.bed
rm GSM1663033_034_iOSE11_FAIRE_combined_sorted_replicates.hg19_rCRSchrm.narrowPeak

cat GSM1663035_iOSE11_H3K27ac_rep1.hg19.narrowPeak GSM1663036_iOSE11_H3K27ac_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663035_036_iOSE11_H3K27ac_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663035_036_iOSE11_H3K27ac_combined_sorted_replicates.hg19.narrowPeak > GSM1663035_036_iOSE11_H3K27ac_merged_replicates.hg19.bed
rm GSM1663035_036_iOSE11_H3K27ac_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663037_iOSE11_H3K4me1_rep1.hg19.narrowPeak GSM1663038_iOSE11_H3K4me1_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663037_038_iOSE11_H3K4me1_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663037_038_iOSE11_H3K4me1_combined_sorted_replicates.hg19.narrowPeak > GSM1663037_038_iOSE11_H3K4me1_merged_replicates.hg19.bed
rm GSM1663037_038_iOSE11_H3K4me1_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663041_iOSE4_FAIRE_rep1.hg19.narrowPeak GSM1663042_iOSE4_FAIRE_rep2.hg19.narrowPeak GSM1663043_iOSE4_FAIRE_rep3.hg19.narrowPeak GSM1663044_iOSE4_FAIRE_rep4.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663041_042_043_044_iOSE4_FAIRE_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663041_042_043_044_iOSE4_FAIRE_combined_sorted_replicates.hg19.narrowPeak > GSM1663041_042_043_044_iOSE4_FAIRE_merged_replicates.hg19.bed
rm GSM1663041_042_043_044_iOSE4_FAIRE_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663045_iOSE4_H3K27ac_rep1.hg19.narrowPeak GSM1663046_iOSE4_H3K27ac_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663045_046_iOSE4_H3K27ac_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663045_046_iOSE4_H3K27ac_combined_sorted_replicates.hg19.narrowPeak > GSM1663045_046_iOSE4_H3K27ac_merged_replicates.hg19.bed
rm GSM1663045_046_iOSE4_H3K27ac_combined_sorted_replicates.hg19.narrowPeak

cat GSM1663047_iOSE4_H3K4me1_rep1.hg19.narrowPeak GSM1663048_iOSE4_H3K4me1_rep2.hg19.narrowPeak | sort -k 1,1 -k2,2n > GSM1663047_048_iOSE4_H3K4me1_combined_sorted_replicates.hg19.narrowPeak
mergeBed -i GSM1663047_048_iOSE4_H3K4me1_combined_sorted_replicates.hg19.narrowPeak > GSM1663047_048_iOSE4_H3K4me1_merged_replicates.hg19.bed
rm GSM1663047_048_iOSE4_H3K4me1_combined_sorted_replicates.hg19.narrowPeak

cat ENCFF025JZZ-HCT116-CTCF-MYERS-IDR-HG19-narrowPeak.bed ENCFF340EDA-HCT116-CTCF-STAM-IDR-HG19-narrowPeak.bed ENCFF364QXM-HCT116-CTCF-STAM-IDR-HG19-narrowPeak.bed ENCFF917ZPO-HCT116-CTCF-MYERS-IDR-HG19-narrowPeak.bed | sort -k 1,1 -k2,2n > ENCODE_HCT116_CTCF_combined_sorted_replicates.hg19.bed
mergeBed -i ENCODE_HCT116_CTCF_combined_sorted_replicates.hg19.bed > ENCODE_HCT116_CTCF_merged_replicates.hg19.bed
rm ENCODE_HCT116_CTCF_combined_sorted_replicates.hg19.bed

cat ENCFF282NVG-HCT116-H3K27ac-BERNSTEIN-replicated-HG19-narrowPeak.bed ENCFF590YIB-HCT116-H3K27ac-FARNHAM-replicated-HG19-narrowPeak.bed | sort -k 1,1 -k2,2n > ENCODE_HCT116_H3K27ac_combined_sorted_replicates.hg19.bed
mergeBed -i ENCODE_HCT116_H3K27ac_combined_sorted_replicates.hg19.bed > ENCODE_HCT116_H3K27ac_merged_replicates.hg19.bed
rm ENCODE_HCT116_H3K27ac_combined_sorted_replicates.hg19.bed

mv GSM1663008_iEEC16_FAIRE_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663009_iEEC16_FAIRE_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663010_iEEC16_H3K27ac_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663011_iEEC16_H3K27ac_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663012_iEEC16_H3K4me1_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663013_iEEC16_H3K4me1_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663016_iFTSEC246_FAIRE_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663017_iFTSEC246_FAIRE_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663018_iFTSEC246_H3K27ac_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663019_iFTSEC246_H3K27ac_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663020_iFTSEC246_H3K4me1_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663021_iFTSEC246_H3K4me1_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663024_iFTSEC33_FAIRE_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663025_iFTSEC33_FAIRE_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663026_iFTSEC33_H3K27ac_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663027_iFTSEC33_H3K27ac_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663028_iFTSEC33_H3K4me1_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663029_iFTSEC33_H3K4me1_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663033_iOSE11_FAIRE_rep1.hg19_rCRSchrm.narrowPeak ../merged_replicates
mv GSM1663034_iOSE11_FAIRE_rep2.hg19_rCRSchrm.narrowPeak ../merged_replicates
mv GSM1663035_iOSE11_H3K27ac_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663036_iOSE11_H3K27ac_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663037_iOSE11_H3K4me1_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663038_iOSE11_H3K4me1_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663041_iOSE4_FAIRE_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663042_iOSE4_FAIRE_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663043_iOSE4_FAIRE_rep3.hg19.narrowPeak ../merged_replicates
mv GSM1663044_iOSE4_FAIRE_rep4.hg19.narrowPeak ../merged_replicates
mv GSM1663045_iOSE4_H3K27ac_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663046_iOSE4_H3K27ac_rep2.hg19.narrowPeak ../merged_replicates
mv GSM1663047_iOSE4_H3K4me1_rep1.hg19.narrowPeak ../merged_replicates
mv GSM1663048_iOSE4_H3K4me1_rep2.hg19.narrowPeak ../merged_replicates
mv ENCFF025JZZ-HCT116-CTCF-MYERS-IDR-HG19-narrowPeak.bed ../merged_replicates
mv ENCFF340EDA-HCT116-CTCF-STAM-IDR-HG19-narrowPeak.bed ../merged_replicates
mv ENCFF364QXM-HCT116-CTCF-STAM-IDR-HG19-narrowPeak.bed ../merged_replicates
mv ENCFF917ZPO-HCT116-CTCF-MYERS-IDR-HG19-narrowPeak.bed ../merged_replicates
mv ENCFF282NVG-HCT116-H3K27ac-BERNSTEIN-replicated-HG19-narrowPeak.bed ../merged_replicates
mv ENCFF590YIB-HCT116-H3K27ac-FARNHAM-replicated-HG19-narrowPeak.bed ../merged_replicates
