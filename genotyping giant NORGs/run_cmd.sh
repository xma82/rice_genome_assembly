
#Step1. Counting read number over breakpoints or integration sites.
#First parameter represents inputfile.
#Second parameter represents type of NORGs relative to reference genome:ins or del.
#Third parameter represents BAM file directory. 
python s1_ReadNumber.py norg_del.txt del ref bam_path
python s1_ReadNumber.py norg_del.txt del nonref bam_path

python s1_ReadNumber.py norg_ins.txt ins ref bam_path
python s1_ReadNumber.py norg_ins.txt ins nonref bam_path

#input files (norg_del.txt or norg_ins.txt) format as follows:
# Chr1    51174   53134   1961    -2      1963    Chr1    58647   58644   NORG_n1    N_0     del     
# Chr1    54581   54581   0    5001      5001    Chr1    60078   65078   NORG_n2     N_0    ins



#Step2. According to read number, genotyping  giant NORGs
python s2_genotype.py del del_ref_total_read.txt  del_nonref_total_read.txt del_geno.txt

python s2_genotype.py ins ins_ref_total_read.txt ins_nonref_total_read.txt ins_geno.txt