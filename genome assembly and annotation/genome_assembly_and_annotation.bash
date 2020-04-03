# genome assembly 
canu -p barthii -d barthii_corERate0.060_GS411m_blasr5.3.0MinReadLen500  genomeSize=349m minThreads=36 \
minMemory=110 minReadLength=500 corMinCoverage=0 MhapSensitivity=high correctedErrorRate=0.060 -pacbio-raw barthii_all_subreads.fq

# gene predictions
augustus --strand=both --genemodel=complete --singlestrand=false --UTR=on --noInFrameStop=true --protein=on --introns=on --start=on \
--stop=on --cds=on --codingseq=on --alternatives-from-sampling=true --gff3=on --outfile=barthii_sm.gff --species=rice barthii.sm.fa

glimmerhmm ./barthii_Chr1.fa -d /data4/pub/tool/GlimmerHMM/trained_dir/rice/ -g -f -n 1 -o barthii_Chr1.sm_glimmerhmm


snap /pub/tool/snap/HMM/O.sativa.hmm ../rice barthii.sm.fa -gff -aa rice barthii.sm.pep -tx \
rice barthii.sm.gene.transcript > rice barthii.sm.gene.gff


exonerate --model protein2genome --ryo '#PIS\t%pi\t%ps\t%qal\t%ql\n' --score 800 --saturatethreshold 100 --dnahspthreshold 60 \
--dnawordlen 14 --forwardcoordinates FALSE --softmasktarget TRUE --showtargetgff TRUE --showquerygff FALSE --showalignment FALSE \
--showvulgar FALSE --showcigar FALSE  R498_CoreSet.pros.Chr1.fasta barthii_Chr1_sm.fa > barthii_Chr1_R498_CoreSet.exonerate

#TE annotation
RepeatMasker -pa 24 -div 22 -cutoff 250 -gff -lib barthiiSpecificTelib.fa -nolow -norna -xsmall barthii.fa
