# alevin em
/mnt/scratch5/laraib/alevin2/salmon/bin/salmon alevin -lISR -1 ./4k_simulation/hg_100_S1_L001_R1_001.fastq.gz -2 ./4k_simulation/hg_100_S1_L001_R2_001.fastq.gz -o ./em -i /mnt/scratch7/avi/snare-seq/indices/human_index --tgMap /mnt/scratch7/avi/snare-seq/reference/human/txp2gene.txt --chromium --whitelist ./4k_simulation/minnow_whitelist.txt

# seurat mapping
for i in {1..30}
do
	Rscript --vanilla run_seurat.R
	mkdir priorNum"$i"
	mv ./*.txt priorNum"$i"
	cd priorNum"$i"

	paste <(cut -f1 -d'_' refs.txt) <(cut -f1 -d'_' query.txt) > mappings.txt
	awk '($1 != $2)' mappings.txt > mappings_uniq.txt
	awk '($1 != $2)' mappings.txt | cut -f1 | sort | uniq > whitelist.txt
	cd ../

done

# instead of average, use all mappings together
cat priorNum{1..30}/mappings_uniq.txt | sort | uniq > mappings.txt
mkdir prior
python generate_prior.py

## alevin vbem
/mnt/scratch5/laraib/alevin2/salmon/bin/salmon alevin -lISR -1 ./4k_simulation/hg_100_S1_L001_R1_001.fastq.gz -2 ./4k_simulation/hg_100_S1_L001_R2_001.fastq.gz -o ./vbem -i /mnt/scratch7/avi/snare-seq/indices/human_index --tgMap /mnt/scratch7/avi/snare-seq/reference/human/txp2gene.txt --chromium --whitelist ./4k_simulation/minnow_whitelist.txt --vbemPrior prior --vbemNorm 20055723

