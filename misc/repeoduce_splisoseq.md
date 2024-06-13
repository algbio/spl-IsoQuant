## Data simulation

Transcript quantification (script taken from [https://github.com/andrewprzh/lrgasp-simulation](https://github.com/andrewprzh/lrgasp-simulation)):

```
seqkit sample -p 0.012 4337_combined.fastq.gz > 4337.1M.fq
```

```
python3 spl-IsoQuant/misc/quantify.py -t 16 \
--reference_transcripts gencode.v35.basic.transcripts.fa \
--fastq 4337.1M.fq --output 4337.counts.tsv
```

cDNA template generation:

```
python3 spl-IsoQuant/misc/generate_umi_list.py > simulated_umis.tsv
```

```
python3 spl-IsoQuant/misc/simulate_barcoded.py \
--transcriptome gencode.v35.basic.transcripts.fa \
--counts 4337.counts.tsv --template_count 5000000 --mode spatial \
--barcodes A0079_044_BeadBarcodes.txt --umis simulated_umis.tsv -o Human.A0079_044.templates
```


Modified NanoSim version is available here [https://github.com/andrewprzh/lrgasp-simulation](https://github.com/andrewprzh/lrgasp-simulation).

NanoSim model generation:

```
python3 lrgasp-simulation/src/NanoSim/src/read_analysis.py transcriptome \
-rg GRCh38.chr.fa  -rt gencode.v36.basic.transcripts.fa \
-annot gencode.v36.basic.annotation.gtf -a minimap2 -i  4337.1M.fq \
--no_intron_retention -t 16 -o ONT.4337/model 
```


Simulation: 

```
python3 lrgasp-simulation/src/NanoSim/src/simulator.py transcriptome \
--ref_t Human.A0079_044.templates.fasta --exp 4337.counts.tsv \
--model_prefix ONT.4337/model -n 20000000 -b guppy -r cDNA_1D \
--no_model_ir --aligned_only --truncation_mode ont_r9 --output ONT.4337_simulated
```


## Running Spl-IsoQuant


Simulated data:

```
python3 spl-IsoQuant/splisoquant.py --reference GRCh38.chr.fa \
--complete_genedb --genedb gencode.v35.annotation.gtf \
--barcode_whitelist A0079_044_BeadBarcodes.txt \
--fastq ONT.4337_simulated_aligned_reads.fasta \
-t 16 -d ont -p ONT.Simulated.4337 -o Simulated.A0079_044.4337
 
```


Real data:

```
python3 spl-IsoQuant/splisoquant.py --reference GRCh38.chr.fa \
--complete_genedb --genedb gencode.v35.annotation.gtf \
--barcode_whitelist <BeadBarcodes.txt> \
--fastq <samples.fastq.gz> \
-t 16 -d ont -p <sample_name> -o <SplIsoQuant.sample_name>
 
```


