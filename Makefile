# TODO: Add a bismark-ify step (and stopping if column not present)

SCRIPTS := /project/silvers_dev/iSA/scripts
OUTDIR := /project/silvers_dev/iSA/results
LOGS := $(OUTDIR)/logs
INPUT ?= /project/silvers_dev/iSA/seqdir
NCORES ?= 8

BISMARK_METHYLATION_EXTRACTOR := /package/sequencer/bin/bismark_methylation_extractor

# Assuming process all BAM files in INPUT directory
BAM_FILES := $(wildcard $(INPUT)/*.bam)
WAT_BAM := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.bam, $(BAM_FILES))
CRI_BAM := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.bam, $(BAM_FILES))
WAT_BED := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.bed, $(BAM_FILES))
CRI_BED := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.bed, $(BAM_FILES))
ISA_BED := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.iSA.bed, $(BAM_FILES))
ISA_REPORT := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.iSA.report.txt, $(BAM_FILES))
HEADER := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.header, $(BAM_FILES))
WAT_LN := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.LN, $(BAM_FILES))
CRI_LN := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.LN, $(BAM_FILES))
WAT_ISA_BAM := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.iSA.bam, $(BAM_FILES))
CRI_ISA_BAM := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.iSA.bam, $(BAM_FILES))
WAT_EXTRACT := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.extract/.done, $(BAM_FILES))
CRI_EXTRACT := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.extract/.done, $(BAM_FILES))
WAT_ISA_IDS := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Wat.ID, $(BAM_FILES))
CRI_ISA_IDS := $(patsubst $(INPUT)/%.bam, $(OUTDIR)/%.Cri.ID, $(BAM_FILES))

.PHONY: all clean

all: directories divide_alignments write_alignment_beds annealing_strands extract_paired_alignments extract_methylation_calls

directories: $(OUTDIR) $(LOGS)

$(OUTDIR) $(LOGS):
	mkdir -p $@

# Divide alignments for strands

.PHONY: divide_alignments

divide_alignments: $(WAT_BAM) $(CRI_BAM)

# -- bismark --
# $(OUTDIR)/%.Wat.bam: $(INPUT)/%.bam
# 	samtools view -h $< | awk '$$3!="chrM" && $$3!="ChrM" && $$16!="XG:Z:GA"' | samtools view -bh -@ $(NCORES) - > $@

# $(OUTDIR)/%.Cri.bam: $(INPUT)/%.bam
# 	samtools view -h $< | awk '$$3!="chrM" && $$3!="ChrM" && $$16!="XG:Z:CT"' | samtools view -bh -@ $(NCORES) - > $@

# -- bsmap --
$(OUTDIR)/%.Wat.bam: $(INPUT)/%.bam
	samtools view -h $< | awk '$$1 ~ /^@/ || $$3!="chrM" && ($$15=="ZS:Z:++" || $$15=="ZS:Z:+-")' | samtools view -bh -@ $(NCORES) - > $@

$(OUTDIR)/%.Cri.bam: $(INPUT)/%.bam
	samtools view -h $< | awk '$$1 ~ /^@/ || $$3!="chrM" && ($$15=="ZS:Z:-+" || $$15=="ZS:Z:--")' | samtools view -bh -@ $(NCORES) - > $@

# Writing bed files for strands

.PHONY: write_alignment_beds

write_alignment_beds: $(WAT_BED) $(CRI_BED)

$(OUTDIR)/%.Wat.bed: $(OUTDIR)/%.Wat.bam
	bamToBed -bedpe -i $< | awk -v t1=0 -v t2=0 '{if($$2<t1){$$9=0} else{$$9=$$2-t1}; print $$1,$$9,$$6+t2,NR,$$7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -u -S 64G > $@;

$(OUTDIR)/%.Cri.bed: $(OUTDIR)/%.Cri.bam
	bamToBed -bedpe -i $< | awk -v t1=0 -v t2=0 '{if($$2<t2){$$9=0} else{$$9=$$2-t2}; print $$1,$$9,$$6+t1,NR,$$7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -u -S 64G > $@;

# Annealing Watson and Crick strands

.PHONY: annealing_strands

annealing_strands: $(ISA_BED) $(ISA_REPORT)

$(OUTDIR)/%.iSA.bed: $(OUTDIR)/%.Wat.bed $(OUTDIR)/%.Cri.bed
	intersectBed -a $(OUTDIR)/$*.Wat.bed -b $(OUTDIR)/$*.Cri.bed -wa -wb -f 1 -F 1 -sorted | cut -f 1-5,9,10 > $@;

$(OUTDIR)/%.iSA.report.txt: $(OUTDIR)/%.iSA.bed
	$(SCRIPTS)/prepare_report.sh $(OUTDIR)/$*;

# Extracting paired alignments

.PHONY: extract_paired_alignments

extract_paired_alignments: $(HEADER) $(WAT_LN) $(CRI_LN) $(WAT_ISA_BAM) $(CRI_ISA_BAM)

$(OUTDIR)/%.header: $(OUTDIR)/%.Wat.bam
	samtools view -H $< > $@;

$(OUTDIR)/%.Wat.LN: $(OUTDIR)/%.iSA.bed
	cut -f 4 $< | awk '{print $$1*2-1"\n"$$1*2}' > $@;

$(OUTDIR)/%.Cri.LN: $(OUTDIR)/%.iSA.bed
	cut -f 6 $< | awk '{print $$1*2-1"\n"$$1*2}' > $@;

# -- bismark --
# $(OUTDIR)/%.Wat.iSA.bam: $(OUTDIR)/%.Wat.LN $(OUTDIR)/%.header
# 	samtools view $(OUTDIR)/$*.Wat.bam | awk 'FNR==NR {h[$$1];next} (FNR in h)' $< - | cat $(OUTDIR)/$*.header - | samtools view -bh -@ $(NCORES) - > $@;

# $(OUTDIR)/%.Cri.iSA.bam: $(OUTDIR)/%.Cri.LN $(OUTDIR)/%.header
# 	samtools view $(OUTDIR)/$*.Cri.bam | awk 'FNR==NR {h[$$1];next} (FNR in h)' $< - | cat $(OUTDIR)/$*.header - | samtools view -bh -@ $(NCORES) - > $@;

# -- bsmap --
$(OUTDIR)/%.Wat.iSA.bam: $(OUTDIR)/%.Wat.LN $(OUTDIR)/%.header
	samtools view $(OUTDIR)/$*.Wat.bam | awk 'FNR==NR {h[$$1];next} (FNR in h)' $< - | cat $(OUTDIR)/$*.header - | samtools view -bh -@ $(NCORES) - > $@.temp && \
	samtools view -h $@.temp | awk 'BEGIN {OFS="\t"} {if ($$1 ~ /^@/) {print $$0} else {split($$1, a, ":"); $$1=a[1]":"a[2]":"a[3]":"a[4]; print $$0}}' | samtools view -b -@ $(NCORES) -o $@ && \
	rm $@.temp;

$(OUTDIR)/%.Cri.iSA.bam: $(OUTDIR)/%.Cri.LN $(OUTDIR)/%.header
	samtools view $(OUTDIR)/$*.Cri.bam | awk 'FNR==NR {h[$$1];next} (FNR in h)' $< - | cat $(OUTDIR)/$*.header - | samtools view -bh -@ $(NCORES) - > $@.temp && \
	samtools view -h $@.temp | awk 'BEGIN {OFS="\t"} {if ($$1 ~ /^@/) {print $$0} else {split($$1, a, ":"); $$1=a[1]":"a[2]":"a[3]":"a[4]; print $$0}}' | samtools view -b -@ $(NCORES) -o $@ && \
	rm $@.temp;

# Extracting DNA methylation calls

.PHONY: extract_methylation_calls

extract_methylation_calls: $(WAT_EXTRACT) $(CRI_EXTRACT) $(WAT_ISA_IDS) $(CRI_ISA_IDS)

$(OUTDIR)/%.Wat.extract/.done: $(OUTDIR)/%.Wat.iSA.bam
	@mkdir -p $(dir $@) && $(BISMARK_METHYLATION_EXTRACTOR) -p --no_header --gzip $< --output $(dir $@) && touch $@;

$(OUTDIR)/%.Cri.extract/.done: $(OUTDIR)/%.Cri.iSA.bam
	@mkdir -p $(dir $@) && $(BISMARK_METHYLATION_EXTRACTOR) -p --no_header --gzip $< --output $(dir $@) && touch $@;

# Pairing DNA methylation calls between Watson and Crick
# -- bismark --
# $(OUTDIR)/%.Wat.ID: $(OUTDIR)/%.iSA.bed
# 	awk '{print $$5"\t"NR}' $< | sort -k1,1 > $@;

# $(OUTDIR)/%.Cri.ID: $(OUTDIR)/%.iSA.bed
# 	awk '{print $$7"\t"NR}' $< | sort -k1,1 > $@;

# -- bsmap --
$(OUTDIR)/%.Wat.ID: $(OUTDIR)/%.iSA.bed
	awk 'BEGIN {OFS="\t"} {if ($$5 ~ /^@/) {print $$0} else {split($$5, a, ":"); $$5=a[1]":"a[2]":"a[3]":"a[4]; print $$0}}' $< | awk '{print $$5"\t"NR}' | sort -k1,1 > $@;

$(OUTDIR)/%.Cri.ID: $(OUTDIR)/%.iSA.bed
	awk 'BEGIN {OFS="\t"} {if ($$5 ~ /^@/) {print $$0} else {split($$5, a, ":"); $$5=a[1]":"a[2]":"a[3]":"a[4]; print $$0}}' $< | awk '{print $$5"\t"NR}' | sort -k1,1 > $@;

$(OUTDIR)/%.Wat.CG.me: $(OUTDIR)/%.Wat.ID $(OUTDIR)/%.Wat.extract/.done
	gzip -cd CpG_OT_$1.Wat.iSA.txt.gz | sort -k1,1 | join -j 1 - $< | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4,$3,$4-1,$4+1,$5}' OFS='\t' | sort -k1,1 > $@;

$(OUTDIR)/%.Wat.CHG.me: $(OUTDIR)/%.Wat.ID $(OUTDIR)/%.Cri.extract/.done
	gzip -cd CHG_OT_$1.Wat.iSA.txt.gz | sort -k1,1 | join -j 1 - $< | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4,$3,$4-1,$4+2,$5}' OFS='\t' | sort -k1,1 > $@;

clean:
	rm -rf $(WAT_BAM) $(CRI_BAM) $(WAT_BED) $(CRI_BED) $(ISA_BED) $(ISA_REPORT) \
	$(HEADER) $(WAT_LN) $(CRI_LN) $(WAT_ISA_BAM) $(CRI_ISA_BAM) $(WAT_EXTRACT) \
	$(CRI_EXTRACT) $(WAT_ISA_IDS) $(CRI_ISA_IDS)