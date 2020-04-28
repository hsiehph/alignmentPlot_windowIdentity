import os
import itertools

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
configfile: "config.yaml"
shell.prefix("source config.sh; set -eo pipefail ; ")

FASTA = config["FASTA"]
add_line_onY = float(config["addline_onY"])
window_step = [(2000, 100), (500,100), (1000, 100)]
runByPairFASTA = config["runByPairFASTA"]

if not os.path.exists("log"):
    os.makedirs("log")
CUTOFF = 0.998



dict_ID_fasta = {}
if not runByPairFASTA:
	for ID in config["FASTA"]:
		dict_ID_fasta[ID] = config["FASTA"][ID]
else:
	if not os.path.exists("pairedFASTA"):
		os.mkdir("pairedFASTA")

	for a,b in itertools.combinations(config["FASTA"].keys(), 2):
		ID = a + "_" + b
		outFile = os.path.join("pairedFASTA/", ID + ".fa")
		if os.path.exists(outFile):
			dict_ID_fasta[ID] = outFile
			continue
		with open(outFile, "w") as fout, open(config["FASTA"][a]) as f1, open(config["FASTA"][b]) as f2:
			for line in f1:
				fout.write(line)
			for line in f2:
				fout.write(line)
			dict_ID_fasta[ID] = outFile
		

rule all:
	input: expand("{ID}/fracIdentity_w{win_step[0]}s{win_step[1]}.HighIDregions.bed", win_step = window_step, ID=dict_ID_fasta.keys()), 
			expand("plots/{ID}_frac_identity_plot.pdf", ID=dict_ID_fasta.keys()),


rule plot_fraction_Identity:
	input: "{ID}/fracIdentity.cat.txt"
	output: "plots/{ID}_frac_identity_plot.pdf"
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00", addline_y=add_line_onY
	shell:
		"Rscript {SNAKEMAKE_DIR}/plot_identity.r {input} {params.addline_y} {output}"


rule catIdentityFiles:
	input: ["{ID}/fracIdentity_w%ss%s.txt" % (x,y) for x,y in window_step]
	output: "{ID}/fracIdentity.cat.txt"
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00"
	shell: """ cat {input} >> {output} """

rule bedfile_HighIDregions:
	input: ["{ID}/fracIdentity_w%ss%s.txt" % (x,y) for x,y in window_step]
	output: ["{ID}/fracIdentity_w%ss%s.HighIDregions.bed" % (x,y) for x,y in window_step]
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00", cutoff=CUTOFF
	run:
		for pair in window_step:
			shell(""" awk '{{OFS="\\t"}} {{if($3 > {params.cutoff}) print "chrAA", $1, $2, $3, $4}}' {wildcards.ID}/fracIdentity_w%ss%s.txt > {wildcards.ID}/tmp_fracIdentity_w%ss%s.HighIDregions.bed """ % (pair[0], pair[1], pair[0], pair[1]))
			with open("%s/tmp_fracIdentity_w%ss%s.HighIDregions.bed" % (wildcards.ID, pair[0], pair[1])) as f:
				while True:
					line = f.readline()
					if line !=  "":
						shell(""" bedtools merge -i {wildcards.ID}/tmp_fracIdentity_w%ss%s.HighIDregions.bed -c 4 -o max > {wildcards.ID}/fracIdentity_w%ss%s.HighIDregions.bed """ % (pair[0], pair[1], pair[0], pair[1]))
						shell(""" rm {wildcards.ID}/tmp_fracIdentity_w%ss%s.HighIDregions.bed """ % (pair[0], pair[1]))
						break
					else:
						shell(""" mv {wildcards.ID}/tmp_fracIdentity_w%ss%s.HighIDregions.bed {wildcards.ID}/fracIdentity_w%ss%s.HighIDregions.bed """ % (pair[0], pair[1], pair[0], pair[1]))
						break

rule calcIdentity:
	input: "{ID}/{ID}.aln.fa"
	output: ["{ID}/fracIdentity_w%ss%s.txt" % (x,y) for x,y in window_step]
	params: sge_opts = "-l mfree=4G -l h_rt=24:00:00"
	run: 
		for pair in window_step:
			shell("python {SNAKEMAKE_DIR}/calc_windowIdentity.py {input} --window {pair[0]} --step {pair[1]} > {wildcards.ID}/fracIdentity_w%ss%s.txt" % (pair[0], pair[1]))

rule rename:
	input: "{ID}/{ID}.best.fas"
	output: "{ID}/{ID}.aln.fa"
	params: sge_opts = "-l mfree=1G -l h_rt=24:00:00"
	run: 
		shell("mv {input} {output} ; touch {input} ; touch {output}")

rule alignment:
	input: lambda wildcards: dict_ID_fasta[wildcards.ID]
	output: "{ID}/{ID}.best.fas"
	params: sge_opts = "-l mfree=8G -pe serial 6 -l h_rt=24:00:00", prefix="{ID}"
	run: 
#		shell("prank -d={input} -F -o={ID}/{params.prefix}")
		shell("mafft --thread $NSLOTS --memsave --maxiterate 100 {input} > {output}")
