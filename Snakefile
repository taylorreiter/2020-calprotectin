##########################
## IMPORTANT NOTE
##########################

# This snakefile is built on intermediate
# outputs from github.com/taylorreiter/ibd
# The rules that generate these intermediate 
# outputs are not included in this snakefile,
# and include preprocessing and sourmash
# signature generation. 
# This Snakefile is a proof-of-concept for
# correlating k-mers with quantitative/metabolite
# measurements. 

import pandas as pd
import feather
from sourmash import signature
import glob
import os
from collections import Counter

m = pd.read_csv("inputs/working_metadata_fecal_calprotectin.tsv", sep = "\t", header = 0)
SAMPLES = m['sample'].unique().tolist()

rule all:
    input:
        "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        "outputs/vita_rf/pomona_install.txt",
        #"outputs/gather/vita_vars.csv",

########################################
## PREPROCESSING
########################################

# already done by Snakefile in ibd workflow
# See github.com/taylorreiter/ibd

########################################
## Filtering and formatting signatures
########################################


# made new repo. copy over signatures from HMP and libs to to sigs 
# folder. This way, the outputs structure can remain the same 
# and not conflict/overwrite existing files, and combine the libs and
# hmp vars into one workflow.

rule get_greater_than_1_filt_sigs:
    input: expand("outputs/sigs/{sample}.sig", sample = SAMPLES) 
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    run:
        # Determine the number of hashes, the number of unique hashes, and the number of
        # hashes that occur once across 954 IBD/control gut metagenomes (excludes the 
        # iHMP). Calculated for a scaled of 2k. 9 million hashes is the current 
        # approximate upper limit with which to build a sample vs hash abundance table 
        # using my current methods.

        files = input

        all_mins = []
        for file in files:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[1]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes 
                all_mins += mins

        counts = Counter(all_mins) # tally the number of hashes

        # remove hashes that occur only once
        for hashes, cnts in counts.copy().items():
            if cnts < 2:
                counts.pop(hashes)

        # write out hashes to a text file
        with open(str(output), "w") as f:
            for key in counts:
                print(key, file=f)


rule convert_greater_than_1_hashes_to_sig:
    input: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{sample}.sig"
    output: "outputs/filt_sigs/{sample}_filt.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{sample}_filt.sig"
    output: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.sample}_filt
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{sample}_filt_named.csv"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{sample}_filt_named.csv", sample = SAMPLES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/normalized_hash_abund_long.R"

rule make_hash_abund_table_wide:
    input: "outputs/hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    run:
        import pandas as pd
        import feather
        
        ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
        ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
        ibd_wide = ibd_wide.fillna(0)
        ibd_wide['sample'] = ibd_wide.index
        ibd_wide = ibd_wide.reset_index(drop=True)
        ibd_wide.columns = ibd_wide.columns.astype(str)
        ibd_wide.to_feather(str(output)) 


########################################
## Random forests & optimization
########################################

rule install_pomona:
    input: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    output:
        pomona = "outputs/vita_rf/pomona_install.txt"
    conda: 'rf.yml'
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/vita_rf/vita_rf.RDS",
        vita_vars = "outputs/vita_rf/vita_vars.txt",
        ibd_novalidation = "outputs/vita_rf/ibd_novalidation_filt.csv",
        ibd_novalidation_diagnosis = "outputs/vita_rf/bd_novalidation_filt_diagnosis.txt",
        ibd_validation = "outputs/vita_rf/ibd_validation_filt.csv"
    conda: 'rf.yml'
    script: "scripts/vita_rf.R"

#rule tune_rf:
#    input:
#        ibd_novalidation = "outputs/vita_rf/ibd_novalidation_filt.csv",
#        ibd_novalidation_diagnosis = "outputs/vita_rf/bd_novalidation_filt_diagnosis.txt" 
#    output:
#        optimal_rf = "outputs/optimal_rf/optimal_ranger.RDS",
#        pred_test = "outputs/optimal_rf/pred_test_tab.txt",
#        pred_train = "outputs/optimal_rf/pred_train_tab.txt"
#    conda: 'rf.yml'
#    script: "scripts/tune_rf.R"
#
#rule validate_rf:
#    input:
#        optimal_rf = "outputs/optimal_rf/optimal_ranger.RDS",
#        ibd_validation = "outputs/vita_rf/ibd_validation_filt.csv",
#        info = "inputs/working_metadata.tsv" 
#    output:
#        pred_srp = "outputs/rf_validation/pred_srp057027.txt",
#        pred_srp_df = "outputs/rf_validation/pred_srp057027.csv",
#        pred_prjna = "outputs/rf_validation/pred_prjna385949.txt",
#        pred_prjna_df = "outputs/rf_validation/pred_prjna385949.csv"
#    conda: "rf.yml"
#    script: "scripts/validate_rf.R"



########################################
## Predictive hash characterization
########################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf/vita_vars.txt"
    output: "outputs/vita_rf/vita_vars.sig"
    conda: "sourmash.yml"
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name vita_vars --filename {input} {input}
    '''

rule download_gather_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/5jyzr/download
    '''

rule untar_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.sbt.json"
    input: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/3vebw/download
    '''

rule untar_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.sbt.json"
    input: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_nayfach:
    output: "inputs/gather_databases/nayfach-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/y3vwb/download
    '''

rule untar_nayfach:
    output: "inputs/gather_databases/nayfach-k31.sbt.json"
    input: "inputs/gather_databases/nayfach-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz
    '''

rule untar_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.sbt.json"
    input:  "inputs/gather_databases/genbank-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule gather_vita_vars:
    input:
        sig="outputs/vita_rf/vita_vars.sig",
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json"
    output: 
        csv="outputs/gather/vita_vars.csv",
        matches="outputs/gather/vita_vars.matches",
        un="outputs/gather/vita_vars.un"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db4} {input.db3} {input.db2}
    '''

