#! /usr/bin/python 

#### Moducles, variables, paths, ...
import os,sys,pdb,re,subprocess
import pandas as pd

from bx.intervals.cluster import ClusterTree


# Your working directory
path = os.getcwd() # Path to directory where input file lives and all output files are written


####### GWAS summary statistics file settings

label = "EA2_EduYears_pooled_single_gc.meta" # Enter the filename of your gwas summary statistics file (without the file extension)
gwas_filename = "%s"%(label) 
independent_snps_filename = "%s_independent_snps.tab"%(label) 
independent_loci_filename = "%s_independent_loci.tab"%(label) 
pvalue_col = 10 # NB: Counting starts from 0, ie. first columns is referred to as '0'
marker_col = 1 # Format: <chr:pos>, ie. '6:2321'.  If this column does not exist chr_col and pos_col will be used and this column should be set to 'None'
chr_col = None # Set ot 'None' if this column does not exist
pos_col = None # Set ot 'None'if this column does not exist
plink_clumping_snp_column = "rsID"
plink_clumping_pvalue_column = "Pvalue"


####### PLINK settings

plink_genotype_data_plink_prefix = "1000_genomes_project_phase1_EUR/CEU_GBR_TSI_unrelated.phase1_dup_excluded" # Change path the the prefix of you 1000 Genomes data
plink_binary = "%s/plink_mac/plink"%path # Change to your PLINK executable
plink_clumping_pvalue =  5e-8
plink_clumping_distance = 500
plink_clumping_r2 = 0.1
plink_ld_snp_list_file = "{}/{}.snps".format(path,label)
plink_ld_window_kb = 1000
plink_ld_window = 99999
plink_ld_window_r2 = 0.2
plink_output_prefix = "{}/{}_ldpartners".format(path,label)


####### Locus merging settings
merging_distance_kb = 250000
r2_threshold = 0.8

# SNPsnap collection used to define locus boundaries and genes in associated loci
collection_file = "%s/ld0.6_collection.tab.gz"%path # Change path the SNPsnap collection (Download from http://www.broadinstitute.org/mpg/snpsnap/database_download.html)

# Gene identifier mapping file
ensembl_identifier_mapping_file = "mapping_ENSG_HGNC_Entrez_v73.tab.gz"


############################################################################################
#### Helper functions for running PLINK

####### Run PLINK

def run_plink_clumping(plink_genotype_data_plink_prefix, plink_binary, plink_clumping_pvalue, plink_clumping_distance, plink_clumping_r2, plink_clumping_snp_column, plink_clumping_pvalue_column, path, gwas_filename, out_dir):
	"""
    	Function that calls PLINK to get clump GWAS SNPs
	"""
	cmd = [plink_binary, 
		"--bfile", plink_genotype_data_plink_prefix,
		"--clump-p1", str(plink_clumping_pvalue),
		"--clump-kb", str(plink_clumping_distance),
		"--clump-r2", str(plink_clumping_r2),
		"--clump-snp-field", plink_clumping_snp_column,
		"--clump-field", plink_clumping_pvalue_column,
		"--clump", path+"/"+gwas_filename,
		"--out", out_dir+"/"+gwas_filename
	]
	print "making plink pruning call {}".format( " ".join(cmd)  )
	return subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

def run_plink_ld(plink_binary, plink_genotype_data_plink_prefix, plink_ld_snp_list_file, plink_ld_window_kb, plink_ld_window, plink_ld_window_r2, path, gwas_filename, plink_output_prefix):
    """
    Function that calls PLINK to get LD buddies
    """
    cmd = [plink_binary, 
        "--bfile", plink_genotype_data_plink_prefix,
        "--r2",
        "--ld-snp-list",  str(plink_ld_snp_list_file),
        "--ld-window-kb", str(plink_ld_window_kb), 
        "--ld-window", str(plink_ld_window),
        "--ld-window-r2", str(plink_ld_window_r2),
        "--out", plink_output_prefix
    ]
    
    print "making plink ld call {}".format( " ".join(cmd)  )
    
    return subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()


####### Read PLINK results

def get_plink_clumping(path,label):    
    df = pd.read_csv("{}/{}.clumped".format(path,label),delimiter=r"\s+")
    return df


############################################################################################
#### Run the analysis

####### Define linkage-disequilibrium (LD)-independent SNPs as those with low LD (r2 < 0.1) using a window of 500-kb (using PLINK's clumping option).

(std, err) = run_plink_clumping(plink_genotype_data_plink_prefix, plink_binary, plink_clumping_pvalue, plink_clumping_distance, plink_clumping_r2, plink_clumping_snp_column, plink_clumping_pvalue_column, path, gwas_filename, path)
index_snps_df = get_plink_clumping(path,label)
index_snps_df.set_index(index_snps_df.CHR.astype(str) + ":" + index_snps_df.BP.astype(str),inplace=True)

###### Add new column with clumped SNPs sorted by LD
pd.Series(index_snps_df.SNP).to_csv(plink_ld_snp_list_file, index=False, header=False)
(std, err) = run_plink_ld(plink_binary, plink_genotype_data_plink_prefix, plink_ld_snp_list_file, plink_ld_window_kb, plink_ld_window, plink_ld_window_r2, path, gwas_filename, plink_output_prefix)
ld_partners_df = pd.read_csv("{}.ld".format(plink_output_prefix),delimiter=r"\s+")
def get_ld_partners(row):
    snp_partners_df = ld_partners_df.ix[ld_partners_df.SNP_A==row.SNP,:].sort(['R2'], ascending=[0])
    clump = row.SP2.replace("(1)","").split(",")
    return (",".join(["{} ({})".format(snp_b,r2) for (snp_b,r2) in zip(snp_partners_df.SNP_B,snp_partners_df.R2.map(lambda x: "{0:.2f}".format(x))) if row.SNP != snp_b and snp_b in clump]))
index_snps_df['SP2_ld'] = index_snps_df.apply(get_ld_partners,axis=1)


# ###### Define associated loci (upstream and downstream boundaries of SNPs from step #1). We use the precomputed SNPsnap collections to directly extract genes and locus boundaries
# * Get boundaries (LD r2>0.6)
# * Get genes (to be added)

collection = pd.read_csv(collection_file, index_col=0, header=0, delimiter="\t", compression = 'gzip')
results_df = index_snps_df.join(collection,how="left")
results_snps_df = results_df.ix[:,['SNP','CHR','BP','P','SP2_ld','loci_upstream','loci_downstream']]
results_snps_df.rename(columns={'SNP': 'snp_name', 'CHR': 'chr', 'BP': 'pos','P': 'pvalue', 'SP2_ld': 'plink_ld_partners', 'loci_upstream': 'locus_upstream_boundary', 'loci_downstream': 'locus_downstream_boundary'}, inplace=True)
results_snps_df['proxy_snp'] = None


####### Find LD-partners of SNPs that are not found in collection

missing_from_collection = results_snps_df.ix[results_snps_df.locus_upstream_boundary.isnull(),:]

def find_ld_partner(ld_partner_str, column_id):
    for (rs_id, r2) in [ ( x.split("(")[0], float(x.split("(")[1].replace(")","") ) ) for x in ld_partner_str.split(',')]:
        snp_index = collection.index[collection.rsID == rs_id]
        if r2 > r2_threshold and len(snp_index):
            return collection.ix[snp_index,column_id].values[0]
    return None


missing_from_collection['proxy_snp'] = missing_from_collection.plink_ld_partners.apply(find_ld_partner,args=('rsID',))
missing_from_collection['locus_upstream_boundary'] = missing_from_collection.plink_ld_partners.apply(find_ld_partner,args=('loci_upstream',))
missing_from_collection['locus_downstream_boundary'] = missing_from_collection.plink_ld_partners.apply(find_ld_partner,args=('loci_downstream',))
results_snps_df.update(missing_from_collection)
results_snps_df['chr'] = results_snps_df['chr'].astype('int')
results_snps_df['pos'] = results_snps_df['pos'].astype('int')
results_snps_df['locus_downstream_boundary'] = results_snps_df['locus_downstream_boundary'].astype('int')
results_snps_df['locus_upstream_boundary'] = results_snps_df['locus_upstream_boundary'].astype('int')


####### Merge associated loci within 250 kb of each other.

trees = {}
min_intervals = 0
for i in range(1, 23, 1):
    trees[i] = ClusterTree(merging_distance_kb,min_intervals)


for i, (index, row) in enumerate(results_snps_df.iterrows()):
    if row.chr in range(1,23,1):
        trees[row.chr].insert(row.locus_upstream_boundary,row.locus_downstream_boundary,i)


results_loci_df = pd.DataFrame(columns=['snp_name','chr','pos','pvalue','locus_upstream_boundary','locus_downstream_boundary'])
results_snps_df['locus'] = None
counter = 0
for chrom in trees:
    for (start, end, loci) in trees[chrom].getregions():
        
        # Single SNP locus
        if len(loci) == 1:
            results_snps_df.ix[loci[0],'locus'] = counter
            results_loci_df.loc[counter] = results_snps_df.ix[loci[0],['snp_name','chr','pos','pvalue','locus_upstream_boundary','locus_downstream_boundary']]

        # Multi SNP locus
        else:    
            min_start = float('-inf')
            max_end  = float('inf')
            min_pvalue = 1
            min_snp_name = None
            min_pos = None
            for locus in loci:
                results_snps_df.ix[locus,'locus'] = counter
                min_start = start if start < min_start else None
                max_end = end if end > max_end else None
                if results_snps_df.ix[locus,'pvalue'] < min_pvalue:
                    min_pvalue = results_snps_df.ix[locus,'pvalue']
                    min_pos = results_snps_df.ix[locus,'pos']
                    min_snp_name = results_snps_df.ix[locus,'snp_name']
            results_loci_df.loc[counter] = pd.Series(data=[min_snp_name, chrom, min_pos, min_pvalue, start, end], index = ['snp_name','chr','pos','pvalue','locus_upstream_boundary','locus_downstream_boundary'])
        
        counter += 1
results_loci_df.set_index(results_loci_df.chr.astype('int').astype('str') + ":" + results_loci_df.pos.astype('int').astype('str'),inplace=True)


# Make sure types are correct
results_loci_df['chr'] = results_loci_df['chr'].astype('int')
results_loci_df['pos'] = results_loci_df['pos'].astype('int')
results_loci_df['locus_downstream_boundary'] = results_loci_df['locus_downstream_boundary'].astype('int')
results_loci_df['locus_upstream_boundary'] = results_loci_df['locus_upstream_boundary'].astype('int')


# Add genes in locus
results_loci_genes_df = results_loci_df.join(collection.ix[results_loci_df.index,['ID_genes_in_matched_locus']],how="left")
results_loci_genes_df.rename(columns={'ID_genes_in_matched_locus': 'genes_in_locus_ensembl'}, inplace=True)
ensembl_raw_df = pd.read_csv(ensembl_identifier_mapping_file, index_col=0, sep="\t", compression='gzip')
grouped = ensembl_raw_df.groupby(level=0)
ensembl_df = grouped.last()

def ensembl_to_hgnc(ensembl_gene_list):
    hgnc_list = []
    if not isinstance(ensembl_gene_list, float):
        for e in ensembl_gene_list.split(";"):
            if e in ensembl_df.index and not isinstance(ensembl_df.ix[e,'HGNC symbol'],float):
                hgnc_list.append(ensembl_df.ix[e,'HGNC symbol'])
            else:
                hgnc_list.append(e)

    return ";".join(hgnc_list) if len(hgnc_list) > 0 else None


results_loci_genes_df['genes_in_locus_hgnc'] = results_loci_genes_df.genes_in_locus_ensembl.apply(ensembl_to_hgnc)



############################################################################################
####### Save to output files

results_snps_df.sort(['locus','pos'],inplace=True)
results_snps_df.to_csv(independent_snps_filename, header=True, sep="\t",index=False, na_rep='NA',float_format="%10.2e")
results_loci_genes_df.to_csv(independent_loci_filename, header=True, sep="\t",index=False, na_rep='NA', float_format="%10.2e")

