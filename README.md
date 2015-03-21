# gwas-snps-loci.py
Script to compute independent SNPs and loci based on GWAS data.

## Dependencies
* Mac OS X, or UNIX operating system (Microsoft Windows is not supported)
* Python version 2.7 (or higher)
  * [Python.org](https://www.python.org/downloads/)
* PIP (used for install Python libraries)
  * `sudo easy_install pip` 
* Python-bx (on Mac OS X you may be prompted to install XCode)
  * `sudo pip install bx-python`   
* Pandas (version 0.15.2 or higher)
  * `sudo pip install pandas`
* PLINK (only needed if you want to construct loci yourself instead of using the precomputed onces for this example)
  * [PLINK version 1](http://pngu.mgh.harvard.edu/~purcell/plink/) or [PLINK version 2](https://www.cog-genomics.org/plink2/) 
 * Make sure you have PLINK binary genotype formated 1000 Genomes Project genotype data.  You can download [preformated data here](http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz).
 * Please download the SNPsnap collection file, which provides you with precomputed LD r2 boundaries for each 1000 Genome Project phase 3 SNPs. Files can be [downloaded here](http://www.broadinstitute.org/mpg/snpsnap/database_download.html).

## Quick start
  1. Set `label` to the name of your GWAS summary statistics file (leave out the file extension)
  2. Set `plink_genotype_data_plink_prefix` to the path and prefix of you 1000 Genomes Project genotype data.
  3. Set `collection_file` to the filename of the SNPsnap collection file.
  4. Run the script.  The default settings will compute
    * Indpendent SNPs, that is SNPs with low linkage disequilibrium (r2 < 0.1) to a more significantly associated SNP within a 500-kb window.
    * Indpendent loci, that is loci containing all SNPs correlated at r2 > 0.6 with any other associated SNP.  Associated loci closer to 250 kb to each other are merged. 
  
