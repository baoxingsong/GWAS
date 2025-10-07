## Prepare the genotype files
```sh
plink --vcf hwc358.snp.chr.mis9_maf5.bg.vcf.gz --maf 0.05 --make-bed --out snp   #filter MAF and reformat
mv snp.bim snp.bim.backup
cat snp.bim.backup | awk '{print($1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6)}' > snp.bim   #give each SNP an ID
##go to the matrix folder to filter genomicly super similar individuls
plink --bfile ../snp  --geno 0.2 --maf 0.05 --out snp --make-bed   # fillter by MAF and missing
plink --bfile snp --freq --missing  # this and the another two lines, filter MAF and missing again by count. Parameters could be set via the source code of perl excludeList.pl
perl excludeList.pl
plink --bfile snp --exclude excludeSnps --make-bed --out snp_gwas
plink --bfile snp_gwas --indep-pairwise 10000 5000 0.90  ### the line and the next line: filter SNPs via LD
plink --bfile snp_gwas --extract plink.prune.in --make-bed --out snp_matrix   

#for matrix plot begin
plink --bfile snp_matrix --cluster --matrix #make IBS kinship matrix
Rscript	plotMatrix.R   # plot IBS matrix
#for matrix plot end

perl filterIndividual.pl | sort | uniq > ../remove   # generate list with geneomicly duplicated individuls, parameters could be set via filterIndividual.pl

cd ../ 

plink --vcf hwc358.snp.chr.mis9_maf5.bg.vcf.gz --remove remove --maf 0.05 --geno 0.1 --make-bed --out hwc.snp.chr.mis9_maf5.ibs  # filter duplicated individuls and MAF and missiong ratio,  please pay attention to the input file, is the original VCF file.

cp hwc.snp.chr.mis9_maf5.ibs.fam hwc.snp.chr.mis9_maf5.ibs.fam.backup
cat hwc.snp.chr.mis9_maf5.ibs.fam.backup | awk '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t10")}' > hwc.snp.chr.mis9_maf5.ibs.fam    # fill a non-sense phenotype
```

## GWAS
```sh
gemma-0.98.5 -bfile hwc.snp.chr.mis9_maf5.ibs -gk 1 -o hwc.snp.chr.mis9_maf5.ibs  # prepare kinship matrix using gemma for GWAS analysis
gemma-0.98.5 -bfile hwc.snp.chr.mis9_maf5.ibs -gk 2 -o hwc.snp.chr.mis9_maf5.ibs  # prepare kinship matrix, using a different method implemented in gemma, for GWAS analysis

# prepare phenotype begin   # highlikely, this part should be reimplemented
perl preparePhenotype.pl > hwc.snp.chr.mis9_maf5.ibs.fam.t  # likely need to reimplement this part
cat hwc.snp.chr.mis9_maf5.ibs.fam.backup | awk '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t10")}' > hwc.snp.chr.mis9_maf5.ibs.fam
cat hwc.snp.chr.mis9_maf5.ibs.fam.t | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > hwc.snp.chr.mis9_maf5.ibs.fam
# prepare phenotype end

gemma-0.98.5 -bfile ./hwc.snp.chr.mis9_maf5.ibs -k ./output/hwc.snp.chr.mis9_maf5.ibs.cXX.txt -lmm 4 -o FW_C.result



```

## Plot the Manhattan and QQ plot
```
cd output
Rscript resultVirtulization_emmax.R FW_C.result.assoc.txt    # There are four methods to calculate p-values, here I only ploted one of them. for input files et al. parameters should be set via the source code

```