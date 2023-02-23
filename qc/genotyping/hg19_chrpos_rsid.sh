# Assign rsid
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/

dbsnp_hg19=~/shared-gandalm/refGenomes/hg19/dbSNP

# - Take biallelic SNPs and take the chr, pos, ID
bcftools view -m2 -M2 -v snps ${dbsnp_hg19}/00-All.vcf.gz | \
bcftools query -f '%CHROM\t%POS\t%ID\n' > \
${dbsnp_hg19}/00-All_biallelic_rsid.txt

# - Format chr and pos into chr:pos
sed 's/\t/:/' ${dbsnp_hg19}/00-All_biallelic_rsid.txt > \
${dbsnp_hg19}/00-All_coord_rsid_format_unique3.txt

# - Swap around columns, in order to sort and filter the unique entries by the final column (uniq -u -f1)
awk -F'\t' '{print $2,$1}' OFS='\t' ${dbsnp_hg19}/00-All_coord_rsid_format_unique3.txt > \
${dbsnp_hg19}/00-All_coord_rsid_format_unique3_changeorder.txt

# https://www.redhat.com/sysadmin/uniq-command-lists

# - Sort the file and filter the unique entries by the final column (uniq -u -f1)
sort -k2,2 ${dbsnp_hg19}/00-All_coord_rsid_format_unique3_changeorder.txt | \
uniq -u -f1 > \
${dbsnp_hg19}/00-All_coord_rsid_format_unique4_changeorder.txt

# - Reformat into the column order required for plink
awk -F'\t' '{print $2,$1}' OFS='\t' ${dbsnp_hg19}/00-All_coord_rsid_format_unique4_changeorder.txt > \
${dbsnp_hg19}/00-All_coord_rsid_format_unique4_changeorder_out.txt
