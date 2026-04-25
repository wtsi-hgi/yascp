vireo_gt_vcf=$1
out_dir=$2
sorted_vcf=$out_dir/sorted_vcf.vcf.gz
vireo_fixed_vcf=$out_dir/headfix_srt_vcf.vcf.gz
# fix header of vireo VCF
bcftools view -h $vireo_gt_vcf > init_head.txt
sed -i '/^##fileformat=VCFv.*/a ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' init_head.txt
head -n -1 init_head.txt > header.txt
echo '##INFO=<ID=AD,Number=A,Type=Integer,Description="alternative allele  (variant-by-cell) of reads">' >> header.txt
echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
echo '##INFO=<ID=PL,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
echo '##INFO=<ID=OTH,Number=1,Type=Integer,Description="????">' >> header.txt
echo '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="???">' >> header.txt
echo '##FORMAT=<ID=AD,Number=G,Type=Integer,Description="????n">' >> header.txt
echo '##FORMAT=<ID=DP,Number=G,Type=Integer,Description="????n">' >> header.txt
tail -n1 init_head.txt >> header.txt

# sort VCF file (bcftools sort bails out with an error)

bcftools view $vireo_gt_vcf |     awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' |     bcftools view -Oz -o $sorted_vcf -

bcftools reheader -h header.txt $sorted_vcf | \
bcftools view -Oz -o $vireo_fixed_vcf
bcftools index $vireo_fixed_vcf
rm header.txt
rm init_head.txt