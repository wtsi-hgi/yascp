input_dir=$1
encryptor_output=$2
shopt -s nullglob
for dir in ./$input_dir/*/
do
    part2=$(basename "$dir")
    java -jar /opt/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar -i "$dir" -o $encryptor_output/"$part2"
done