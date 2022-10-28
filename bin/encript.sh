input_dir=$1
encryptor_output=$2
version="2.0.0"
shopt -s nullglob
for dir in ./$input_dir/*/
do
    part2=$(basename "$dir")
    java -jar /opt/EGA-Cryptor-$version/ega-cryptor-$version.jar -i "$dir" -o $encryptor_output/"$part2"
done
return $version