IN="$1"
IFS=' ' read -ra ADDR <<< "$IN"
mails=$(echo $IN | tr " " "\n")
for addr in $mails
do
    echo $addr >> ./fofn_vcfs.txt
done