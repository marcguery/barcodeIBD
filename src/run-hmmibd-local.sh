#!/bin/bash

# The first argument should be the binary of hmmIBD
# The second argument should be the path to the input and output directory
# The third argument should be "TRUE" if you want hmmIBD to run in a background process

$1 &> /dev/null
if [ "$?" -eq 127 ]; then
	echo "The binary $1 is not executable."
	exit 1
fi

if [ ! -d "$2" ]; then
	echo "The directory $2 does not exist."
	exit 1
fi
if [ ! -f "$2/hmmIBD-barcodes.txt" ]; then
	echo "The file hmmIBD-barcodes.txt does not exist in $2"
	exit 1
fi

# Running hmmIBD
echo "Running hmmIBD, this could take a while..."

if [ "$3" == "TRUE" ]; then
	$1 -i "$2/hmmIBD-barcodes.txt" -o "$2/barcodes-IBD" -m 150 &
	wait %1
else
	$1 -i "$2/hmmIBD-barcodes.txt" -o "$2/barcodes-IBD" -m 150
fi

zip -j "$2/barcodes-IBD.hmm_fract.txt.zip" "$2/barcodes-IBD.hmm_fract.txt" && rm "$2/barcodes-IBD.hmm_fract.txt"
zip -j "$2/barcodes-IBD.hmm.txt.zip" "$2/barcodes-IBD.hmm.txt" && rm "$2/barcodes-IBD.hmm.txt"

echo "hmmIBD finished running."
