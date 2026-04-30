  set -e   
  # expects the target directory to be passed in in which to find the .pbs file 
  
target="$1"

if [ -z "$target" ]; then
    echo "Usage: $0 <target_dir>"
    exit 1
fi

# Check files exist
for f in genray.in genray_ece.in genr_sam.pbs; do
    if [ ! -f "$target/$f" ]; then
        echo "Missing file: $target/$f"
        exit 1
    fi
done

# Swap genray.in and genray_ece.in using genray_LH.in as temporary
mv "$target/genray.in"     "$target/genray_LH.in"
mv "$target/genray_ece.in" "$target/genray.in"

# Submit job
genray_id=$(sbatch --parsable "$target/genr_yuri.pbs")
echo "Submitted job $genray_id"

# Restore originals
#mv "$target/genray.in"     "$target/genray_ece.in"
#mv "$target/genray_LH.in"  "$target/genray.in"