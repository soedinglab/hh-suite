#!/bin/sh -e
INPUT="$1"
OUTPUT="$2"
THREADS=""
if [ -n "$3" ]; then
  THREADS="--threads $3"
fi

MPI="$4"
#command -v "$(echo "$MPI" | cut -d " " -f1)" >/dev/null 2>&1
#if [ "$?" -eq "1" ]; then
#  echo "MPI will not be used"
#  unset MPI
#fi

if [ "$#" -lt 2 ]; then
  echo "createdb.sh input output [threads] '[mpirun --mpi-params]'"
  echo "Create HH-suite datebase from MMseqs2 input"
  echo "Example:"
  echo "  conda create -n hhdb mmseqs2 famsa hhsuite"
  echo "  conda activate hhdb"
  echo "  mmseqs createdb input.fas input"
  echo "  mmseqs cluster input clu tmp"
  echo "  mmseqs createseqfiledb input clu seqfiledb"
  echo "  createdb.sh seqfiledb hhdatabase"
  exit 1
fi

> "make_msa.sh" cat <<EOF
#!/bin/sh -e
famsa  -t 1 STDIN STDOUT | hhconsensus -maxres 65535 -i stdin -o stdout -v 0 2> /dev/null
EOF

chmod +x make_msa.sh

echo "calculate alignments"
$MPI mmseqs apply "${INPUT}" "${OUTPUT}_a3m" $THREADS -- "./make_msa.sh" 2> /dev/null
rm -f "make_msa.sh"

echo "calculate hidden markov models of large MSAs"
$MPI mmseqs apply "${OUTPUT}_a3m" "${OUTPUT}_hhm_sizes" $THREADS -- awk '/^>/{cnt++;} cnt>51 {print "1"; exit }' 
awk 'FNR==NR && $3 > 1 {f[$1]=1; next} $1 in f {print}' "${OUTPUT}_hhm_sizes.index" "${OUTPUT}_a3m.index" > "${OUTPUT}_a3m_large.index"
ln -fs "${OUTPUT}_a3m" "${OUTPUT}_a3m_large"
$MPI mmseqs apply "${OUTPUT}_a3m_large" "${OUTPUT}_hhm" $THREADS -- hhmake -i stdin -o stdout -v 0
rm -f "${OUTPUT}_a3m_large" "${OUTPUT}_a3m_large.index" "${OUTPUT}_hhm_sizes.ffindex"

echo "calculate context states"
ln -fs "${OUTPUT}_a3m" "${OUTPUT}_a3m.ffdata"
ln -fs "${OUTPUT}_a3m.index" "${OUTPUT}_a3m.ffindex"
if [ "$MPI" = "" ]; then
  cstranslate -x 0.3 -c 4 --ffindex -I a3m -i "${OUTPUT}_a3m" -o "${OUTPUT}_cs219" > /dev/null
else
  $MPI cstranslate_mpi -x 0.3 -c 4 -I a3m -i "${OUTPUT}_a3m" -o "${OUTPUT}_cs219" 2> /dev/null
fi
ln -f "${OUTPUT}_cs219.ffdata" "${OUTPUT}_cs219"
ln -f "${OUTPUT}_cs219.ffindex" "${OUTPUT}_cs219.index"

echo "reorder databases for faster access"
sort -k 3 -n "${OUTPUT}_cs219.index" | cut -f1 > "${OUTPUT}_sorted.tsv"
for type in a3m hhm; do
  mmseqs createsubdb "${OUTPUT}_sorted.tsv" "${OUTPUT}_${type}" "${OUTPUT}_${type}_opt" -v 0
  mv -f "${OUTPUT}_${type}_opt" "${OUTPUT}_${type}.ffdata"
  LC_ALL=C sort -k1,1 "${OUTPUT}_${type}_opt.index" > "${OUTPUT}_${type}.ffindex"
  rm -f "${OUTPUT}_${type}_opt.index"
done
rm -f "${OUTPUT}_sorted.tsv"
