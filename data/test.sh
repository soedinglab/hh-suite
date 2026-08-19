#!/bin/bash -e

rm -f single* search_* blits_*

if $(command -v "ffindex_apply_mpi" >/dev/null 2>&1); then
    MPI=1
else
    unset MPI
fi

hhalign -i query.a3m -t query.a3m

ffindex_build -s single.ffdata single.ffindex query.a3m

APPLY="ffindex_apply"
if [ -n "$MPI" ]; then
    APPLY="mpirun -np 2 ffindex_apply_mpi"
fi

$APPLY single.ffdata single.ffindex -d single_a3m_cons.ffdata -i single_a3m_cons.ffindex \
        -- hhconsensus -i stdin -oa3m stdout -M a3m -v 0

$APPLY single_a3m_cons.ffdata single_a3m_cons.ffindex -d single_a3m.ffdata -i single_a3m.ffindex \
        -- hhfilter -i stdin -o stdout -diff 1000 -v 0

$APPLY single_a3m.ffdata single_a3m.ffindex -d single_hhm.ffdata -i single_hhm.ffindex \
        -- hhmake -i stdin -o stdout -v 0

if [ -n "$MPI" ]; then
    mpirun -np 2 cstranslate_mpi -i single -o single_cs219 -b -x 0.3 -c 4 -I a3m
else
    cstranslate -i single -o single_cs219 -b -x 0.3 -c 4 -I a3m -f
fi

hhblits -i query.a3m -d single -blasttab blits_app_res -n 1
hhblits_omp -i single -d single -blasttab blits_omp_res -n 1
diff <(tr -d '\000' < blits_omp_res.ffdata) blits_app_res
if [ -n "$MPI" ]; then
    mpirun -np 2 hhblits_mpi -i single -d single -blasttab blits_mpi_res -n 1
    diff <(tr -d '\000' < blits_mpi_res.ffdata) blits_app_res
fi

hhsearch -i query.a3m -d single -blasttab search_app_res
hhsearch_omp -i single -d single -blasttab search_omp_res
diff <(tr -d '\000' < search_omp_res.ffdata) search_app_res

if [ -n "$MPI" ]; then
    mpirun -np 2 hhsearch_mpi -i single -d single -blasttab search_mpi_res
    diff <(tr -d '\000' < search_mpi_res.ffdata) search_app_res
fi

diff <(cut -f 1-10,12 blits_app_res) <(cut -f 1-10,12 search_app_res)

# Entry names are scoped to a database. Verify that MAC realignment does not
# combine identically named entries from two databases and reuse the wrong HMM.
MULTIDB_TEST_DIR=$(mktemp -d "${TMPDIR:-/tmp}/hhsuite-multidb.XXXXXX")
trap 'rm -rf "$MULTIDB_TEST_DIR"' EXIT

mkdir -p "$MULTIDB_TEST_DIR"/{db1,db2}/{a3m,hhm,cs219}
awk 'NR <= 2' query.a3m > "$MULTIDB_TEST_DIR/query.a3m"
awk 'NR == 1 { print ">db1_target" } NR == 2 { print }' query.a3m \
    > "$MULTIDB_TEST_DIR/db1/a3m/collision"
awk 'NR == 1 { print ">db2_target" } NR == 2 { print substr($0, 1, 300) }' query.a3m \
    > "$MULTIDB_TEST_DIR/db2/a3m/collision"

for db in db1 db2; do
    hhmake -i "$MULTIDB_TEST_DIR/$db/a3m/collision" \
        -o "$MULTIDB_TEST_DIR/$db/hhm/collision" -v 0
    cstranslate -i "$MULTIDB_TEST_DIR/$db/a3m/collision" \
        -o "$MULTIDB_TEST_DIR/$db/cs219/collision" -I a3m -x 0.3 -c 4

    for kind in a3m hhm cs219; do
        ffindex_build -s "$MULTIDB_TEST_DIR/${db}_${kind}.ffdata" \
            "$MULTIDB_TEST_DIR/${db}_${kind}.ffindex" \
            "$MULTIDB_TEST_DIR/$db/$kind"
    done
done

hhblits -i "$MULTIDB_TEST_DIR/query.a3m" \
    -d "$MULTIDB_TEST_DIR/db1" -d "$MULTIDB_TEST_DIR/db2" \
    -oa3m "$MULTIDB_TEST_DIR/result.a3m" -o "$MULTIDB_TEST_DIR/result.hhr" \
    -n 1 -e 10 -maxfilt 100 -realign_max 100 -min_prefilter_hits 1 \
    -premerge 0 -cpu 1 -v 1

test -s "$MULTIDB_TEST_DIR/result.a3m"
awk '$2 == "db2_target" && $8 == 300 { found = 1 } END { exit !found }' \
    "$MULTIDB_TEST_DIR/result.hhr"
