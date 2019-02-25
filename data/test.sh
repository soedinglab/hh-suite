#!/bin/bash -e

rm -f single* search_* blits_*

hhalign -i query.a3m -t query.a3m

ffindex_build -s single.ffdata single.ffindex query.a3m

mpirun -np 2 \
    ffindex_apply_mpi single.ffdata single.ffindex -d single_a3m_cons.ffdata -i single_a3m_cons.ffindex \
        -- hhconsensus -i stdin -oa3m stdout -M a3m -v 0

mpirun -np 2 \
    ffindex_apply_mpi single_a3m_cons.ffdata single_a3m_cons.ffindex -d single_a3m.ffdata -i single_a3m.ffindex \
        -- hhfilter -i stdin -o stdout -diff 1000 -v 0

mpirun -np 2 \
    ffindex_apply_mpi single_a3m.ffdata single_a3m.ffindex -d single_hhm.ffdata -i single_hhm.ffindex \
        -- hhmake -i stdin -o stdout -v 0

mpirun -np 2 cstranslate_mpi -i single -o single_cs219 -b -x 0.3 -c 4 -I a3m

hhblits -i query.a3m -d single -blasttab blits_app_res -n 1
hhblits_omp -i single -d single -blasttab blits_omp_res -n 1
mpirun -np 2 hhblits_mpi -i single -d single -blasttab blits_mpi_res -n 1

diff <(tr -d '\000' < blits_omp_res.ffdata) blits_app_res
diff <(tr -d '\000' < blits_mpi_res.ffdata) blits_app_res

hhsearch -i query.a3m -d single -blasttab search_app_res
hhsearch_omp -i single -d single -blasttab search_omp_res
mpirun -np 2 hhsearch_mpi -i single -d single -blasttab search_mpi_res

diff <(tr -d '\000' < search_omp_res.ffdata) search_app_res
diff <(tr -d '\000' < search_mpi_res.ffdata) search_app_res

diff <(cut -f 1-10,12 blits_app_res) <(cut -f 1-10,12 search_app_res)
