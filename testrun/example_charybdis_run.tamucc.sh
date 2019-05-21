# An example charybdis run, using the data in this charybdis/testrun folder
# See charybdis/charybdis_man.pdf for creating the needed databases

bash ../charybdis_generic.sh \
    -p TestData -i . -o out -n 20 \
    -x 313 -g ../bin \
    -t /work/hobi/GCL/db/TAXO \
    -b /work/hobi/GCL/db/NCBI_BLAST_DBs/nt \
    -d blast_ignore_stub.txt \
    -c chimera_db_stub.txt

