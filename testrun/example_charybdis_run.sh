# An example charybdis run, using the data in this charybdis/testrun folder
# Assumes that you have set up the databases in the data directory, as described in the README
# Assumes you are running with within the testrun directory

bash ../charybdis_generic.sh \
    -p TestData -i . -o out -n 20 \
    -x 313 -g ../bin \
    -t ../data/taxo \
    -b ../data/blastdb/nt \
    -v ../data/vsearchdb_coi_clean.fasta \
    -d blast_ignore_stub.txt \
    -T dc-megablast \
    -c chimera_db_stub.txt

