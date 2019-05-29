# An example charybdis run, using the data in this charybdis/testrun folder
# See charybdis/charybdis_man.pdf for creating the needed databases
# You will need to manually set the database paths below

bash ../charybdis_generic.sh \
    -p TestData -i . -o out -n 20 \
    -x 313 -g ../bin \
    -t ../data/taxo \
    -b ../data/blastdb/nt \
    -d blast_ignore_stub.txt \
    -c chimera_db_stub.txt

