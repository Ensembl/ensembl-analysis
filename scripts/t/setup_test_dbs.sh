mysql -u eat_tester -papwd1234 -e "DROP DATABASE IF EXISTS eat_rat_core_db"
mysql -u eat_tester -papwd1234 -e "DROP DATABASE IF EXISTS eat_production_db"
mysql -u eat_tester -papwd1234 -e "CREATE DATABASE eat_production_db"

mysql -u eat_tester -papwd1234 eat_production_db < $EAT_DIR/ensembl-analysis/modules/t/test-genome-DBs/rattus_norvegicus/production/prod_tables_init.sql
