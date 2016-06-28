 #!/bin/bash

#TODO - add port option
#TODO - add argument checks

mysql -u $1 -p$2 -h $3 -e "DROP DATABASE IF EXISTS $USER_eat_rat_core_db"
mysql -u $1 -p$2 -h $3 -e "CREATE DATABASE $USER_eat_rat_core_db"
mysql -u $1 -p$2 -h $3 $USER_eat_rat_core_db < $EAT_DIR/ensembl-analysis/modules/t/farm/test-genome-DBs/pre_peats.sql
