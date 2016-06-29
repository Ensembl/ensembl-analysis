 #!/bin/bash

display_usage()
{
	echo "The user will need GRANT and DROP privileges on the specified host."
	echo -e "\nUsage:\n$0 username password hostname portnumber\n"
}


# if less than two arguments supplied, display usage
if [  $# -le 1 ]
then
    display_usage
    exit 1
fi

# check whether user had supplied -h or --help . If yes display usage
if [[ ( $# == "--help") ||  $# == "-h" ]]
then
    display_usage
    exit 0
fi

dbname=$USER'_eat_rat_core_db'

echo Dropping, then creating, $dbname  on MySQL host $3 as user $1

mysql -u $1 -p$2 -h $3 -P $4 -e "DROP DATABASE IF EXISTS $dbname"
mysql -u $1 -p$2 -h $3 -P $4 -e "CREATE DATABASE $dbname"
mysql -u $1 -p$2 -h $3 -P $4 $dbname < $EAT_DIR/ensembl-analysis/modules/t/farm/test-genome-DBs/pre_peats.sql

exit 0
