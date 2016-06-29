 #!/bin/bash

display_usage()
{
	echo "The user will need GRANT and DROP privileges on the specified host."
	echo -e "\nUsage:\n$0 username password hostname portnumber read_only_username\n"
}


# if less than five arguments supplied, display usage
if [  $# -le 4 ]
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


farm_test_dir=$EAT_DIR'/ensembl-analysis/modules/t/farm'

#read in values
db_user=$1  # for EAT_DB_USER
db_pass=$2  # for EAT_DB_PASS
db_host=$3  # for EAT_DB_HOST
db_port=$4  # for EAT_DB_PORT
db_user_ro=$5  # for EAT_DB_USER_RO

db_name=$USER'_eat_rat_core_db'  # for EAT_DB_NAME

#get a list of template files and update
echo Working in $farm_test_dir
for f in  $farm_test_dir/*.template ;
do
 filename=${f##*/}
 output_file=${filename%.*}
 echo Using template $filename to generate $output_file
 gawk '{ gsub("EAT_DB_USER", "'$db_user'"); print }' $f | \
     gawk '{ gsub("EAT_DB_PASS", "'$db_pass'"); print }' | \
     gawk '{ gsub("EAT_DB_PORT", "'$db_port'"); print }' | \
     gawk '{ gsub("EAT_DB_NAME", "'$db_name'"); print }' | \
     gawk '{ gsub("EAT_DB_USER", "'$db_user_ro'"); print }' | \
     gawk '{ gsub("EAT_DB_HOST", "'$db_host'"); print }' > $farm_test_dir/$output_file
done

echo Working in $farm_test_dir/conf
for f in  $farm_test_dir/conf/*.template ;
do
 filename=${f##*/}
 output_file=${filename%.*}
 echo Using template $filename to generate $output_file
done
