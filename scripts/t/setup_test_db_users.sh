mysql -u root -p%1 -e "CREATE USER 'eat_tester'@'localhost' IDENTIFIED BY 'apwd1234'"
mysql -u root -p%1 -e "GRANT ALL PRIVILEGES ON * . * TO 'eat_tester'@'localhost'"
