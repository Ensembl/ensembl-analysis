#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Connection to the target DB
my $host = 'mysql-rfam-public.ebi.ac.uk';
my $port = '4497';
my $user = 'rfamro';
my $pass;
my $dbname = 'Rfam';
my $output_dir;
my $clade;

&GetOptions (
            'h|host|dbhost=s' => \$host,
            'P|port|dbport=s' => \$port,
            'u|user|dbuser=s' => \$user,
            'p|pass|dbpass=s' => \$pass,
            'd|dbname=s'      => \$dbname,
            'o|output_dir=s'  => \$output_dir,
            'c|clade=s'       => \$clade,
        );

my $query = "
SELECT 
  rfam_acc 
FROM 
  (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff 
FROM 
  full_region fr, rfamseq rf, taxonomy tx, family f 
WHERE 
  rf.ncbi_id = tx.ncbi_id AND 
  f.rfam_acc = fr.rfam_acc AND 
  fr.rfamseq_acc = rf.rfamseq_acc AND 
  LOWER(tx.tax_string) LIKE '%$clade%' AND 
  (f.type LIKE '%snRNA%' OR 
    f.type LIKE '%rRNA%' OR 
    LOWER(f.rfam_id) LIKE '%rnase%' OR 
    LOWER(f.rfam_id) LIKE '%vault%' OR 
    LOWER(f.rfam_id) LIKE '%y_rna%' OR 
    f.rfam_id LIKE '%Metazoa_SRP%') AND 
    is_significant = 1) AS TEMP
  WHERE
    rfam_id NOT LIKE '%bacteria%' AND 
    rfam_id NOT LIKE '%archaea%' AND 
    rfam_id NOT LIKE '%microsporidia%'
";

my $command = "mysql -u $user -h $host -P $port -D $dbname -NB -e \" $query \" > " . $output_dir . "/accessions.txt";

system($command);

