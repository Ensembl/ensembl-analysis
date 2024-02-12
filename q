[1mdiff --git a/modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm b/modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm[m
[1mindex 3ade44042..32dc55db4 100644[m
[1m--- a/modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm[m
[1m+++ b/modules/Bio/EnsEMBL/Analysis/Hive/Config/LayerAnnotationStatic.pm[m
[36m@@ -473,6 +473,8 @@[m [msub _master_config {[m
 	  {[m
 	      ID         => 'LAYER3',[m
 	      BIOTYPES   => [[m
[32m+[m		[32m  'self_pe12_tr_3',[m
[32m+[m		[32m  'self_pe12_tr_4'[m
 		  'human_pe12_sp_3',[m
 		  'human_pe12_sp_4',[m
 		  'fish_pe12_sp_3',[m
[36m@@ -480,7 +482,9 @@[m [msub _master_config {[m
 		  'vert_pe12_sp_1',[m
                   'vert_pe12_sp_2',[m
 		  'fish_pe12_tr_3',[m
[32m+[m		[32m  'human_pe12_tr_3'[m
 		  'genblast_rnaseq_medium',[m
[32m+[m		[32m  'rnaseq_tissue_5',[m
 		  ],[m
 	      FILTER_AGAINST => ['LAYER1','LAYER2'],[m
 	      DISCARD    => 0,[m
[36m@@ -496,14 +500,24 @@[m [msub _master_config {[m
 	      FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],[m
 	      DISCARD    => 0,[m
 	  },[m
[31m-[m
[32m+[m[32m          {[m
[32m+[m[32m              ID         => 'LAYER5',[m
[32m+[m[32m              BIOTYPES   => [[m
[32m+[m		[32m  'vert_pe12_tr_1',[m
[32m+[m		[32m  'vert_pe12_tr_2',[m
[32m+[m		[32m  'vert_pe12_tr_3',[m
[32m+[m		[32m  'rnaseq_tissue_6',[m
[32m+[m[32m                  ],[m
[32m+[m[32m              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],[m
[32m+[m[32m              DISCARD    => 0,[m
[32m+[m[32m          },[m
 	  {[m
[31m-	      ID         => 'LAYER5',[m
[32m+[m	[32m      ID         => 'LAYER6',[m
 	      BIOTYPES   => [[m
 		  'cdna',[m
 		  'rnaseq_tissue',[m
 		  ],[m
[31m-	      FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],[m
[32m+[m	[32m      FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],[m
 	      DISCARD    => 0,[m
 	  },[m
 [m
