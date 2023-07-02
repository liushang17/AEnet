Rscript ../bin/01.merge.sj.r --sj_dir ../resource/SJ.out.tab --outdir ../result/03.DESJ --min_read 10 --min_cell 4 --cpu 4
perl ../bin/02.Junction.ann.pl ../resource/gencode.v27.annotation.gtf ../result/03.DESJ/Alljunction.filter.list.xls ../resource/SJ.out.tab ../result/03.DESJ
Rscript ../bin/03.normalize.r ../result/03.DESJ/merge.count.txt ../resource/Cellinfo.xls ../result/03.DESJ
Rscript ../bin/04.splicing.pattern.v3.r ../result/03.DESJ/Normalize.step3.rds ../result/03.DESJ/Alljunction.filter.list.ann.onegene.xls ../result/04.AEN/gene_pattern
perl ../bin/05.merge.pattern1.pl ../result/04.AEN/gene_pattern ../result/04.AEN/Cell.chose.isoform.xls
Rscript ../bin/06.trans.xls2rds.r ../result/04.AEN/Cell.chose.isoform.xls
Rscript ../bin/07.DE.pattern.v2.r ../resource/Neu.exp.newname.rds ../result/04.AEN/Cell.chose.isoform.rds ../result/04.AEN/DE_pattern 1 4800
perl ../bin/08.merge.DE.all.v2.pl ../result/04.AEN/DE_pattern > ../result/04.AEN/AS_DE.all.xls
Rscript ../bin/09.ASP_module.phenotype.r ../result/04.AEN/AS_DE.all.xls 30 20 30 30 ../result/04.AEN
Rscript ../bin/10.downstream.R ../result/04.AEN/AS_DE.all.xls ../resource/splicing.factor.csv ../result/04.AEN/Phenotype.v1.xls ../result/04.AEN/ASP.v1.xls ../result/04.AEN