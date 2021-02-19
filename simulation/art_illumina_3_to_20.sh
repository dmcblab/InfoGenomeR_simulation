#!/bin/bash
reference="/DASstorage6/leeyh/HG19_Broad_variant_for_simulation/hg19.fa"

j=1;
	mkdir set5
	cd set5
	art_illumina -d 5a -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 2  -m 350 -s 20 -o normal2 -rs $RANDOM &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 5b -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 5  -m 350 -s 20 -o normal5.1 -rs $RANDOM &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 5c -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 5  -m 350 -s 20 -o normal5.2 -rs $RANDOM &
        pids[${j}]=$!;
        j=$(($j+1));
	cd ../


       mkdir set0
       cd set0
        art_illumina -d 0a -ss HS20 -i ../simulated.fa -p -l 100 -f 1.8  -m 350 -s 20 -o tumor1.8 -rs $RANDOM && bwa-0.7.15 mem -t 2 $reference tumor1.81.fq tumor1.82.fq | samtools view -bS > tumor1.8.bam && samtools sort --threads 2 tumor1.8.bam > tumor1.8_sorted.bam && samtools index tumor1.8_sorted.bam && rm tumor1.8.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 0b -ss HS20 -i ../simulated.fa -p -l 100 -f 0.45  -m 350 -s 20 -o tumor0.45.1 -rs $RANDOM && bwa-0.7.15 mem -t 1 $reference tumor0.45.11.fq tumor0.45.12.fq | samtools view -bS > tumor0.45.1.bam && samtools sort --threads 1 tumor0.45.1.bam > tumor0.45.1_sorted.bam && samtools index tumor0.45.1_sorted.bam && rm tumor0.45.1.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 0c -ss HS20 -i ../simulated.fa -p -l 100 -f 0.45  -m 350 -s 20 -o tumor0.45.2 -rs $RANDOM && bwa-0.7.15 mem -t 1 $reference tumor0.45.21.fq tumor0.45.22.fq | samtools view -bS > tumor0.45.2.bam && samtools sort --threads 1 tumor0.45.2.bam > tumor0.45.2_sorted.bam && samtools index tumor0.45.2_sorted.bam  && rm tumor0.45.2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 0d -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 0.45  -m 350 -s 20 -o normal0.45.1 -rs $RANDOM && bwa-0.7.15 mem -t 1 $reference normal0.45.11.fq normal0.45.12.fq | samtools view -bS > normal0.45.1.bam && samtools sort --threads 1 normal0.45.1.bam > normal0.45.1_sorted.bam && samtools index normal0.45.1_sorted.bam &&  rm normal0.45.1.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d 0e -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 0.45  -m 350 -s 20 -o normal0.45.2 -rs $RANDOM && bwa-0.7.15 mem -t 1 $reference normal0.45.21.fq normal0.45.22.fq | samtools view -bS > normal0.45.2.bam  && samtools sort --threads 1 normal0.45.2.bam > normal0.45.2_sorted.bam && samtools index normal0.45.2_sorted.bam &&  rm normal0.45.2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
       art_illumina -d 0f -ss HS20 -i ../germ_simulated.fa -p -l 100 -f 0.3  -m 350 -s 20 -o normal0.3 -rs $RANDOM && bwa-0.7.15 mem -t 1 $reference normal0.31.fq normal0.32.fq | samtools view -bS > normal0.3.bam && samtools sort --threads 1 normal0.3.bam > normal0.3_sorted.bam && samtools index normal0.3_sorted.bam && rm normal0.3.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        cd ../

for i in `seq 1 4`; do
	mkdir set$i
	cd set$i
	art_illumina -d $i\a -ss HS20 -i ../simulated.fa -p -l 100 -f 3  -m 350 -s 20 -o tumor3 -rs $RANDOM &&  bwa-0.7.15 mem -t 4 $reference tumor31.fq tumor32.fq | samtools view -bS > tumor3.bam  && samtools sort --threads 4 tumor3.bam > tumor3_sorted.bam && samtools index tumor3_sorted.bam && rm tumor3.bam &
        pids[${j}]=$!;
	j=$(($j+1));
        art_illumina -d $i\b -ss HS20 -i ../simulated.fa -p -l 100 -f 0.75  -m 350 -s 20 -o tumor0.75.1 -rs $RANDOM &&  bwa-0.7.15 mem -t 1 $reference tumor0.75.11.fq tumor0.75.12.fq | samtools view -bS > tumor0.75.1.bam  && samtools sort --threads 1 tumor0.75.1.bam > tumor0.75.1_sorted.bam && samtools index tumor0.75.1_sorted.bam && rm tumor0.75.1.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d $i\c -ss HS20 -i ../simulated.fa -p -l 100 -f 0.75  -m 350 -s 20 -o tumor0.75.2 -rs $RANDOM &&  bwa-0.7.15 mem -t 1 $reference tumor0.75.21.fq tumor0.75.22.fq | samtools view -bS > tumor0.75.2.bam && samtools sort --threads 1 tumor0.75.2.bam > tumor0.75.2_sorted.bam && samtools index tumor0.75.2_sorted.bam && rm tumor0.75.2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d $i\d -ss HS20 -i ../germ_simulated.fa -p -l 100  -f 0.75  -m 350 -s 20 -o normal0.75.1 -rs $RANDOM &&  bwa-0.7.15 mem -t 1 $reference normal0.75.11.fq normal0.75.12.fq | samtools view -bS > normal0.75.1.bam &&  samtools sort --threads 1 normal0.75.1.bam > normal0.75.1_sorted.bam && samtools index normal0.75.1_sorted.bam && rm normal0.75.1.bam &
        pids[${j}]=$!;
        j=$(($j+1));
        art_illumina -d $i\e -ss HS20 -i ../germ_simulated.fa -p -l 100 -f 0.75  -m 350 -s 20 -o normal0.75.2 -rs $RANDOM &&  bwa-0.7.15 mem -t 1 $reference normal0.75.21.fq normal0.75.22.fq | samtools view -bS > normal0.75.2.bam && samtools sort --threads 1 normal0.75.2.bam > normal0.75.2_sorted.bam && samtools index normal0.75.2_sorted.bam && rm normal0.75.2.bam &
        pids[${j}]=$!;
        j=$(($j+1));
       art_illumina -d $i\f -ss HS20 -i ../germ_simulated.fa -p -l 100 -f 0.5  -m 350 -s 20 -o normal0.5 -rs $RANDOM &&  bwa-0.7.15 mem -t 1 $reference normal0.51.fq normal0.52.fq | samtools view -bS > normal0.5.bam &&  samtools sort --threads 1 normal0.5.bam > normal0.5_sorted.bam && samtools index normal0.5_sorted.bam && rm normal0.5.bam &
        pids[${j}]=$!;
        j=$(($j+1));
	cd ../

done

for pid in ${pids[*]}; do
        wait $pid;
done


cd set5
bwa-0.7.15 mem -t 40 $reference normal21.fq normal22.fq | samtools view -bS > normal2.bam && samtools sort --threads 40 normal2.bam > normal2_sorted.bam && samtools index normal2_sorted.bam && rm normal2.bam 
bwa-0.7.15 mem -t 40 $reference normal5.11.fq normal5.12.fq | samtools view -bS > normal5.1.bam && samtools sort --threads 40 normal5.1.bam > normal5.1_sorted.bam && samtools index normal5.1_sorted.bam && rm normal5.1.bam
bwa-0.7.15 mem -t 40 $reference normal5.21.fq normal5.22.fq | samtools view -bS > normal5.2.bam && samtools sort --threads 40 normal5.2.bam > normal5.2_sorted.bam && samtools index normal5.2_sorted.bam && rm normal5.2.bam
cd ../





rm set0/*.aln
rm set1/*.aln
rm set2/*.aln
rm set3/*.aln
rm set4/*.aln
rm set5/*.aln
