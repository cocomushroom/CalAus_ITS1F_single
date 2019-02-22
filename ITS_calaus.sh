### ITS workflow for ITS1F/ITS4 sequenced data, use Australia/California MiSeq data 301PE here ##https://benjjneb.github.io/dada2/ITS_workflow.html

## cutadapt ##

for file in `find . -name "*R1*fastq.gz"`; do fastaFile=${file}; cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATATCAATAAGCGGAGGA -e 0.2 -o cut_R1/${fastaFile}_trimmed.fastq.gz ${fastaFile} >>log_R1; done


for file in `find . -name "*R1*fastq.gz"`; do fastaFile=${file}; cutadapt -g GTGARTCATCGAATCTTTG -a GCATATCAATAAGCGGAGGA -e 0.2 -o cut_R1/${fastaFile}_trimmed.fastq.gz ${fastaFile} >>log_R1; done



for i in .; do /Users/liaolab/bin/FastQC.app/Contents/MacOS/fastqc *; done

multiqc .   



for file in `find . -name "*R1*fastq.gz"`; do fastaFile=${file}; java -jar /Users/liaolab/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 ${fastaFile} trimmo/${fastaFile}_trimo.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36>>trimmolog_R1; done

for file in `find . -name "*R1*fastq.gz"`; do fastaFile=${file}; wc -l ${fastaFile} >> line_cutadapt_trimmo; done



for f in `find . -name "*R1*fastq.gz"`; do mv "$f" $(echo "$f" | sed 's/\.fastq\.gz_trimmed\.fastq\.gz//g'); done

for f in `find . -name "*R1*fastq.gz"`; do mv "$f" $(echo "$f" | sed 's/_B_/_/g'); done

## test rename sometime



# For some reason still find high percentage of illumina adaptor (>70%). Use Trimmomatic as a 2nd clean up##




