## Análise Somática com GATK 4 Mutect2 Somático dos cromossomos 13 e 19 para as amostras: ##

WP190 (tumor) e WP191 (normal)
WP017 (tumor) e WP018 (normal)


**Download chr19**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
```

**Download chr13**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
```

**Concatenar os arquivos .fa.gz**

```bash
zcat chr13.fa.gz chr19.fa.gz | sed -e "s/chr//g" > hg19.fa
```

**Gerar o index do arquivo hg19.fa**
```bash
samtools faidx hg19.fa
```

**Alterar nome do header: DE: >chr19 para >19 e >chr13 para >13**
> Essa alteração é necessária pois no BAM a referência não tinha `>chr` era apenas `>19`.

```bash
zcat hg19.fa.gz | sed -e "s/chr//g" > hg19.fa
```

**Verificar se alteração foi feita com o comando `head`**

```bash
head hg19.fa
```

**samtools install**

```bash
brew install samtools 
```

**samtools faidx e index**
>cria um dicionário do arquivo .fa

```bash
samtools faidx chr9.fa
```

## GATK4

> Version: 4.2.2.0

Genome Analysis Toolkit - Variant Discovery in High-Throughput Sequencing Data. https://gatk.broadinstitute.org/



### GATK4 install

GATK4 install Mac e Linux

* Download

```bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
```

* Descompactar

```bash
unzip gatk-4.2.2.0.zip 
```

* Testando gatk

```bash
./gatk-4.2.2.0/gatk
```


### GATK4 .dict
# Cria o dicionário do arquivo com informações gerais, por exemplo, nucleotide length (LN)

```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
```



### GATK4 intervals
#Dispõe intervalos no arquivo, excluindo as bases substituídas por N

```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R hg19.fa -O hg19.interval_list -OT ACGT 
```


## Mutect2
# Chamada de variantes por comparação de haplótipos 

### Mutect2 Tumor e Normal
#Uso da referência hg19 para comparação, sinalizando que a amostra normal é a WP018 ou WP191

```bash
./gatk-4.2.2.0/gatk Mutect2 \
	-R hg19.fa \                                             ##referencia a ser usada
	-I tumor_wp017.bam \                                     ##input tumoral
	-I normal_wp018.bam \                                    ##input normal
	-normal WP018 \                                          ##argumento para identificar a amostra normal
	--germline-resource af-only-gnomad-chr13-chr19.vcf.gz \  ##frequencia alelica do Gnomad que será usada para determinar se a variante é germinatia
	-O somatic_wp017.vcf.gz \                                ##arquivo gerado
	-L hg19.interval_list                                    ##quais regiões devem ser consideradas para a chamada de variantes
```

```bash
./gatk-4.2.2.0/gatk Mutect2 \
	-R hg19.fa \
	-I tumor_wp190.bam \ 
	-I tumor_wp191.bam \ 
	-normal WP191 \  
	--germline-resource af-only-gnomad-chr13-chr19.vcf.gz \  
	-O somatic_wp190.vcf.gz \
	-L hg19.interval_list  
```

## Calcular Contaminação



### GetPileupSummaries

* GetPileupSummaries para tecido tumoral

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_wp017.bam \                      #input 
	-V af-only-gnomad-chr13-chr19.vcf.gz \    #vcf 
	-L hg19.interval_list \                   #
	-O tumor_wp017.table                      #arquivo gerado
```

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_wp190.bam \
	-V af-only-gnomad-chr13-chr19.vcf.gz \  (referência da frequência alélica)
	-L hg19.interval_list \
	-O tumor_wp190.table
```

* GetPileupSummaries para tecido normal

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_wp018.bam \
	-V af-only-gnomad-chr13-chr19.vcf.gz \  
	-L hg19.interval_list \
	-O tumor_wp018.table
```

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_wp191.bam \
	-V af-only-gnomad-chr13-chr19.vcf.gz \  
	-L hg19.interval_list \
	-O tumor_wp191.table
```


### CalculateContamination

Compara quantitativamente a frequência alélica de todas as variantes encontradas com a frequência alélica descrita pelo Gnomad, neste caso, gerando score e desvio padrão. 
Para o CalculateContamination o score acima de 2 é contaminação. 
Importante lembrar que estamos olhando só um gene e esse cálculo pode ser pouco eficiente pela região ser pequena

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_wp019.table \                       #amostra de interesse 
	-matched tumor_wp018.table \                 #amostra normal para comparação 
	-O contaminationwp17.table
```

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_wp190.table \                       
	-matched tumor_wp191.table \                  
	-O contaminationwp190.table
```



### FilterMutectCalls

Considera diversos fatores como qualidade da base, da read, do alinhamento, das regiões ao redor, se está presente nas duas fitas, se está presente só na amostra normal ou no tumor,
para determinar se é artefato ("artefact") ou variante ("PASS") 

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R hg19.fa \                                 #referência
	-V somaticwp017.vcf.gz \                     #arquivo vcf a ser usado
	--contamination-table contamination.table \  #tabela de contaminação 
	-O filtered.vcf.gz                           #arquivo gerado
```
##* é interessante abrir o VCF, bam e bai das amostras no IGV para visualizar as regiões filtradas


