[![Gitpod ready-to-code](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/renatopuga/somatico)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1r4LDUiQqirUFQT2nIrhVW_AHLaBhmS7y?usp=sharing)


# Análise Somática
GATK 4 Mutect2 Somático

# T12020 - Aula Prática

![image](https://user-images.githubusercontent.com/8321336/130251648-7ad77cae-435f-49be-950f-b7af5f59fd7a.png)

## Pipeline

* GATK4 - Mutect2
* Gene JAK2
* Referência chr9
  * Sobre as versões do Genoma Humano: https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies#grch37
* Amostras: 
  * WP043 (tumor)
  * WP044 (normal)
* https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-


## Amostras Extras

- WP190 (tumor) e WP191 (normal)
- WP017 (tumor) e WP018 (normal)

**Nota 1:** Utilizar o af-gnomad chr13 e chr19. 

**Nota 2:** Será preciso baixar o chr13 e chr19 da UCSC.

**Download chr19**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
```

**Download chr13**
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
```



**Concatenar os arquivos .fa.gz**
> Dica 1: zcat lê arquivos compactados .gz e zip

```bash
zcat chr13.fa.gz chr19.fa.gz | sed -e "s/chr//g" > hg19.fa
```

**Gerar o index do arquivo hg19.fa**
```bash
samtools faidx hg19.fa
```

 ____________________ 
< Primeiro é o chr9 >
 -------------------- 
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||



## Workflow
Os arquivos BAM (tumor e normal) já foram gerados e estão prontos para a chamada de variates (ver Anexo 1). Então, agora vamos fazer download da referência `chr9` e gerar o index com o programa `samtools`.

### Download da Referência - chr9

* Download 

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
```

* Alterar nome do header: DE: >chr9 para >9
> Essa alteração é necessária pois no BAM a referência não tinha `>chr` era apenas `>9`.

```bash
zcat chr9.fa.gz | sed -e "s/chr//g" > chr9.fa
```

* Verificar se alteração foi feita com o comando `head`

```bash
head chr9.fa
```
> O comando `head` lê as 10 primeiras linha de um arquivo texto


## samtools

Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:



### samtools install. Escolher um dos instaladores (todo terminal virtual é Ubuntu)

* samtools install (Mac)

```bash
brew install samtools 
```

* samtools install (Ubuntu)

```bash
sudo apt-get install samtools 
```

* samtools install (Docker)

```bash
docker pull biocontainers/samtools
```



### samtools faidx e index (cria um dicionário do arquivo chr9.fa)

* samtools faidx

```bash
samtools faidx chr9.fa
```


## GATK4

> Version: 4.2.2.0

Genome Analysis Toolkit - Variant Discovery in High-Throughput Sequencing Data. https://gatk.broadinstitute.org/



### GATK4 install

GATK4 install Docker

```bash
docker pull broadinstitute/gatk:4.2.2.0
```



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

```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary -R chr9.fa -O chr9.dict
```



### GATK4 intervals

```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs -R chr9.fa -O chr9.interval_list -OT ACGT ##ACTG é a descrição de quais bases devem ser selecionadas no novo arquivo, o que exclui os Ns
```



## Mutect2

Call somatic SNVs and indels via local assembly of haplotypes



### Mutect2 Tumor e Normal
> O comando: `samtools view -H normal_JAK2.bam` você consegue pegar o campo SM: que contém o ID da amostra normal.

```bash
./gatk-4.2.2.0/gatk Mutect2 \
	-R chr9.fa \ ##referencia a ser usada
	-I tumor_JAK2.bam \ ##input 
	-I normal_JAK2.bam \  ##outro input
	-normal WP044 \  ## argumento "-normal" direciona qual a amostra normal
	--germline-resource af-only-gnomad-chr9.vcf.gz \  ##frequencia alelica vinda do Gnomad para saber se é germinativo ou não
	-O somatic.vcf.gz \
	-L chr9.interval_list  ##quais regiões devem ser consideradas para chamar variante
```
##extra 

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

Tabulates pileup metrics for inferring contamination

* GetPileupSummaries Tumor

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I tumor_JAK2.bam \
	-V af-only-gnomad-chr9.vcf.gz \  (referência da frequência alélica)
	-L chr9.interval_list \
	-O tumor_JAK2.table
```

* GetPileupSummaries Normal

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I normal_JAK2.bam \  
	-V af-only-gnomad-chr9.vcf.gz \
	-L chr9.interval_list \
	-O normal_JAK2.table
```



### CalculateContamination

Calculate the fraction of reads coming from cross-sample contamination (compara quantitativamente a frequência alélica de todas as variantes encontradas com a freq. alélica descrita pelo gnomad, neste caso, gerando score e desvio padrão. Para o CalculateContamination o score acima de 2 é contaminação. Importante lembrar que estamos olhando só um gene e esse cálculo pode ser pouco eficiente pela região ser pequena)

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_JAK2.table \ ##arquivo de interesse 
	-matched normal_JAK2.table \  #com qual será comparado
	-O contamination.table
```



### FilterMutectCalls

Filter somatic SNVs and indels called by Mutect2 (considera diversos fatores como qualidade da base, da read, do alinhamento, das regiões ao redor, se está presente nas duas fitas, se está presente só na amostra normal ou no tumor para determinar se é artefato "artefact" ou variante mesmo "PASS") 

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R chr9.fa \
	-V somatic.vcf.gz \  ##arquivo vcf a ser usado
	--contamination-table contamination.table \  ##table de contaminação 
	-O filtered.vcf.gz
```
##* é interessante abrir o VCF, bam e bai das amostras no IGV para visualizar as regiões filtradas


## VEP ensembl - Anotação

### VEP install

```bash
docker pull ensemblorg/ensembl-vep
```

* Criar diretório de output do VEP com permissão total (aplicado apenas no gitpod)

Voltar para a casa
```
cd
```

Verificar se o diretório é o /home/gitpod
```
pwd
```
> ex.: /home/gitpod

Copiar o arquivo filtered.vcf.gz
```
cp /workspace/somatico/filtered.vcf.gz .
```

Copiar o arquivo chr9.fa
```
cp /workspace/somatico/chr9.fa .
```

Copiar o arquivo chr9.fa.fai
```
cp /workspace/somatico/chr9.fa.fai .
```


Criar o diretorio vep_output
```
mkdir -p vep_output
```

Modificar a permissao do diretorio vep_output
```
chmod 777 vep_output
```

* Aplicar apenas no Google Colab

```bash
mkdir -p vep_output
chmod 777 vep_output
```

# rodar o vep

```bash
docker run -it --rm  -v $(pwd):/data ensemblorg/ensembl-vep ./vep  \
	-i /data/filtered.vcf.gz  \
	-o /data/vep_output/filtered.vep.tsv \
	--database --assembly GRCh37 --refseq  \
	--pick --pick_allele --force_overwrite --tab --symbol --check_existing \
	--fields "Location,SYMBOL,Consequence,Feature,Amino_acids,CLIN_SIG" \
	--fasta /data/chr9.fa \
	--individual all 
```



## Panel of Normal (PoN)


# Panel of Normal (PoN)

GATK Best Practices - Exome PoN


* vcf

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf
```

* vcf.idx

```bash
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx
```

* Mutect2

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  --panel-of-normals Mutect2-exome-panel.vcf \
  -L hg19.interval_list \
  -O WP190.somatic.pon.vcf.gz
  
```

* CalculateContamination somente com o table do tumor (ex.: wp190)

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
	-I tumor_wp190.table \
	-O WP190.contamination.pon.table
```

* FilterMutectCalls

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R hg19.fa \
	-V WP190.somatic.pon.vcf.gz \
	--contamination-table WP190.contamination.pon.table \
	-O WP190.filtered.pon.vcf.gz
```





# Anexo 1

### Converter BAM para FASTQ

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP043/WP043.bam 9:5021937-5126899 | samtools bam2fq -1 tumor_R1.fq -2 tumor_R2.fq - 
```

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP044/WP044.bam 9:5021937-5126899 | samtools bam2fq -1 normal_R1.fq -2 normal_R2.fq - 
```



### Converter BAM para JAK2 BAM

Aqui estão os comandos que foram utilizados para gerar os BAMs intermediários e como gerar arquivos FASTQs de regiões específicas do seu arquivo BAM completo.
> Essas etapas não precisam ser executadas nesse pipeline

* BAM para BAM

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP043/WP043.bam 9:5021937-5126899 > tumor_JAK2.bam
```

```bash
samtools view -h -b /Volumes/Seagate\ Expansion\ Drive/data-lpfap10/projects/proadi/exome/bam/WP044/WP044.bam 9:5021937-5126899 > normal_JAK2.bam
```

* Gerar index do BAM (.BAI)

```bash
samtools index tumor_JAK2.bam 
```

```bash
samtools index normal_JAK2.bam 
```



### af-only-gnomad.vcf.gz (apenas região JAK2)
> [Google Clou af-only-gnomad.vcf.gz](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37;tab=objects?project=broad-dsde-outreach&pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&forceOnObjectsSortingFiltering=false&pli=1)

* Header do VCF

```bash
zgrep -w "\#" af-only-gnomad.raw.sites.chr.vcf.gz > header
```

* Apenas Região do Gene JAK2

```bash
zgrep -w "^chr9" af-only-gnomad.raw.sites.chr.vcf.gz  | awk '$2>=5021937 && $2<=5126899' > JAK2.region
```

* Concatenar header + vcf

```bash
cat header JAK2.region > af-only-gnomad-chr9.vcf
```

* Compactar

```bash
bgzip af-only-gnomad-chr9.vcf
```

* Index do VCF

```bash
tabix -p vcf af-only-gnomad-chr9.vcf.gz 
```

# Anexo 2

Rodar amostras extras pareadas:

* Com par:
```bash
sh run.chr13-chr19.sh tumor_JAK2.bam normal_JAK2.bam WP043 WP044
sh run.chr13-chr19.sh tumor_wp017.bam tumor_wp018.bam WP017 WP018
sh run.chr13-chr19.sh tumor_wp190.bam tumor_wp191.bam WP190 WP191
```

* Com PoN:
```bash
sh run.chr13-chr19.pon.sh tumor_JAK2.bam normal_JAK2.bam WP043 WP044 
sh run.chr13-chr19.pon.sh tumor_wp017.bam tumor_wp018.bam WP017 WP018 
sh run.chr13-chr19.pon.sh tumor_wp190.bam tumor_wp191.bam WP190 WP191 

```


**Código Fonte (run.chr13-chr19.pon.sh e run.chr13-chr19.sh)**

```bash

# variaveis fixas
gnomad="af-only-gnomad-chr13-chr19.vcf.gz"
interval="hg19.interval_list"
genome="hg19.fa"

# bam e ids 
tumor=$1
normal=$2
id_tumor=$3
id_normal=$4

# rodando mutect2
./gatk-4.2.2.0/gatk Mutect2 \
	-R $genome \
	-I $tumor \
	-I $normal \
	-normal $id_normal \
	--germline-resource $gnomad \
	-O $id_tumor.somatic.vcf.gz \
	-L $interval

# getPileup tumor
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $tumor \
	-V $gnomad \
	-L $interval \
	-O $id_tumor.table

# getPileup normal
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $normal \
	-V $gnomad \
	-L $interval \
	-O $id_normal.table

# CalculateContamination
./gatk-4.2.2.0/gatk CalculateContamination \
	-I $id_tumor.table \
	-matched $id_normal.table \
	-O $id_tumor.contamination.table

# FilterMutectCalls
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R $genome \
	-V $id_tumor.somatic.vcf.gz \
	--contamination-table $id_tumor.contamination.table \
	-O $id_tumor.filtered.vcf.gz

```

Rodar amostras extras com Panel of Normal (broad institute).

```bash

# variaveis fixas
gnomad="af-only-gnomad-chr13-chr19.vcf.gz"
interval="hg19.interval_list"
genome="hg19.fa"
pon="Mutect2-exome-panel.vcf"

# bam e ids 
tumor=$1
id_tumor=$3

# rodando mutect2
./gatk-4.2.2.0/gatk Mutect2 \
	-R $genome \
	-I $tumor \
	--germline-resource $gnomad \
	--panel-of-normals $pon \
	-O $id_tumor.somatic.pon.vcf.gz \
	-L $interval

# getPileup tumor
./gatk-4.2.2.0/gatk GetPileupSummaries \
	-I $tumor \
	-V $gnomad \
	-L $interval \
	-O $id_tumor.pon.table

# CalculateContamination
./gatk-4.2.2.0/gatk CalculateContamination \
	-I $id_tumor.pon.table \
	-O $id_tumor.contamination.pon.table

# FilterMutectCalls
./gatk-4.2.2.0/gatk FilterMutectCalls \
	-R $genome \
	-V $id_tumor.somatic.pon.vcf.gz \
	--contamination-table $id_tumor.contamination.pon.table \
	-O $id_tumor.filtered.pon.vcf.gz
```
