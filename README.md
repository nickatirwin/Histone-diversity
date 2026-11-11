# Histone-diversity
Associated files and methodological details for the analysis presented in:

Irwin, N. A. T, Casado, D. K., Cho, C. H. 2026. The origin and diversification of the nucleosome.
## Methods
We aimed to examine the diversity of histone proteins across diverse eukaryotic species. To do this, we identified all histones representing H2A, H2A.Z, H2B, H3, CenH3, and H4 using hidden Markov models, curated them phylogenetically, and compared them to their respective canonical human histone ortholog (H2A: H2AC4, H2A.Z: H2AZ2, H2B: H2BC21, H3: H3-3A, CenH3: CENPA, H4: H4C1) using BLAST.

#### 0. Requirements
```
conda install conda-forge::ete3
conda install bioconda::diamond
conda install bioconda::mafft
conda install bioconda::trimal
conda install bioconda::iqtree
```

#### 1. Collect genome-predicted proteomes for representative taxa
| Species | Taxa ID | Source | Accession |
| ------------- | ------------- | ------------- | ------------- |
| Ectocarpus siliculosus | 2880 | NCBI | GCA_000310025.1
| Oxytricha trifallax | 1172189 | NCBI | GCA_000295675.1
| Symbiodinium microadriaticum | 2951 | SRA | ERR3711255
| Micromonas pusilla CCMP1545 | 564608 | NCBI | GCF_000151265.4
| Cyanidioschyzon merolae strain 10D | 280699 | NCBI | GCF_000091205.1
| Phaeodactylum tricornutum CCAP 1055/1 | 556484 | NCBI | GCF_000150955.2
| Capsaspora owczarzaki ATCC 30864 | 595528 | NCBI | GCF_000151315.2
| Ostreococcus tauri | 70448 | NCBI | GCF_000214015.3
| Drosophila melanogaster | 7227 | NCBI | GCF_000001215.4
| Phytophthora sojae | 67593 | EukProt | EP00647
| Ceratopteris richardii | 49495 | NCBI | GCA_020310875.1
| Arabidopsis thaliana | 3702 | NCBI | GCF_000001735.4
| Bigelowiella natans | 227086 | EukProt | EP00466
| Chlamydomonas reinhardtii | 3055 | NCBI | GCF_000002595.2
| Salpingoeca rosetta | 946362 | NCBI | GCF_000188695.1
| Planoprotostelium fungivorum | 1890364 | NCBI | GCA_003024175.1
| Giardia duodenalis | 5741 | EukProt | EP00701
| Toxoplasma gondii VEG | 432359 | NCBI | GCA_000150015.2
| Plasmodium falciparum 3D7 | 36329 | NCBI | GCF_000002765.6
| Acanthamoeba castellanii str. Neff | 1257118 | NCBI | GCF_000313135.1
| Saccharomyces cerevisiae S288C | 559292 | NCBI | GCF_000146045.2
| Schizosaccharomyces pombe | 4896 | NCBI | GCF_000002945.2
| Nematostella vectensis | 45351 | NCBI | GCF_932526225.1
| Homo sapiens | 9606 | NCBI | GCF_000001405.40
| Chrysochromulina tobinii | 1460289 | NCBI | GCA_001275005.1
| Klebsormidium nitens | 105231 | NCBI | GCA_000708835.1
| Spironucleus salmonicida | 348837 | EukProt | EP00702
| Trichomonas vaginalis | 5722 | EukProt | EP00705
| Naegleria gruberi | 5762 | JGI | Naegr1
| Marchantia polymorpha subsp. ruderalis | 1480154 | NCBI | GCA_037833965.1
| Amphimedon queenslandica | 400682 | NCBI | GCF_000090795.2
| Tetrahymena thermophila SB210 | 312017 | NCBI | GCF_000189635.1
| Danio rerio | 7955 | NCBI | GCF_000002035.6
| Caenorhabditis elegans | 6239 | NCBI | GCF_000002985.6
| Trypanosoma brucei brucei TREU927 | 185431 | NCBI | GCF_000002445.2
| Amphidinium carterae | 2961 | SRA | SRR31183665
| Mucor circinelloides 1006PhL | 1220926 | NCBI | GCA_000401635.1
| Volvox carteri f. nagariensis | 3068 | NCBI | GCF_000143455.1
| Dictyostelium discoideum AX4 | 352472 | NCBI | GCF_000004695.1
| Emiliania huxleyi CCMP1516 | 280463 | NCBI | GCF_000372725.1
| Chondrus crispus | 2769 | NCBI | GCF_000350225.1
| Guillardia theta CCMP2712 | 905079 | NCBI | GCF_000315625.1
| Physcomitrium patens | 3218 | NCBI | GCF_000002425.4
| Plasmodiophora brassicae | 37360 | EukProt | EP00473
| Cryptococcus bacillisporus CA1280 | 1296109 | NCBI | GCF_000836335.1
| Leishmania major strain Friedlin | 347515 | NCBI | GCF_000002725.2
| Cafeteria roenbergensis | 33653 | NCBI | GCA_008330645.1

#### 2. Identify histones using hidden Markov models and phylogenetically curate them
Identify histone homologs using HMMs generated in Irwin & Richards (2024)
```
cat proteomes/*fasta > eukaryotes.fasta

for i in `ls *hmm | cut -f 1 -d '.'`; do hmmsearch -E 1e-5 --incE 1e-5 --domE 1e-5 --cpu 2 --tblout $i.hmmhits $i.hmm eukaryotes.fasta; done

for i in *hmmhits; do faSomeRecords eukaryotes.fasta $i $i.fasta; done
```
Phylogenetically annotate and extract histones from each family.
Note: All python scripts are available from: https://github.com/nickatirwin/phylogenomic_scripts
```
# make phylogenies and curate into families
for i in *hmmhits.fasta; do mafft --thread 32 $i > $i.aln; done
for i in *aln; do trimal -in $i -out $i.trimal -gt 0.5; done
ls *trimal | parallel -j 4 "iqtree -m LG4M+R5 -s {} -fast --runs 10 --alrt 1000"
# annotate
for i in *hmmhits.fasta; do python ncbi_taxid_to_annotation.py $i; done
for i in *hmmhits.fasta; do diamond blastp --query $i --db uniprot_sprot.dmnd --outfmt 6 --threads 32 --sensitive --out $i.blastp --evalue 1e-5; done
for i in *hmmhits.fasta; do python blast_annotation.py $i $i.blastp; done

# extract families
python remove_and_extract.py H2A.hmmhits.fasta H2A.hmmhits.fasta.aln.trimal.treefile.coloured NA 0000ff
mv H2A.hmmhits.fasta.extract H2A.all.fasta
python remove_and_extract.py H2A.hmmhits.fasta H2A.hmmhits.fasta.aln.trimal.treefile.coloured NA ff0000
mv H2A.hmmhits.fasta.extract H2AZ.all.fasta
python remove_and_extract.py H2B.hmmhits.fasta H2B.hmmhits.fasta.aln.trimal.treefile.coloured NA 0000ff
mv H2B.hmmhits.fasta.extract H2B.all.fasta
python remove_and_extract.py H3.hmmhits.fasta H3.hmmhits.fasta.aln.trimal.treefile.coloured NA 0000ff
mv H3.hmmhits.fasta.extract H3.all.fasta
python remove_and_extract.py H3.hmmhits.fasta H3.hmmhits.fasta.aln.trimal.treefile.coloured NA ff0000
mv H3.hmmhits.fasta.extract CenH3.all.fasta
python remove_and_extract.py H4.hmmhits.fasta H4.hmmhits.fasta.aln.trimal.treefile.coloured NA 0000ff
mv H4.hmmhits.fasta.extract H4.all.fasta
```
#### 3. Create reference histone sequences for comparison
```
# H2A
>9606.P04908
MSGRGKQGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYSERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAIRNDEELNKLLGRVTIAQGGVLPNIQAVLLPKKTESHHKAKGK
# H2AZ
>9606.NP_036544.1
MAGGKAGKDSGKAKAKAVSRSQRAGLQFPVGRIHRHLKTRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIGKKGQQKTA
# H2B
>9606.NP_003519.1
MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSIYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK
# H3
>9606.NP_001365972.1
MARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA
# CenH3
>9606.P49450
MGPRRRSRKPEAPRRRSPSPTPTPGPSRRGPSLGASSHQHSRRRQGWLKEIRKLQKSTHLLIRKLPFSRLAREICVKFTRGVDFNWQAQALLALQEAAEAFLVHLFEDAYLLTLHAGRVTLFPKDVQLARRIRGLEEGLG
# H4
>9606.NP_001029249.1
MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG
```
#### 4. Calculate similarities between all histones and the human orthologs
```
for i in *human.fasta; do diamond makedb --in $i --db $i.dmnd; done
for i in `ls *all.fasta | cut -f 1 -d '.'`; do diamond blastp --query $i.all.fasta --db $i.human.fasta.dmnd --outfmt 6 --max-target-seqs 100 --threads 32 --ultra-sensitive --out $i.comparison.blastp; done
```
Extract maximum similarity and minimum similarity values for each species
```
from glob import glob

# record dictionary
d = {} # species:{H2A:[min,max],H2B:[min,max],H3:[min,max],H4:[min,max]} - min = min percent identity, max = max percent identity
# record all species
fasta = open('eukaryotes.fasta','r').read().split('>')[1:]
for seq in fasta:
    d[seq.split('.')[0]] = {'H2A':[100,0],'H2AZ':[100,0],'H2B':[100,0],'CenH3':[100,0],'H3':[100,0],'H4':[100,0]}

for fname in glob('*comparison.blastp'):
    blast = open(fname,'r').readlines()
    histone = fname.split('.')[0]
    for line in blast:
        query = line.split('.')[0]
        similarity = float(line.split('\t')[2])
        if d[query][histone][0] > similarity:
            d[query][histone][0] = similarity
        if d[query][histone][1] < similarity:
            d[query][histone][1] = similarity

# output dataframe
out = open('Histone.comparison.tab','w')
out.write('species\thistone\tvalue\n')
for sp in d:
    for h in d[sp]:
        for v in d[sp][h]:
            out.write(sp+'\t'+h+'\t'+str(v)+'\n')
out.close()
```



