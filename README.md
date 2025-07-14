# FUSE-PhyloTree
*FUSE-PhyloTree enables to identify sequence conservation modules and associate them with functional annotations through phylogenomic analysis.*

## Phylogenetic Prediction of Functional Sequence Modules

Introducing an original approach to characterizing functional motifs. This methodology encompasses:
1. Detection of conserved sequence modules (using Partial Local Multiple Alignment)
2. Phylogenetic inference of species/genes/modules/functions evolutionary histories
3. Identification of co-appearances of modules and functions

The process considers protein sequences and their associated functional annotations as input. It then returns the presence of the conserved sequence modules and of the functional annotations in each ancestral gene along the gene phylogeny.


![Pipeline](img/fig1.png)

<details><summary><strong>Detailed Pipeline</strong></summary>

Detailed flowchart of the pipeline:

```mermaid id="pipeline-diagram"
---
config:
  look: classic
  theme: default
  themeVariables:
    fontSize: 50px

---

flowchart TB

U@{ shape: manual-input, label: "User Input"} --> A
U --> B
U --> C

A@{ shape: doc, label: "sequences.fasta" } --> m1@{ shape: rectangle, label: "modules_segm.py: segmentation_and_module_phylo()"}
m1 --> m2@{ shape: docs, label: "module_segm_dir_seq/" }
m2 --> m3@{ shape: rectangle, label: "tools.py: segmentation()"}
m5@{ shape: flag, label: "paloma-D" } o--o|call| m3
m3 --> m4@{ shape: doc, label: "seq_filename.dot" }
m4 --> m6@{ shape: rectangle, label: "tools.py: modules_fasta()"}
m6 --> m7@{ shape: docs, label: "module_seq/*.fasta" }
m7 --> m8@{ shape: rectangle, label: "tools.py: all_phylo()"}
m9@{ shape: flag, label: "PhyML" } o--o|call| m8
m8 --> m12@{ shape: docs, label: "module_seq/*.tree" }
m12 --> m10@{ shape: rectangle, label: "modules_segm.py: correct_modules_tree()"}
m11@{ shape: flag, label: "TreeFix" } o--o|call| m10
m10 --> m13@{ shape: docs, label: "module_seq/*.treefix.tree" }
m10 --> m14@{ shape: doc, label: "modules_path_modules.txt" }
m13 o--o|describe paths of| m14

B@{ shape: doc, label: "sequences.tree" } --> g1{gene tree as input ?}
g1 -->|yes| g2@{ shape: doc, label: "rooted binary gene tree .tree" }
g1 -->|no| g3@{ shape: rectangle, label: "gene_phylo.py: whole_phylo()"}
A -.-> g3
g4@{ shape: flag, label: "Muscle" } o--o|call| g3
g5@{ shape: flag, label: "TrimAl" } o--o|call| g3
g6@{ shape: flag, label: "PhyML" } o--o|call| g3
g7@{ shape: flag, label: "TreeFix" } o--o|call| g3
g8@{ shape: flag, label: "PhyML (branch len)" } o--o|call| g3
g3 --> g2

A --> s1@{ shape: rectangle, label: "species_phylo.py: taxid_from_fasta()"}
s2@{ shape: cyl, label: "NCBI Taxonomy" } o--o|use| s1
s1 --> s3@{ shape: doc, label: "species.tree" }

A --> d1@{ shape: rectangle, label: "tools.py: known_domains()"}
d1 --> d2@{ shape: doc, label: "domains.csv" }

g2 --> r1@{ shape: rectangle, label: "tools.py: seadog_md()"}
m14 --> r1
s3 --> r1
r2@{ shape: flag, label: "Seadog-MD" } o--o|call| r1
r1 --> r3@{ shape: doc, label: "seadogMD.output" }
r3 --> r4@{ shape: rectangle, label: "integrates_3phylo.py: write_sp_gene_event()"}
r4 --> r5@{ shape: doc, label: "seadogMD_sp_gene_event.csv" }
r4 --> r6@{ shape: doc, label: "seadogMD_gene.tree" }

C@{ shape: doc, label: "functional annotations.csv" } --> a1@{ shape: rectangle, label: "ances_scenario.py: acs_inference()"}
r5 --> a1
r6 --> a1
a1 --> a2@{ shape: doc, label: "pastml_seadogMD.csv" }
a2 --> a3@{ shape: rectangle, label: "tools.py: pastml()"}
a4@{ shape: flag, label: "PastML" } o--o|call| a3
a3 --> a5@{ shape: doc, label: "pastml_seadogMD_combined_ancestral_states.tab" }
a3 --> a6@{ shape: docs, label: "acs_dir_seadogMD/" }

r3 --> i1@{ shape: rectangle, label: "integrate_3phylo.py: __main__()"}
g2 --> i1
a5 --> i1
d2 --> i1

i1 --> i2@{ shape: doc, label: "0_gene.tree" }
i1 --> i3@{ shape: doc, label: "1_modules_and_functions_evolution.csv" }
i1 --> i4@{ shape: doc, label: "2_module_descriptions.csv" }
i1 --> i5@{ shape: docs, label: "3_visuReconc/" }
i1 --> i6@{ shape: docs, label: "working_dir/" }

linkStyle 0 stroke:#670067, stroke-width:25px
linkStyle 1 stroke:#670067, stroke-width:25px
linkStyle 2 stroke:#670067, stroke-width:25px
linkStyle 3 stroke:#EE6677, stroke-width:20px
linkStyle 4 stroke:#EE6677, stroke-width:20px
linkStyle 5 stroke:#EE6677, stroke-width:20px
linkStyle 6 stroke:#EE6677, stroke-width:20px
linkStyle 7 stroke:#EE6677, stroke-width:20px
linkStyle 8 stroke:#EE6677, stroke-width:20px
linkStyle 9 stroke:#EE6677, stroke-width:20px
linkStyle 10 stroke:#EE6677, stroke-width:20px
linkStyle 11 stroke:#EE6677, stroke-width:20px
linkStyle 12 stroke:#EE6677, stroke-width:20px
linkStyle 13 stroke:#EE6677, stroke-width:20px
linkStyle 14 stroke:#EE6677, stroke-width:20px
linkStyle 15 stroke:#EE6677, stroke-width:20px
linkStyle 16 stroke:#EE6677, stroke-width:20px
linkStyle 17 stroke:#EE6677, stroke-width:20px
linkStyle 18 stroke:#BBBBBB, stroke-width:20px
linkStyle 19 stroke:#BBBBBB, stroke-width:20px
linkStyle 20 stroke:#BBBBBB, stroke-width:20px
linkStyle 21 stroke:#BBBBBB, stroke-width:20px
linkStyle 22 stroke:#BBBBBB, stroke-width:20px
linkStyle 23 stroke:#BBBBBB, stroke-width:20px
linkStyle 24 stroke:#BBBBBB, stroke-width:20px
linkStyle 25 stroke:#BBBBBB, stroke-width:20px
linkStyle 26 stroke:#BBBBBB, stroke-width:20px
linkStyle 27 stroke:#BBBBBB, stroke-width:20px
linkStyle 28 stroke:#66CCEE, stroke-width:20px
linkStyle 29 stroke:#66CCEE, stroke-width:20px
linkStyle 30 stroke:#66CCEE, stroke-width:20px
linkStyle 31 stroke:brown, stroke-width:20px
linkStyle 32 stroke:brown, stroke-width:20px
linkStyle 33 stroke:gold, stroke-width:20px
linkStyle 34 stroke:gold, stroke-width:20px
linkStyle 35 stroke:gold, stroke-width:20px
linkStyle 36 stroke:gold, stroke-width:20px
linkStyle 37 stroke:gold, stroke-width:20px
linkStyle 38 stroke:#4477AA, stroke-width:20px
linkStyle 39 stroke:#4477AA, stroke-width:20px
linkStyle 40 stroke:#4477AA, stroke-width:20px
linkStyle 41 stroke:#4477AA, stroke-width:20px
linkStyle 42 stroke:#4477AA, stroke-width:20px
linkStyle 43 stroke:#4477AA, stroke-width:20px
linkStyle 44 stroke:#4477AA, stroke-width:20px
linkStyle 45 stroke:#4477AA, stroke-width:20px
linkStyle 46 stroke:#4477AA, stroke-width:20px
linkStyle 47 stroke:#4477AA, stroke-width:20px
linkStyle 48 stroke:#4477AA, stroke-width:20px
linkStyle 49 stroke:gold, stroke-width:20px
linkStyle 50 stroke:#BBBBBB, stroke-width:20px
linkStyle 51 stroke:#4477AA, stroke-width:20px
linkStyle 52 stroke:#228833, stroke-width:20px
linkStyle 53 stroke:#228833, stroke-width:20px
linkStyle 54 stroke:#228833, stroke-width:20px
linkStyle 55 stroke:#228833, stroke-width:20px
linkStyle 56 stroke:#228833, stroke-width:20px
linkStyle 57 stroke:#228833, stroke-width:20px

classDef tools fill:#E0B0FF,stroke:#000,stroke-width:4px;
class m5,m9,m11,g4,g5,g6,g7,g8,s2,r2,a4,l4 tools

classDef files fill:#fff,stroke:#000,stroke-width:2px;
class m2,m4,m7,m12,m13,m14,g2,s3,r3,r5,r6,a2,a5,a6,d2,l2 files

classDef input stroke:#670067,stroke-width:15px;
class U,A,B,C,i2,i3,i4,i5,i6,l1 input

classDef scripts fill:#9BD184,stroke:#000,stroke-width:2px;
class m1,m3,m6,m8,m10,g3,s1,d1,r1,r4,a1,a3,i1,l3 scripts

```
<a href="img/detailed-pipeline.pdf" target="_blank">Detailed pipeline as a PDF</a>


<details><summary><weak>Legend of the detailed pipeline</weak></summary>

```mermaid id="legend-diagram"
---
config:
  look: classic
  theme: default
  themeVariables:
    fontSize: 50px

---

flowchart LR

l1@{ shape: doc, label: "input/output file" } --- l2@{ shape: doc, label: "intermediate file" }
l2 --- l3@{ shape: rectangle, label: "python3_script.py: function()"}
l3 --- l4@{ shape: flag, label: "external software" }

classDef legend1 fill:#E0B0FF,stroke:#000,stroke-width:4px
class l4 legend1

classDef legend2 fill:#fff,stroke:#000,stroke-width:2px;
class l2 legend2

classDef legend3 stroke:#670067,stroke-width:15px;
class l1 legend3

classDef legend4 fill:#9BD184,stroke:#000,stroke-width:2px;
class l3 legend4

```

</details>
</details>

> ❓ **Differences from PhyloCharMod:** FUSE-PhyloTree is an actively maintained fork of [PhyloCharMod](https://github.com/OcMalde/PhyloCharMod_publ), originally developed for studying the ADAMTS-TSL protein family.

## :books: Citation

[Application Note on bioRxiv](https://doi.org/10.1101/2025.03.28.645931)

For a more detailed explanation of the methodology, refer to the following article [[Dennler et al. 2023]](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011404), or the following PhD thesis (only [available in French](https://www.theses.fr/2022REN1B079)).

# :gear: Simple Usage

To run the full pipeline, simply use: 
```
python3 fuse-phylotree.py <sequences.fasta> <annotations.csv> <gene_tree.tree>
```

>To install and run this command, you need a proper environment. We strongly recommend using the available [Docker](https://github.com/OcMalde/fuse-phylotree/tree/main?tab=readme-ov-file#-docker-image-usage) image (recommended for personal computers) or [Singularity](https://github.com/OcMalde/fuse-phylotree/tree/main?tab=readme-ov-file#-singularity-image-usage) image (recommended for high-performance computing clusters).

## Input 
1. ```<sequences.fasta>```:
   
    Each sequence header must satisfies the format: ```>SeqID_taxid```, where ```SeqID``` represents the unique sequence identifier and ```taxid``` is the NCBI species taxid (e.g., ```>NP031426.2_10090```)

    ```fasta
    >NP00001.1_9913
    MKAAPTKARREGALVV
    >XP00002.1_9606
    MKVAPTTKAKRGVLVA
    >NP00003.1_7955
    MTKARRMTGGIW
    >NP00004.3_9606
    MTKARRMTGGLLAV
    ```
   
    Refer to [this file](https://github.com/OcMalde/fuse-phylotree/blob/main/data/analyse_fibulin/run_singularity_fibulin/fibulin_59.fasta) for an example
   
    Please refrain from using special characters in the header (e.g.```, |,()`":;```). Use only ```_``` as a separator.

   >These files can be generated using orthogroups and GFF files, all of which are included in the Docker and Singularity images for nine species (*Homo sapiens, Mus musculus, Bos taurus, Gallus gallus, Xenopus tropicalis, Danio rerio, Ciona intestinalis, Drosophila melanogaster, Caenorhabditis elegans*). You only need to provide a file with a list of RefSeq of interest. For detailed instructions, please refer to [To build a sequence dataset based on orthogroups](https://github.com/OcMalde/fuse-phylotree/tree/main?tab=readme-ov-file#to-build-a-sequence-dataset-based-on-orthogroups).


2. ```<annotations.csv>```:

   This file contains the functional annotations associated with the different sequences (e.g., PPI)

   Each line should be formatted as: ```SeqID,Annotation_1|Annotation_2```. Here, ```SeqID``` is the unique sequence identifier. It is separated from the list of annotations by a comma (```,```), and individual annotations are separated by pipe symbols (```|```) (e.g., ```NP_620594.1,P00451_F8|P04275_VWF```). Annotations must have more than one character in length.

    ```csv
    NP_00001.1,Interact_with_A|Interact_with_B
    XP_00002.1,Interact_with_A|Interact_with_B
    NP_00003.1,Interact_with_A
    NP_00004.3,Interact_with_A    
    ```

   Refer to [this file](https://github.com/OcMalde/fuse-phylotree/blob/main/data/analyse_fibulin/run_singularity_fibulin/ppi_shared.csv) for an example


3. ```<gene_tree.tree>```:

   This file contains a binary rooted gene tree of the sequences in newick format.

   ```newick
    ((NP_00001.1:0.16,XP_00002.1:0.23):0.87,(NP_00003.1:0.35,NP_00004.3:0.42):0.47);
    ```

   Refer to [this file](https://github.com/OcMalde/fuse-phylotree/blob/main/data/analyse_fibulin/run_singularity_fibulin/fibulin_tree_root.tree) for an example


<details><summary><strong>:bulb: For other option</strong></summary>

>Pre-computed phylogenetic trees or/and paloma module decompositions can be use, as long as they respect the required header format. 

```
usage: fuse-phylotree.py [-h] [--output_directory OUTPUT_DIRECTORY] [--iter ITER] [--ml_thres ML_THRES] [--pres_fthr PRES_FTHR] [--gain_fthr GAIN_FTHR] [--lost_fthr LOST_FTHR]
                         [--species_tree SPECIES_TREE] [--infer_gene_tree] [--plma_file PLMA_FILE] [--user_pastml_csv USER_PASTML_CSV] [--reconc_domains] [--paloma_args PALOMA_ARGS]
                         [--phyml_args PHYML_ARGS] [--treefix_args TREEFIX_ARGS] [--raxml_args RAXML_ARGS] [--seadog_args SEADOG_ARGS] [--pastml_args PASTML_ARGS] [--itol] [--itol_api ITOL_API]
                         [--itol_project_name ITOL_PROJECT_NAME]
                         multi_fasta_file [leaf_functions_csv] [gene_tree]

positional arguments:
  multi_fasta_file      Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)
  leaf_functions_csv    csv file containing for each of our sequence, the list of his functions (ex : XP_012810820.2, P59509 | P999999)
  gene_tree             Gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, must correspond to the sequences in the multi fasta file !)

optional arguments:
  -h, --help            show this help message and exit
  --output_directory OUTPUT_DIRECTORY
                        output directory name
  --iter ITER           Number of times the whole module evolution inference will be performed, ie: module tree inference; their corrections; DGS reconciliation (default: 10)
  --ml_thres ML_THRES   Module length thresold: minimum module length, shorter modules will be systematically filtered out (default: 5)
  --pres_fthr PRES_FTHR
                        Sets the presence frequency threshold for modules. A module is considered present at a given gene only if its presence frequency is strictly greater than (>) the specified
                        threshold. Default: 0.0.
  --gain_fthr GAIN_FTHR
                        Sets the gained frequency threshold for modules. A module is considered gained at a given gene only if its gain frequency is strictly greater than (>) the specified threshold.
                        Default: 0.0.
  --lost_fthr LOST_FTHR
                        Sets the lost frequency threshold for modules. A module is considered lost at a given gene only if its lost frequency is strictly greater than (>) the specified threshold.
                        Default: 0.0.
  --species_tree SPECIES_TREE
                        Species tree to use as a support for the reconciliations (WARNING, must correspond to the taxid use in the other files !)
  --infer_gene_tree     Infer gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, user should check it and reroot it - we advise to only use it if you know what you
                        are doing !)
  --plma_file PLMA_FILE
                        Paloma-D output file (.agraph format, .dot, or .oplma format)
  --user_pastml_csv USER_PASTML_CSV
                        PastML full input file, corresponding full custom states to use for the different sequence id (.csv format); eg, header: 'id,P59509,P999999', data: 'XP_012810820.2,1,0' or
                        'NP_001278744.1,0,,' ; unknown states (empty) will be inferred based on known states; sequence id will be converted to fit the reconcilied gene tree ids;
  --reconc_domains      Do a DGS reconciliation with known modules (pfam / prosite) ; not tested
  --paloma_args PALOMA_ARGS
                        Custom arguments to pass to paloma-D (e.g, --paloma_args "--thr 5 --min-size 5")
  --phyml_args PHYML_ARGS
                        Custom arguments to pass to PhyML for module trees inference (e.g, --phyml_args "--model JTT")
  --treefix_args TREEFIX_ARGS
                        Custom arguments to pass to TreeFix for gene-modules (e.g, --treefix_args "--niter 100 -D 1 -L 1" - corresponds to options inside -E from treefix or to --niter)
  --raxml_args RAXML_ARGS
                        Custom arguments to pass to RaxML (for TreeFix) for gene-modules (e.g, --raxml_args "-m PROTGAMMAJTT" - corresponds at options inside -e from treefix)
  --seadog_args SEADOG_ARGS
                        Custom arguments to pass to SEADOG-MD (e.g, --seadog_args "--DD 5 --DL 1 --DTA 20 --GD 2 --GL 1")
  --pastml_args PASTML_ARGS
                        Custom arguments to pass to PastML (e.g, --pastml_args "--prediction_method ACCTRAN -m JTT")
  --itol                Upload directly on my itol account (need --itol_api and --itol_project_name)
  --itol_api ITOL_API   User iTOL api key for batch upload
  --itol_project_name ITOL_PROJECT_NAME
                        iTOL project name where to upload
```

>**:warning: Infering the gene tree using the pipeline option ```--infer_gene_tree```** When running the analysis without a gene tree as input (with ```--infer_gene_tree```), a default rooted tree will be generated. However, for optimal results, **it is strongly recommended to infer a properly rooted gene tree prior to analysis and use it as the input for the gene_tree option**. The gene phylogenetic tree serves as a critical template for the entire analysis, thus it is essential that a reliable and accurately rooted gene tree is prepared and utilized.

</details>

## Output
The main workflow output is the list of modules and functions that are present, gained or lost at the different ancestral genes. This output is presented as a table in the file ```1_modules_and_functions_evolution.csv``` (Example [here](https://github.com/OcMalde/fuse-phylotree/blob/main/data/analyse_fibulin/run_singularity_fibulin/1_modules_and_functions_evolution.csv)). It is strongly advised to also look at the final gene tree (with internal node names) ```0_gene.tree``` to visualise the annotated gene nodes. The file ```0_gene.tree``` corresponds to the input phylogenetic tree and assigns an ancestor name  to each internal node. The description of all the modules is available in ```2_module_descriptions.csv``` and it enables to get module segments (sequences and positions) based on module names. For an interactive visualisation of these data, various iTOL files are generated in ```3_visuReconc/``` and compressed in ```3_visuReconc.zip``` for batch upload on iTOL.

0. ```0_gene.tree```:

      This file contains the binary rooted gene tree of the sequences in newick format with internal node names.

1. ```1_modules_and_functions_evolution.csv```:

      This file contains the list of modules and functions that are present, gained or lost at the different ancestral and actual genes—i.e., at all nodes in the gene tree—and records the frequency with which each presence, gain, or loss is observed across iterations.

   | gene             | modules_present   | function_present | modules_gained | function_gained | modules_lost | function_lost |
   |------------------|-------------------|------------------|----------------|-----------------|--------------|---------------|
   | XP_0000001_9606  | B23:1.0 &#124; B27:1.0 | FBN1 &#124; MFAP2 | B24:0.8 | – | – | – |
   | XP_0000002_10090 | B25:1.0 &#124; B30:1.0 | ACAN &#124; VCAN | B31:0.9 | VCAN | B23:0.5 | ACAN |
   | XP_0000003_7227  | B10:1.0           | –                | B11:1.0       | –               | B10:0.4     | – |



3. ```2_module_descriptions.csv```:

      This file contains all module descriptions: all their segments, proteins, and positions in the protein sequences.

4. ```3_visuReconc/```:

      This directory contains various iTOL files for the interactive visualisation.

5. ```X_module_function_assoc_stats.csv```:

       This file contains different co-presence statistics for every module-function pairs. Each row lists a module–function pair, flags whether they co-emerged, and reports Jaccard overlap, Fisher exact p-value, and odds ratio comparing genes that contain the module with genes annotated for the function.
   

   | module | function | co_emergence | co-presence_jaccard_index | co-presence_FET_p_value | co-presence_FET-odds_ratio |
   |--------|----------|--------------|---------------------------|-------------------------|---------------------------|
   | B734   | MFAP3    | False        | 0.11764705882352941       | 0.18209949232049583     | 3.4482758620689653        |
   | B734   | COL8A1   | True         | 0.07692307692307693       | 0.3016405747582532      | 3.888888888888889         |
   | B734   | FBXW5    | True         | 0.07692307692307693       | 0.3016405747582532      | 3.888888888888889         |


7. *Working directory* - ```working_dir_*```:

      The intermediates outputs obtained at each of the pipeline steps are provided in a directory called ```working_dir_*``` (see [this directory](https://github.com/OcMalde/fuse-phylotree/tree/main/data/analyse_fibulin/run_singularity_fibulin) for an example)

<details><summary><strong>Detailed Output & Working Directory</strong></summary>

```
working_directory
├── 0_gene.tree										-----> file contains the binary rooted gene tree of the sequences in newick format with internal node names
├── 1_modules_and_functions_evolution.csv						-----> file contains the list of modules/functions present/gained/lost at the different ancestral and actual genes, i.e., for all nodes in the gene tree
├── 2_module_descriptions.csv								-----> file contains all module descriptions: all their segments, proteins, and positions in the protein sequences
├── 3_visuReconc									-----> All itol visualisation files
│   ├── geneReconc_seadogMD_*.tree							-----> The "final" gene tree after treefix correction / internal node labelling by seadog mDGS reconciliation / branch length computing by PhyML
│   ├── itolAnnotPresence_seadogMD_*.txt						-----> Annotation presences in leaves as heatmap
│   ├── itolBarModulesNb_seadogMD_*.txt							-----> Module number in leaves as barplot
│   ├── itolDomains_seadogMD_*.txt							-----> Domain decompositions as domain mosaics
│   ├── itolGOt_seadogMD_*.txt								-----> Annotation presence as symbols
│   ├── itolModPresence_seadogMD_*.txt							-----> Module presences in leaves as heatmap
│   ├── itolModTransfer_seadogMD_*.txt							-----> Module transfers as arrows
│   ├── itol_modules_PieGainsLost_seadogMD_*.txt					-----> Number of modules gained / lost at each gene nodes as a pie charts
│   ├── itolModulesPresent_G*_*_*_only_mod.txt						-----> For each gene node * : all module segments in actual gene present at this * (ancestral) gene as domain mosaics (brown squares)
│   ├── itolModules_seadogMD_*.txt							-----> Module decompositions as domain mosaics
│   ├── itolModulesThatChanged_G*_*_*_only_mod.txt					-----> For each gene node * : all module segments in actual gene gained at this * (ancestral) gene as domain mosaics (green squares)
│   ├── itolPopup_seadogMD_*.txt							-----> Gene nodes popup with module and annotation decriptions
│   ├── itol_ppi_PieGainsLost_seadogMD_*.txt						-----> Number of modules gained / lost at each gene nodes as a pie charts
│   └── itolSpGeneEvents_seadogMD_*.txt							-----> Gene - Species reconcilation events as symbols
├── 3_visuReconc.zip									-----> All itol visualisation files compressed for itol batch upload
└── working_dir_*
    ├── gene.fasta									-----> The input fasta file
    ├── *_species.tree									-----> Species tree extracted from ncbi taxonomy
    ├── leaf_Manual.csv									-----> The input with gene functional annotation file
    ├── gene.tree                                                    			-----> The input gene tree file (if given as input)
    │
    ├── acs_dir_seadogMD_*_gene								-----> All ancestral scenario reconstruction files from pastML
    │   ├── leaf_Manual.csv								-----> The input with gene functional annotation file
    │   ├── pastml_seadogMD_*_gene_leaf_Manual_combined_ancestral_states.tab		-----> pastML output regrouping all annotation ancestral states
    │   ├── pastml_seadogMD_*_gene_leaf_Manual.csv					-----> The input annotation formated as a {0,1,?} matrix for pastML use
    │   ├── seadogMD_*_gene_pastml							-----> pastML outputs
    │   │   ├── marginal_probabilities.character_*.model_*.tab
    │   │   ├── params.character_*.method_MPPA.model_*.tab
    │   ├── seadogMD_*_gene.tree							-----> The "final" gene tree after treefix correction / internal node labelling by seadog mDGS reconciliation / branch length computing by PhyML
    │   └── seadogMD_*_sp_gene_event.csv						-----> Gene nodes event from Species - Gene reconciliation (e.g., Gene duplication, Speciation)
    │
    ├── gene_phylo_dir_*								-----> All tree construction files
    │   ├── gene.fasta									-----> The input fasta file
    │   ├── gene.tree									-----> The "final" gene tree after treefix correction / internal node labelling by seadog mDGS reconciliation / branch length computing by PhyML
    │   ├── *_species.tree								-----> Species tree extracted from ncbi taxonomy
    │   ├── muscle_*.fasta								-----> Alignment file from muscle  
    │   ├── phyml_trimal_muscle_*.tree							-----> Tree generated by PhyML after alignment with MUSCLE and filtration with trimAl
    │   ├── phyml_trimal_muscle_*_treefix_dir						-----> Directory containing TreeFix tree correction files
    │   │   ├── phyml_trimal_muscle_*.fasta						-----> Sequence alignment used for tree correction
    │   │   ├── phyml_trimal_muscle_*.tree						-----> Initial tree generated by PhyML before correction
    │   │   └── phyml_trimal_muscle_*.treefix.tree					-----> Tree corrected after treatment with TreeFix
    │   ├── phyml_trimal_muscle_*_treefix_dir.smap					-----> Node correspondence file between initial tree and corrected tree
    │   ├── phyml_trimal_muscle_*_treefix_dir_treeToFixPath.txt				-----> Path of the file containing the tree to be corrected
    │   ├── trimal_muscle_*.fasta							-----> Sequence alignment after filtration with trimAl
    │   ├── trimal_muscle_*.phylip							-----> Phylip file containing filtered alignment
    │   ├── trimal_muscle_*.phylip_phyml_stats.txt					-----> Phylogenetic model statistics used by PhyML
    │   └── trimal_muscle_*.phylip_phyml_tree.txt					-----> Phylogenetic tree generated by PhyML from Phylip alignment
    │
    ├── gene_tree_*
    │   └── gene.tree									-----> The "final" gene tree after treefix correction / internal node labelling by seadog mDGS reconciliation / branch length computing by PhyML
    │
    ├── modulesChange_seadogMD_*.csv							-----> Lists of modules gained / lost at each gene nodes
    ├── modulesCompo_seadogMD_*.csv							-----> Lists of modules presents at each gene nodes
    ├── modules_paths_modules_*_t*m*M*.txt						-----> Define list of module trees to use as seadog input
    ├── modules_segm_dir_*								-----> All modules decomposition files: paloma plma -> for all modules: module fasta -> PhyML tree -> treefix corrected tree
    │   ├── *.dnd
    │   ├── *.fasta
    │   ├── *_t*m*M*.dot								-----> The plma file (from input or computed by paloma)
    │   ├── *_t*m*M*.oplma
    │   ├── *_t*m*M*.plma
    │   ├── *_t*m*M*_q*.agraph								-----> The plma file (from input or computed by paloma)
    │   ├── *_t*M*.afc
    │   ├── *.weights
    │   └── modules_*_t*m*M*
    │       ├── B*.fasta								-----> A module fasta file
    │       ├── B*.phylip								-----> A module in phylip format
    │       ├── B*.phylip_phyml_stats.txt						-----> PhyML output
    │       ├── B*.phylip_phyml_tree.txt						-----> PhyML output
    │       ├── B*.tree									-----> PhyML output
    │       ├── B*_treefix_dir								-----> treefix working directory
    │       │   ├── B*.fasta								-----> Module fasta file 
    │       │   ├── B*.tree								-----> Module PhyML tree
    │       │   └── B*.treefix.tree							-----> Module treefix corrected tree - "final"
    │       ├── B*_treefix_dir.smap							-----> Gene - Module nodes mapping
    │       └── B*_treefix_dir_treeToFixPath.txt					-----> Define treefix input
    │
    ├── domains_*.csv									-----> Known domains/motifs from scans
    ├── complete_functionChange_moduleChange_seadogMD_*.csv				-----> Table of all annotations and modules gained / lost at each gene nodes (actual and ancestral)
    ├── functionChange_moduleChange_expand_seadogMD_*.csv				-----> Table of modules gained / lost at each gene nodes where there is annotation changes (actual and ancestral) with all module segment details in leaf
    ├── functionChange_moduleChange_seadogMD_*.csv					-----> Table of modules gained / lost at each gene nodes where there is annotation changes (actual and ancestral)
    ├── seadogMD_*_gene.tree								-----> Gene tree from seadog mDGS, with internal gene labelled
    ├── seadogMD_*.output								-----> Seadog mDGS output file
    ├── seadogMD_*_sp_gene_event.csv							-----> Gene nodes event from Species - Gene reconciliation (e.g., Gene duplication, Speciation)
    └── specieGeneEvent_seadogMD_*.csv							-----> Gene nodes event from Species - Gene reconciliation (e.g., Gene duplication, Speciation)
```
</details>

# 🐳 Docker Image Usage

## Docker Installation
For ease of use, considering the various software and dependencies required, we strongly recommend using our Docker image. 
After installing [Docker](https://www.docker.com/get-started/), you can pull our Docker image using the following command:

```
docker pull ghcr.io/ocmalde/fuse-phylotree:v1.1.0
```

## Run the Docker Image

First, you need to create a container, linked to your workspace, and connect to it interactively:

```
docker run -it --entrypoint /bin/bash -v "$(pwd):/workspace" -w /workspace ghcr.io/ocmalde/fuse-phylotree:v1.0.0
```

Then you can simply execute the whole pipeline using:

```
python3 /fuse-phylotree/fuse-phylotree.py <sequences.fasta> <annotations.csv> <gene_tree.tree>
```
To quit it, simply type ```exit```

Example with test files provided in the Docker image (you need to connect to the container first):
```
cd test_dir/ && python3 ../fuse-phylotree/fuse-phylotree.py --infer_gene_tree INFER_GENE_TREE 712buddy37seq.fasta leaf_Manual_712.csv
```

If you already made a container, you can obtain its ```<CONTAINER ID>``` using ```docker ps -a```.
Using the ```<CONTAINER ID>```, you can connect to the existing container with:
```
docker start <CONTAINER ID> && docker attach <CONTAINER ID>
```

To use move file from/to the container, use:
```
docker cp <CONTAINER ID>:/path/in/container/ /path/in/local
```

# 🌀 Singularity Image Usage

## Singularity Installation
On most High Performance Computing clusters, you can directly pull our Singularity image using the following command:

```
singularity pull fuse-phylotree.sif https://github.com/OcMalde/fuse-phylotree/releases/download/V1.1.0/fuse-phylotree.sif
```

## Run the Singularity Image

```
singularity exec fuse-phylotree.sif python3 /fuse-phylotree/fuse-phylotree.py <sequences.fasta> <annotations.csv> <gene_tree.tree>
```



# To build a sequence dataset based on orthogroups

Using simple lists of RefSeq IDs proteins and a FUSE-PhyloTree Docker or Singularity image, you can easily generate a formatted FASTA file for the pipeline. This file will contain the longest protein sequence for each homolog (both orthologs and paralogs) in your input family. Simply run one the following script:

[make_orthogroup_fasta.sh](https://github.com/OcMalde/fuse-phylotree/blob/main/helper-scripts/make_orthogroup_fasta.sh)

```Usage: make_orthogroup_fasta.sh <id_seq_file> <name> <docker_image|singularity_image_path>```

> 📑 ```<id_seq_file>``` short example:
```txt
NP_001987.2
NP_006476.2
NP_006477.2
NP_006478.2
...
```
Refer to [this file](https://github.com/OcMalde/fuse-phylotree/blob/main/data/analyse_fibulin/advanced/human_fibulins_ids.txt) for the fibulin file example.


# Standalone Modules

## Modules
Specific [modules of the pipeline](https://github.com/OcMalde/fuse-phylotree/tree/main/fuse-phylotree) can be executed independently (see their ```--help``` for all usage details).

For example, to compute a phylogenetic tree using Muscle/Trimal/PhyML/Treefix: 
```
python3 /fuse-phylotree/gene_phylo.py <fasta_file> <species_tree>
```
To only compute the final integration module on computed files:
```
python3 /fuse-phylotree/integrate_3phylo.py <seadogMD.output> <gene_tree.tree> --pastml_tab <pastml_seadogMD_combined_ancestral_states.tab> --domains_csv <domains.csv>
```

<details><summary>Example: Running only the final integration step on the fibulin working directory</summary>
  
Move in the working directory ([available here](https://github.com/OcMalde/fuse-phylotree/tree/main/data/analyse_fibulin/run_singularity_fibulin/dir_fibuline_phylocharmod)):
            
```cd dir_fibuline_phylocharmod``` 
            
and execute the integration module with:
            
```
python3 /fuse-phylotree/integrate_3phylo.py seadogMD_fibulin59.output gene_tree_fibulin59/fibulin59.tree --pastml_tab acs_dir_seadogMD_fibulin59_gene/pastml_seadogMD_fibulin59_gene_leaf_Manual_fibulin59_combined_ancestral_states.tab --domains_csv domains_fibulin59.csv
```
 
</details>

## Softwares
All the different [included softwares](https://github.com/OcMalde/fuse-phylotree/tree/main#softwares) are usable using the Docker or Singularity image.
For example, paloma-D can be used with:
```
docker start <CONTAINER ID> && docker exec <CONTAINER ID> /bin/bash -c ". ~/.bashrc && paloma-D --help && exit"
``` 
Or simply: ```paloma-D --help``` if connected to a container

```
singularity exec fuse-phylotree.sif paloma-D --help
``` 


# Dependencies

All these programs are mandatory to run *FUSE-PhyloTree*, and must be in ;
```/usr/local/bin/```
If not, their path must be specified in the config file ;
```fuse-phylotree/config.txt```

### Softwares

[Muscle](http://www.drive5.com/muscle/) (v3.8.31), one of the best-performing multiple alignment programs, [Conda package](https://anaconda.org/bioconda/muscle)

[PhyML](https://github.com/stephaneguindon/phyml) (v3.3.20190909), maximum likelihood phylogenetic inference for the gene and the module trees, [Conda package](https://anaconda.org/bioconda/phyml)

[TreeFix](https://www.cs.hmc.edu/~yjw/software/treefix/) (v1.1.10), Statistically Informed Gene Tree Error Correction Using Species Trees, [Conda package](https://anaconda.org/OcMalde/treefix)

[trimAL](https://github.com/inab/trimal) (v1.4.rev15), Alignment filtering, [Conda package](https://anaconda.org/bioconda/trimal) 

[SEADOG-MD](https://compbio.engr.uconn.edu/software/seadog/), for DGS-reconciliation.

[PastML](https://pastml.pasteur.fr/) (v1.9.41), for ancestral characters inference, [Pip package](https://pypi.org/project/pastml/)

[paloma-D](http://tools.genouest.org/tools/protomata/learn/) (v0.1), for sequence segmentation. 

