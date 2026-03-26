# Splice-Altering Variant Resources

A summary of software tools and databases used when identifying and interpreting intronic splice-altering variants.

---

## Software Tools

### [SpliceAI](https://github.com/illumina/spliceAI)

SpliceAI is a deep learning–based tool that predicts the impact of variants on
splicing. It outputs four delta scores per variant:

| Score | Code | Description |
|---|---|---|
| Acceptor Gain | DS_AG | Probability of an acceptor splice site being created/strengthened |
| Acceptor Loss | DS_AL | Probability of an acceptor splice site being lost/weakened |
| Donor Gain | DS_DG | Probability of a donor splice site being created/strengthened |
| Donor Loss | DS_DL | Probability of a donor splice site being lost/weakened |

Scores range from 0.0 – 1.0. Higher scores indicate a stronger predicted impact on splicing. 

> ### Recommended Score Thresholds:    
> 0.20 (High Recall): 🔍 Used for capturing most splicing-altering variants.    
> 0.50 (General/Default): ⚖️ Default used for balanced sensitivity and precision.    
> 0.80 (High Precision): 🎯 Used for identifying high-confidence splice-altering variants.    

> #### SpliceAI-10K Calculator Thresholds (more info below):
> 0.05 (Deep Intronic/Pseudoexon): 🔬 Specific contexts often use lower thresholds, such as 0.05 for pseudoexon gain.   

> **Note:** SpliceAI annotations are now incorporated into [nf-cavalier (dev branch)](https://github.com/bahlolab/nf-cavalier/tree/dev)  
> By default, a minimum score of **≥ 0.20** is used to flag variants for review.  
> Parameter = 'FILTER_SHORT_MIN_SPLICEAI_PP'  
---

### [Pangolin](https://github.com/tkzeng/Pangolin)

Pangolin is a complementary splicing prediction tool that uses a tissue-aware model to predict splice site usage across multiple tissues.    
It provides an additional layer of evidence when evaluating candidate splice-altering variants.

---

## Types of Splice Impacts

### A) Exon Skipping

Exon skipping occurs when a variant weakens the donor AND acceptor splice sites for an annotated exon,    
causing the spliceosome to skip the affected exon entirely during pre-mRNA processing.

**[Click here for scores for an example variant that causes exon skipping](https://spliceailookup.broadinstitute.org/#variant=chr8-140300616-T-G&hg=38&bc=basic&distance=500&mask=0&ra=0) (TRAPPC9:c.1623-2A>C)**

### B) Intron Retention

Intron retention arises when a splice site is weakened, leading to either partial or complete retention of the intron.    
Partial retention occurs when an adjacent splice site is activated, whereas complete retention results when no alternative site is used.    
Both forms can coexist within a transcript, with their relative abundance determined by factors such as the strength of the affected splice site,    
the availability of nearby cryptic sites, and the integrity of local splicing regulatory elements.

**[Click here for scores for an example variant that causes intron retention](https://spliceailookup.broadinstitute.org/#variant=MLH1%3Ac.1667G%3EC&hg=38&bc=basic&distance=500&mask=0&ra=0) (MLH1:c.1667G>C)** 

### C) Pseudoexon Retention

Pseudoexon retention occurs when a deep intronic variant activates a cryptic splice acceptor and/or donor,     
causing a stretch of intronic sequence to be incorporated into the mature transcript as a novel exon.

**[Click here for scores for an example variant that causes pseudoexon retention](https://spliceailookup.broadinstitute.org/#variant=DNM1L%3Ac.456%2B215A%3EG&hg=38&bc=basic&distance=500&mask=0&ra=0) (DNM1L:c.456+215A>G)** 

**Example SpliceAI scores — `DNM1L:c.456+215A>G (chr12:32711230:A>G)`**

| Delta Type | Delta Score | Position |
|---|---|---|
| Acceptor Loss | 0.00 | — |
| Donor Loss | 0.00 | — |
| Acceptor Gain | 0.11 | −50 bp |
| Donor Gain | 0.10 | +42 bp |

### Pseudoexon Filtering Criteria
Score threshold can be lowered to **0.05** when specifically searching for pseudoexon-activating variants, as recommended by the [SpliceAI-10K Calculator](https://github.com/adavi4/SAI-10k-calc).


For a variant to be considered a strong pseudoexon candidate, the following criteria
should ideally be satisfied:

| Criterion | Threshold | DNM1L:c.456+215A>G |
|---|---|---|
| Min. Gain score (AG or DG) | ≥ 0.02 | ✅ DG = 0.10 |
| Max. Gain score (AG or DG) | ≥ 0.05 | ✅ AG = 0.11 |
| AG position relative to DG (strand-aware) | AG < DG | ✅ −50 bp < +42 bp |
| Gained exon size | 25 bp ≤ size ≤ 500 bp | ✅ 93 bp |
| Variant distance from nearest canonical junction | > 50 bp | ✅ Deep intronic (+215) |
| AG & DG positions from nearest canonical junction | > 50 bp each | ✅ |
| AG & DG positions | Within the same intron | ✅ Intron 5 |

> **TODO:** Add more complex SAI-10K rules for more accurate pseudoexon filtering to nf-cavalier.

---

## Visualisation — [SpliceAI Lookup](https://spliceailookup.broadinstitute.org/)

SpliceAI (and Pangolin) predictions can be visualised via the **SpliceAI Lookup**
tool, hosted by the Broad Institute.   
   
1. Enter your variant of interest to see a summary
of delta scores   
2. Click **Show** in the **Visualizations** section to view
predicted splice impacts relative to annotated gene transcripts.

---

## Pseudoexon Sequence Extraction

Please try the script via the Google Colab notebook here:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1VbUHMF_GSOyvUH7g0bSoxf3VLXiGv4JK?usp=sharing)

Once a pseudoexon candidate has been identified, run `extract_pseudoexon_sequence.py` (script in this repository)    
to extract the predicted retained intronic sequence and check whether it introduces a **premature termination codon (PTC)** or causes a frameshift.

**Usage:**
```bash
python extract_pseudoexon_sequence.py chr12:32711230 --ag -50 --dg 42
```

> **Note:** 
> 1. Requires pysam, and biopython (pip install pysam biopython). Currently configured for the hg38 human genome assembly; modify accordingly for other builds.     
> 2. `--ag` and `--dg` accept SpliceAI offset values (in bp) relative to the variant position

**Example output — `DNM1L:c.456+215A>G`:**

```
Input: chr12:32711230  AG=-50  DG=42
Computed retained intron range: chr12:32711180-32711272 (n=93 bp)

> Region: chr12:32711180-32711272 (plus strand)  (length = 93 bp)
Sequence: GATACCACACTCTACAGGTTTATCTCCAGTCTCACTGGCCACTCCTTTTCAGTCTCTTTTGCTGTTTCTTCTCATCTCCCCAATGTCTAAACA

Amino acid length  : 31 aa (full codons: 31)
In-frame?          : in-frame (no partial codon at end)
Contains stop (*)  : True
AA (first 60 aa)   : DTTLYRFISSLTGHSFSVSFAVSSHLPNV*T
```

In this example, the retained pseudoexon is **in-frame** and contains a stop codon. This is likely to trigger **nonsense-mediated decay (NMD)** of the aberrant transcript.

---

## [SpliceVault](https://kidsneuro.shinyapps.io/splicevault/) - Unannotated Splice Junction Database

SpliceVault catalogues low-level, unannotated splice junctions observed in a large
control cohort derived from GTEx and other RNA-seq resources.   

A key principle when interpreting splice-altering variants is that they typically do **not** create
entirely novel splice sites from scratch. Instead, they **amplify the use of
pre-existing cryptic splice events** — naturally occurring splicing "mistakes" already
detectable at trace levels in healthy tissue.   

A pathogenic variant near these latent sites shifts the equilibrium away from canonical splicing,   
icreasing usage of the cryptic junction. 

SpliceVault allows you to confirm whether a predicted cryptic site has any baseline support in normal transcriptomes,    
providing biological plausibility for the predicted event.    

### Worked Example — `DNM1L:c.456+215A>G`

**Variant position:** `chr12:32711230`

**Predicted sites (from SpliceAI):**
- Acceptor Gain: variant − 50 bp → `chr12:32711180`
- Donor Gain: variant + 42 bp → `chr12:32711272`

**Steps:**
1. Enter the variant using its transcript notation:
   `NM_012062.5:c.456+215A>G`
   *(SpliceVault requires an Ensembl or RefSeq transcript identifier)*
2. Click **Show all Events**, then **Generate Table**
3. Search the output for your predicted junction coordinates
   (e.g. `Ctrl+F` → `32711272`)

> **Coordinate note (±1 offset):** A 1 bp discrepancy between SpliceAI and SpliceVault
> positions is expected and does not indicate a mismatch. SpliceAI reports the predicted
> splice site position anchored to the **last exonic base** (for donors) or the **first
> exonic base** (for acceptors), using 1-based coordinates. SpliceVault junction data is
> derived from RNA-seq aligners (e.g. STAR), which report the **first intronic base** as
> the donor boundary — one position downstream. This means a SpliceAI donor prediction
> at `32711272` and a SpliceVault junction at `32711273` are describing the same splice
> event. Always treat positions within ±1 bp as equivalent matches.

---

## [UCSC Genome Browser](https://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=hg38) — Alternate Transcript Support

Search the UCSC Genome Browser to determine whether any alternate annotated transcripts share the     
predicted pseudoexon boundaries, providing independent support for thecryptic site's accessibility.    

**Steps — `DNM1L:c.456+215A>G`:**
1. Search coordinates: `chr12:32711180-32711272`
2. Open **GENCODE V49** (or latest available) track settings
3. Set **Display Mode** to `Full`
4. Enable: **Non-coding genes**, **Splice variants**, **Pseudogenes**
5. Set **Tagged Sets** to `All`
6. Click **Submit**

If any annotated isoforms share one or both predicted pseudoexon junctions, this constitutes independent evidence that the spliceosome     
is capable of recognising that intronic position.
