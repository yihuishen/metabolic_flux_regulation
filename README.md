# metabolic flux regulation and proteome resource allocation
This repository provides Supplementary files and R scripts regarding multi-omics data integration that are part of the manuscripts:

**Proteome capacity constraints favor respiratory ATP generation**
Yihui Shen, Hoang V. Dinh, Edward Cruz, et al. submitted

and part of the preprint (https://doi.org/10.1101/2022.08.10.503479)

## file descriptions
### 1) Isor_SIMMER
Systematic Identification of Meaningful Metabolic Enzyme Regulation in _Issatchenkia orientalis_<br>
integration of fluxomics, metabolomics, and proteomics via Michaelis-Menten kinetics<br>
see also Hackett, S. R. et al. Systems-level analysis of mechanisms regulating yeast metabolic flux. Science 354, aaf2786–aaf2786 (2016).<br>

### 2) Sace_Isor_flux_regulation
Explaining flux difference between yeasts (_S. cerevisiae_ and _I. orientalis_) or flux change across nutrient-limited chemostat conditions

**Dataset included**<br>
    fluxomics (from 13C MFA), LC-MS metabolomics, quantitative proteomics<br>
        _S. cerevisiae_: 12 nutrient-limited chemostats + 1 batch culture<br>
        _I. orientalis_: 15 nutrient-limited chemostats + 1 batch culture<br>
    
### 3) mmTcell_flux_regulation
Explaining flux change upon T cell activation

**Dataset included**<br>
    fluxomics (from 13C MFA), LC-MS metabolomics, quantitative proteomics<br>
        naive and activated mouse CD8+ T cells<br>

### 4) protein_resource_allocation_yeast
Coarse-grained analysis of proteome allocation in _S. cerevisiae_ and _I.orientalis_

**Dataset included**<br>
    protein abundance of _S. cerevisiae_ and _I.orientalis_<br>
        batch culture and nutrient-limited chemostats<br>
        aerobic (+O2), anaerobic (-O2), and with antimycin treatment (+antimycin)<br>
    functional assignment of the proteome
    comparison to reference proteomics from literature

### 5) protein_resource_allocation_mouse
Coarse-grained analysis of proteome allocation in mouse cells and tissues, and proteome efficiency of ATP generation pathways of all systems described in the manuscript

**Dataset included**

    protein abundance of naive and activated mouse CD8+ T cells
    protein abundance of mouse tissues and tumors
        healthy pancreas, GEMM PDAC, flank PDAC
        healthy spleen, leukemic spleen
    functional assignment of mouse proteome
