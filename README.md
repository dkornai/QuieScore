# QuieScore
The QuieScore R package is meant to help identify tomour quiescence, and aid in the classification of quiescence subtypes, as described in [Wiecek et al, Genome Biol 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02963-4) 

# Authors 
Daniel Kornai, Maria Secrier

# Installation
The package can be installed by using "devtools".

# How to run

For a sample script and input data, please see the **sampleScript** folder.

# Output

**q_score_raw** - the raw G0 arrest scores (anything above 0 would be arrested in G0, and if you want to be more stringent you can use a cut-off of 1.5 or 3)
**q_score_normalised** - aligns the scores with respect to TCGA samples from the same cancer type (NB this is only suitable when looking at bulk cancer samples, and only for cancers that are profiled in TCGA)
