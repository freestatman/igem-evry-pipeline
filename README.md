This project was created during the 2015 edition of the iGEM 
competition by the Evry iGEM team.

Our idea was to develop a pipeline allowing us to select the best 
candidates tumoral antigens to use for a vaccine (immunotherapy). 
A good candidate must be tumor-specific and sufficiently expressed 
in tumor cells to be presented to the immune system. Furthermore, 
it must be able to be processed efficiently by the immune system.

We intended to identify relevant targets by:
    - Looking for differentially expressed genes (specifically 
      upregulated genes) in tumoral tissue (vs. normal tissue), 
      by a transcriptomic analysis (such as RNASeq)
    - Looking for genetic mutations only found in tumoral tissue 
      (identifying genetic variants)

Once the targets are identified, the goal is not to express the whole 
corresponding genes but to express short fragments, corresponding to 
putative cleavage sites by the proteasome to link to the MHC-I in 
order to potentially triger an tumor-specific immune response.

We created the whole pipeline using data corresponding to RNAseq reads of melanocytes cell lines. 
All the figures and data presented in the software section were obtained using this dataset.

Experimental procedures for data generation are described in: 
  Vardabasso, C. et al. (2015). Histone Variant H2A.Z.2 Mediates 
  Proliferation and Drug Sensitivity of Malignant Melanoma. Molecular Cell, 59:75-88

Data was accessed from the Sequence Read Archive website (NCBI). Study identifier: SRP057616
You can find it here: http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP057616

Detailed information about each script can be found in a dedicated README file.
