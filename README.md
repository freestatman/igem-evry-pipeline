# Patient-specific neo-antigen identification pipeline

This project was created during the 2015 edition of the iGEM competition by the Evry iGEM team. More information the general project can be found on the [wiki](http://2015.igem.org/Team:Evry).


## Principle
Our idea was to develop a pipeline allowing us to select the best
candidates tumoral antigens to use for a vaccine (immunotherapy).
A good candidate must be tumor-specific and sufficiently expressed
in tumor cells to be presented to the immune system. Furthermore,
it must be able to be processed efficiently by the immune system.

We intended to identify relevant targets by:
* Looking for differentially expressed genes (specifically upregulated genes) in tumoral tissue (vs. normal tissue), by a transcriptomic analysis (such as RNASeq)
* Looking for genetic mutations only found in tumoral tissue (identifying genetic variants)

Once the targets are identified, the goal is not to express the whole
corresponding genes in our system but only short fragments that will trigger an immune response targeted to the cells expressing the 'tumor-specific' antigens. These antigens correspond to short peptides that should be processed efficiently by the immune system (i.e. cleaved by the proteasome and able to bind to the MHC-I) in order to potentially trigger an tumor-specific immune response.


## Used data
We created the whole pipeline using data corresponding to RNAseq reads of melanocytes cell lines. All the figures and data presented in the [software section](http://2015.igem.org/Team:Evry/Software) of our project wiki were obtained using this dataset.

Experimental procedures for data generation are described in:
  _Vardabasso, C. et al. (2015)._ Histone Variant H2A.Z.2 Mediates
  Proliferation and Drug Sensitivity of Malignant Melanoma. _Molecular Cell, 59:75-88_

Data was accessed from the Sequence Read Archive website (NCBI). Study identifier: SRP057616
You can find it [here](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP057616) (SRA, NCBI).


## Documentation

Detailed information about each script of the pipeline can be found in a dedicated README file. Otherwise, more general information about the pipeline found on [our iGEM project's wiki](http://2015.igem.org/Team:Evry/Software).
