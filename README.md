# Breast Cancer PLN Experiment

This repo contains the code for inference experiments using Probabilistic Logic Networks, PLN, on breast cancer dataset. It is based on the work in [singnet/cancer](https://github.com/singnet/cancer)

### Prerequisites

-  [AtomSpace](https://github.com/singnet/atomspace) - we use the AtomSpace to store the patient, gene expression, biological process etc in hypergraph form.
- [PLN](https://github.com/singnet/pln/) - used for inference, deducing patient, [GO](http://geneontology.org/docs/ontology-documentation/) biological process relationships and the intentional similarity between pathways and biological processes
- [Unified Rule Engine](https://github.com/singnet/ure/) - This is an indirect dependency used by PLN.

- [Knowledge Import scripts](https://github.com/mozi-ai/knowledge-import) - used to convert Patient gene expression table to Atomese

### Building & Installation

1. Clone the repo
2. Run `autoreconf -vif` inside the project directory
3. Run the following in the project directory:
    ```sh
    $ mdkir build && cd build
    $ ../configure
    $ make && makie install
    ```

### Datasets

**Note**: All the datasets are assumed to be in `kbs` directory under the root directory of the project. Also create `results` directory to save the results

- GO Files: [GO](https://mozi.ai/datasets/current_2020-10-20/GO_2020-10-21.scm), [GO_annotation](https://mozi.ai/datasets/current_2020-10-20/GO_annotation_2020-10-20.scm)

- Pathway files: [Reactome](https://mozi.ai/datasets/current_2020-10-20/reactome_2020-10-20.scm) , [Reactome_annotation](https://mozi.ai/datasets/current_2020-10-20/NCBI2Reactome_PE_Pathway.txt_2020-10-20.scm)

- [Patient gene expressions](https://mozi.ai/datasets/cancer/patient_gene_expression.zip)


### Running the code

#### 1. Similarity Experiment
- Start `Guile` and run the following 
```sh 
    scheme@(guile-user)> (use-modules (opencog) (opencog exec)
    (opencog bioscience) (opencog ure) (opencog pln) 
    (opencog logger))

    scheme@(guile-user)> (use-modules 
    (pln-bio preprocess) (pln-bio bio-utils) 
    (pln-bio intensional-similarity))

    scheme@(guile-user)> (go-pathway-intentional-similarity (list "kbs/GO_2020-10-20.scm" "kbs/GO_annotation_2020-10-20.scm"
 "kbs/reactome_2020-10-20.scm" "kbs/NCBI2Reactome_PE_Pathway.txt_2020-10-20.scms"))
```


#### 2. Patient Biological process deduction
- Start `Guile` and run the following 
```sh 
    scheme@(guile-user)> (use-modules (opencog) (opencog exec)
    (opencog bioscience) (opencog ure) (opencog pln) 
    (opencog logger))

    scheme@(guile-user)> (use-modules (pln-bio bio-utils) (pln-bio expr) (pln-bio bp-deduction))

    ;;set #t for overexpression #f for underexpression
    scheme@(guile-user)> (run-deduction-expr #t)
```