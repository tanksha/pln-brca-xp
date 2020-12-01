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
