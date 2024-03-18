# nf-core/yascp: Comprehensive Documentation

Welcome to the documentation for `nf-core/yascp`, a versatile pipeline designed to address a broad spectrum of research questions in the field. The pipeline's flexibility and array of features make it a powerful tool, though its complexity can seem daunting at first. Our documentation is organized into specific sections to help you get started, understand the pipeline's output, and explore various usage scenarios.

### Getting Started

- **[Usage](usage.md):** This section offers a detailed guide on operating the pipeline, including an overview of its functionality, instructions for execution, and explanations of command-line options. Whether you're a novice or an experienced user, this guide is designed to provide a solid foundation for running `nf-core/yascp` effectively.

### Understanding Your Results

- **[Output](output.md):** Dive into the results produced by the pipeline with this comprehensive overview. Learn how to interpret the different outputs and the insights they offer into your data.

### Detailed Tutorials

For users looking to tailor the pipeline to specific needs, we offer a series of tutorials covering a wide range of use cases:

- **[Running the Full Pipeline](full_pipeline_tutorial.md):** A step-by-step guide to executing the complete workflow.
- **[Excluding Cellbender Ambient RNA Removal](no_cb_full_pipeline.md):** Opt for this tutorial if GPU resources are unavailable.
- **[Full pipeline with available Genotypes](full_pipeline_GT__tutorial.md):**
- **Focusing on Specific Components:** Tailored instructions for running only certain parts of the pipeline:
    - [Celltype Assignment](celltype.md)
    - [Cellbender Ambient RNA Removal](ambient_rna_removal_tutorial)
    - [Doublet Detection](doublet_detection_tutorial.md)
    - [Integration, Clustering, and Cluster Assignments](cluster_integrate_tutorial.md)
    - [GT (Genotype) Matching](gt_match_tutorial.md)
    - [Cleaning Up Result](clean_up_results.md) Best practice for ensuring that unnecessary intermediate files dont occupy a lot of disk space.

### Additional Resources

For further information on installation, configuration, and general usage of `nf-core` pipelines, please visit our official website: [https://nf-co.re](https://nf-co.re).
