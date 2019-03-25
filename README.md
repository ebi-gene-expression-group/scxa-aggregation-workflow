# Single-cell Expression Atlas aggregation workflow

Aggregating matrices from many thousands of quantification runs is a tricky task, involving as it can tens of thousands of data points for hundreds of thousands of cells.

We use the module [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) to parse the results of methods such as Kalliso, since it provides some useful scaling options. Unfortunately it's not suitable for parsing many thousands of runs at once, and memory usage quickly gets out of hand.

The solution to this conundrum is to run tximport on small chunks of runs at a time, and output a sparse matrix from each. These separate sparse matrices can then be merged with a relatively small memory footprint. This is the process implemented in this workflow. At the moment this just works for Kallisto outputs, but there's not reason it could not be parameterised to other methods when need arises.

## Setup

### Conda/ Bioconda

Workflow dependencies are managed via Conda and Bioconda, so you'll need to set that up, see instructions [here](https://bioconda.github.io/#install-conda). 

### Nextflow

Obviously you'll need Nexflow itself. If you don't have it already you can install via Conda:

```
conda install nextflow
```

You may well want to do this within a Conda environment you create for the purpose.

## Run the workflow

### Inputs

Expected inputs are:

* A reference GTF file, used to generate a transcript-to-gene mapping
* A root directory in which results from Kallisto can be found in subdirectories. 

### Parameters
 
You can copy the default configuration, edit the parameters, and provide it to Nextflow to override any of the settings. See the [Nexflow documentation](https://www.nextflow.io/docs/latest/executor.html) for executor settings.
 
### Execution

The workflow can be run directly from the repository:

```
nextflow run -config <your nextflow.config> ebi-gene-expression-group/scxa-aggregation-workflow --resultsRoot <results dir> --referenceGtf <reference GTF file> 
```

This will download the workflow, create any necessary environments, and run the workflow with the specified inputs. Input files and subdirectories are expected to exist in the specified results directory. Output will be to a directory called 'bundle' in the same location.

Future executions will use a cached copy of the pipeline, should you wish to update the code in future, you can do so like:

```
nextflow pull ebi-gene-expression-group/scxa-aggregation-workflow
```

