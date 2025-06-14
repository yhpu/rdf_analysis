# RDF Analysis

RDF Analysis is a comprehensive tool designed to process output data from various molecular simulation software.  
It calculates the radial distribution function (RDF) to characterize molecular structures and interactions,  
and generates detailed charts and data outputs for further analysis.

---

## Features

- Supports multiple molecular simulation output formats  
- Calculates radial distribution functions (RDF) accurately  
- Outputs both raw data files (.dat) and corresponding figures 

---

## Installation

It is recommended to create a dedicated Conda environment:

```bash
conda create --name rdf python=3.10
conda activate rdf
```

## Usage
Run the program with your trajectory file as input:

```bash
python main.py {Trajectory_File}
```

## Test
``` 
pip install pytest
pytest tests
```