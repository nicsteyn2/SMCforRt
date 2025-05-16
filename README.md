# SMC and epidemic renewal models

Renewal models are popular in statistical epidemiology for their use as a semi-mechanistic model of disease transmission. We demonstrate how sequential Monte Carlo (SMC) methods can be used to perform inference on and generate projections with these models.

Our approach:

👍 Is conceptually simple  
👍 Flexible  
👍 Produces valid credible intervals  
👍 Can accommodate underreporting, missing/incomplete data, imported cases, uncertain serial intervals, temporally aggregated data, multiple data sources, etc - all at the same time!

This repository contains all data, scripts, and notebooks required to follow along with the [primer](https://doi.org/10.48550/arXiv.2503.18875).

We also highlight the accompanying website, [An introduction to sequential Monte Carlo methods for reproduction number estimation](https://nicsteyn2.github.io/SMCforRt/).

The primer and website are designed for researchers with a basic understanding of statistical inference. We walk through the steps of constructing an SMC model to build a gold-standard $R_t$ estimator, fine-tuned for your purposes, from scratch. This balances the convenience of an off-the-shelf method, with the flexibility and deep understanding that comes with constructing a model from first principles.

![Figure showing unknown Rt and reported cases in NZ](assets/root_readme_data.png)


## Getting started

All code for both the [primer](https://doi.org/10.48550/arXiv.2503.18875) and the [website](https://nicsteyn2.github.io/SMCforRt/) is available as interactive notebooks in this repository. We encourage you to check them out on [Binder](https://mybinder.org/v2/gh/nicsteyn2/SMCforRt/main), a free cloud-based Jupyter notebook service.

Binder can take a few minutes to start up, but saves you the hassle of downloading the repository and installing Julia (if you do not already have it). As this is a free service, there are substantial computational limitations, but we have succeeded in running the simpler models (such as visualising the data and running Model 1 from the primer).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nicsteyn2/SMCforRt/main)

To run the more computationally intensive models, we recommend running the code locally. Fortunately, as there are no external dependencies, this is very easy to do. To get started, all you need to do is:

1. Install Julia from [julialang.org](https://julialang.org/downloads/)
2. Clone this repository
3. Run Julia in your terminal and run:

```julia
using Pkg
Pkg.add("IJulia")
using IJulia
notebook()
```

This will open a Jupyter notebook in your browser. Navigate to the local copy of the GitHub repository and open `gettingstarted.ipynb` to get started.

**Alternatively**, you can run the notebooks themselves in VS Code (or your preferred IDE). The VS Code Julia extension is helpful for this, as it provides a nice interface for running Julia code. You can install the extension from the [VS Code Marketplace](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia).

Opening the root directory of this repository in VSCode and then opening any of the notebooks should work out of the box. A short script to install the required packages is included in `gettingstarted.ipynb`.

## Structure of this repository

The files and folders you should care about:
- `gettingstarted.ipynb`: A jupyter notebook to help with setup and navigation
- `/data/`: Contains all example data used in the tutorial
- `/src/`: Contains crticial source code
- `/paper/`: Contains the code, in the form of notebooks, to reproduce the examples in the paper
- `/notebooks/`: Contains the main tutorial notebooks

You can ignore (but feel free to explore):
- `/docs/`: Contains the rendered tutorial
- `/site/`: Contains the quarto website files (these are rendered to `/docs/`)
- `/assets/`: Contains images and other assets used in the tutorials and readme files
- Files like `.gitignore` and `.nojekyll` which are for repo management

