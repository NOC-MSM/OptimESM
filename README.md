# OptimESM
Post-Processing &amp; Analysis Workflows for NOC contributions to the [**OptimESM project**](https://optimesm-he.eu).

---

### Description
The primary goal of **OptimESM** is to develop the next generation of ESMs, bringing together increased model resolution and process realism, and to deliver long-term climate projections that better support policy and societal needs, providing guidance on regional climate change at different levels of global warming, the risk of abrupt Earth system changes at these warming levels and the regional impacts arising from such events.

---

### Getting Started

#### **Downloading OptimESM**

To get started, clone this repository to your local machine or HPC system via GitHub as follows.

```bash
git clone git@github.com:NOC-MSM/OptimESM.git
```

#### **Creating OptimESM Virtual Environment**

If you do not already have an installation of miniconda or miniforge on your local machine, see the instructions below to install miniforge3.

The `OptimESM` GitHub repository includes the `environment.yaml` file required to reproduce the Conda virtual environment needed to run OptimESM diagnostics and perform tipping point analysis.

From your Conda base environment, run the following command to create a new Python virtual environment using the `.yaml` file:

```bash
conda env create --name env_optimesm --file environment.yaml
```

This will take a bit of time to install all of the dependencies from the **conda-forge** channel and then the **nemo_cookbook** package using `pip`.

#### Downloading METRIC
Once you have created an OptimESM virtual environment, next download and use pip to run an editable install of the METRIC package used in WP4:

Assuming your virtual environment is activated and you are in the root directory of your cloned copy of the `OptimESM` repository:
```bash
cd WP4/
git clone git@github.com/AMOCcommunity/metric.git
cd metric
pip install -e .
```

#### Downloading NEMO Pipeline
Finally, we need to download and use pip to run an editable install of the NEMO Pipeline package used in WP5:

Again, assuming your virtual environment is activated and you are in the root directory of your cloned copy of the `OptimESM` repository:
```bash
cd WP5/
git clone git@github.com:NOC-MSM/nemo_pipeline.git
cd nemo_pipeline
pip install -e .
```
---

### Structure

`/data`
* The `/data` directory includes symbolic links to all ESM CMORISED outputs available on the CINECA HPC for:
  - **UKESM1-2-LL**
  - **EC-Earth3-ESM-1**
  - **IPSL-CM6-ESMCO2**
  - **CNRM-ESM2-1**

`/data/{ESM}/Ofx`
* For each of the ESMs above, the `Ofx` directory contains NEMO ocean model `domain_cfg.nc`, `mesh_mask.nc` and regional mask files used to create `NEMODataTrees` and perform diagnostic calculations.

`/WP4`
> WP4 is evaluating the ESMs against observational data. It is compiling novel observational datasets to understand the current status of models. Metrics is developed to evaluate how ESMs simulate the mean state, variability, ongoing changes in the Earth system and particular processes relevant to abrupt change and tipping points.

  * `/WP4/metric`
    - METRIC ([https://github.com/AMOCcommunity/metric](https://github.com/AMOCcommunity/metric)) open-source Python library to perform observational-equivalent AMOC diagnostics.
  * `/WP4/configs`
    - METRIC configuration files for UKESM1-2-LL, EC-Earth3-ESM-1, IPSL-CM6-CO2, CNRM-ESM2-1.
  * `/WP4/data`
    - AMOC + MHT + MFT outputs from METRIC for OptimESM idealised, historical and piControl experiments. 
  * `/WP4/figures`
    - Figures for ocean circulation validation.
  * `/WP4/analysis`
    - Jupyter Notebooks + Scripts to perform ocean circulation validation.

`/WP5`
> WP5 is exploring abrupt changes and tipping points in existing CMIP6 runs and in the new OptimESM-runs; developing new techniques to detect and forewarn of tipping points under transient climate change; and developing a risk landscape to summarise how tipping point risks vary with the rate of global warming. The WP is prioritising potential tipping points that may be triggered before the year 2300 and abrupt changes that could have significant impacts on human and natural systems.

  * `/WP5/nemo_pipeline`
    - NEMO Pipeline ([https://github.com/NOC-MSM/nemo_pipeline](https://github.com/NOC-MSM/nemo_pipeline)) open-source Python library to create reproducible diagnostic pipelines for NEMO OGCM outputs.
  * `/WP5/configs`
    - NEMO Pipeline configuration files for UKESM1-2-LL, EC-Earth3-ESM-1, IPSL-CM6-CO2 diagnostics.
  * `/WP5/data`
    - NEMO Pipeline diagnostic outputs for OptimESM idealised, historical and piControl experiments. 
  * `/WP5/figures`
    - Figures for ocean circulation + sea ice tipping point analysis.
  * `/WP5/src`
    - Jupyter Notebooks + Scripts to perform pre/post-processing and analysis of NEMO Pipeline ocean and sea-ice diagnostics.

`/WP6`
> WP6 is investigating the regional consequences of reaching different warming levels or from the triggering of abrupt climate changes, with a focus on Europe and the polar regions. Changes in extreme event occurrence and shifts in climate types over Europe is attributed to different warming levels. Furthermore, ML based statistical downscaling methods is developed for specific variables to provide high-resolution climate information from our global simulations. The added value of these new tools is assessed in pilot applications over Europe and the downscaled data disseminated to the international impacts community enabling a wider assessment of the regional impacts.

---

### Resources
#### **Installing Miniforge**

**Miniforge** is a minimal conda installer that defaults to the **conda-forge** ecosystem. It is lightweight, open-source, and recommended for scientific Python workflows.

**Download Miniforge**

Follow instructions at [**https://github.com/conda-forge/miniforge**](https://github.com/conda-forge/miniforge).

For Unix-like Platforms (macOS, Linux, WSL): 

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

**Creating the Base Environment**

Run the installer by running the following command:

```bash
bash Miniforge3-$(uname)-$(uname -m).sh
```

It is recommended to accept the default location (`~/miniforge3`) and only accept the changes to your `~/bashrc` file if you are happy for the base environment to be activated every time you log in (this is not recommended for JASMIN use).

**Activating the Base Environment**

Once you have successfully installed miniforge3, activate your base environment:

```bash
source ~/miniforge3/bin/activate 
```
