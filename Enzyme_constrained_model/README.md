# EC-model for S-Clav

The latest updated version of the iDG1237 S. Clavuligerus model was used to built an enzyme-constrained model.
The basis for this model are: kcat-values for each enzyme, MW-values for each enzyme, and a total proteome constraint.

We adopted the Ecmpy methodogy for this:
https://github.com/tibbdc/ECMpy

Next to that we added the option to split the proteome pool in multiple pools, according to custom input.

This read-me will explain how the model was built and how to use it

# Model building
The Ecmpy-version of the updated iDG1237 model was built by adapting the scripts from the Ecmpy github.
Step 1:
NB! Installation of Emcpy is ONLY needed when building a new model. When you want to use an existing ecmpy-model,
you only need to install Cobrapy.

Installation of Ecmpy, following the documentation instruction.
NB: ALWAYS install tools as Ecmpy in a venv or a conda-environment, the dependencies are very specific.

However, some issues were presented, issues could be overcome by:

- Installation of cobrapy via conda-forge (before installing Ecmpy via pip)
- Python environment of 3.11
- Downgrade the version of setuptools: `pip install --force-reinstall -v "setuptools==75.1.0"`
- Downgrade torch to v2.2.2, and import the torch module before any other modules in the main script
- Upgrade numpy

NB: the resulting solution does lead to userwarning (by Numpy), but they do not stop the scripts from running.


The code is described in:
Enzyme_constrained_model/code/run_ecmpy_for_streptomyces.py
NB: most steps were run sequential, but they only needed to be run once. Therefore, many parts are commented out.
In order to reproduce the script, all steps in the main function should be run sequentially.

The helper function: ECMpy_function.py was taken from the Ecmpy-repo. It required some custom adaptations to fix
some specific errors (so it is updated compared to the same script on the Ecmpy-git)

## DLKcat method

The first model was achieved by following the DLKcat method.
However, automatic creation of the ecmpy-xml model leaded to a model with Nan values for the kcat and kcat-mw.
The model will NOT work with any Nan values. They had to be removed. Internally, results from a Brenda-search and
earlier DLKcat were available. These results were used to update the Nan values, and all DLKcat values were also
replaced by Brenda-values if available.

## Autopacman method

A second version of the model was created by following the autopacman method. This method relies more on databases
such as Brenda and Sakbio. The amount of empty kcat values is much lower for this model compared to Ecmpy
NB: No attempt was made yes to compare the kcat values from DLKcat and Autopacman.
Due to the higher level of completeness, the Autopacman model was selected as the model to use in further analyses.
Of course, both models can still be used and compared.

# How to load the model

NB: from this point onwards you can use an environment with Cobrapy, installation of Ecmpy is not required.
The relevant script for this part of the readme is: use_ecmpy_model.py

It is important to load the model via the custom-function, not via cobra.io.
Load the model via: get_enzyme_constraint_model

# How to use and export the model

Use the model as you would for the normal model. Make sure to set the right model medium.
Substrate uptake is defined via the reverse-reactions, e.g. EX_glc_e_reverse
Set the upperbound of these reactions to define the medium.
NB: Phosphate uptake is an exception! In the current version of the model, this does not have a reverse reaction,
as a result you need to set e.g. EX_pi_e.lower_bound = -0.15 to make phosphate available.

## Sampling

Next to FBA, we recommend to use sampling. We adopted a pFBA sampling method (adapted from the CFSA-github <>),
an example usage of the sampling is shown in the function: run_sampling.

# How to update kcats

It is possible to access the kcat and the kcat_MW values of each reaction via the reaction attribute.
For example: `enz_model.reactions.PDH_num1.kcat` and `enz_model.reactions.PDH_num1.kcat_MW`.

The MW-value (or total mass) is not stored in the model itself (but it is part of the Ecmpy output files).
The formula used to calculate each kcat_MW value in the model is this: `kcat_mw = <Kcat value (1/s)> * 3600 * 1000 / total_mass(MW)`

Updating only the kcat in the model will have NO effect to the simulation, it is crucial to ALWAYS also update the kcat_MW
value, as that is the one that is used in the cobra-constrained. You can change the kcat_MW value to any value you want,
and then rerun the simulation to see the effect.

# How to change the available proteome
The available proteome is represented by an upper_bound of the proteome-constraint in de model.
This number reflects the product of two variables:
The amount of protein per gram of dry weight (gr/gDW) and the fraction of protein that is an enzyme (gr/gr).
Together this number forms the fraction of enzymes per gDW in the cell.
We have used the same upper_bound as for the E. coli Ecmpy model: 0.227
You can update the proteome-bound by using the function: update_ub_proteome


# How to split the proteome in multiple pools
Use function: get_enzyme_constraint_sep_protein_pools
Example usage: try_out_separate_pools in use_ecmpy_model.py
Input required:
Tsv-file with the columns: Reaction and Pool
Any reaction can be linked to any defined pool. The size of each pool (including other) should be specified in a dict.
<NB: could later write this more elaborate, but this will do for now>


# To-do: improve the kcats

This following step was not managed yet, but it would be interesting to do:
Take the script: 03.ecModel_calibration.ipynb

And adapt the code to the Streptomyces model.
The direction of the while-loop statement needs to be changed to reflect the fact that the model is growing more
than expected (instead of less, as was the case for E. coli).
Also the required input for the argument `fluxes_infile_ori` and `EC_max_file` needs to be considered.

# NB: movement of files

Please note that I've moved files from my working directory to this github. At the moment of writing this read-me I did
not consider all of the relative paths in the scripts. This is something on my to-do lists. In the mean-time, if you
use any of the scripts, please consider that wrong relative paths to scripts or input files might give errors.
