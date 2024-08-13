# Uncertainty estimation
These are the scripts to reproduce the results in the manuscript *"Uncertainty in cardiovascular digital twins despite non-normal errors in 4D flow MRI: identifying reliable biomarkers such as ventricular relaxation rate"* by Tunedal K, Ebbers T, and Cedersund C.

If you use this implementation in your academic projects, please cite this paper.



This project introduces a method for uncertainty estimation despite non-normal errors, applied to a cardiovascular model. The method is validated by creating 100 bootstraps of simulated 4D flow MRI and cuff pressure data with added sources of errors and comparing the resulting model uncertainty with the true values. The method is also applied to a real clinical dataset.

The method defined in the manuscript consists of three steps: 
1. Estimate the combined data uncertainty,
2. Estimate the optimal model-derived biomarkers by fitting to data, and
3. Estimate the biomarker uncertainty with profile likelihood

Descriptions on how to perform the three steps is found below.

For more information about the cardiovascular model and the clinical dataset, see Tunedal et al 2023 ([10.1113/JP284652](https://doi.org/10.1113/jp284652)). The cardiovascular model file with all model equations is located in Modelfiles/avatar_corr.m.

# Re-create result figures
Re-create the main results: Run the script *createFigures.m*. This will take a while, since it takes time to plot all biomarkers and data sets and save the figures in high resolution. 

# Applying the method to real data: example using the clinical dataset from Tunedal et al 2023 to find reliable biomarkers
The selected subjects from Tunedal *et al* are saved in the files `Data/dataP*.mat`. The structure of those datasets are examples of the structure needed for the scripts to be used on other data.

Below, a short guide to the scripts for the different steps in the method is presented. Scripts to add combined data uncertainty (Step 1) are found in the *Data* folder. All scripts for profile likelihood and parameter estimation (Step 2-3) are found in the folder *Optimization*.

## Step 1 in the method: add combined data uncertainty
An example on how to add estimated data uncertainty distributions (Step 1 in the metod) to new data is shown in *addErrosToRealData.m* (in the *Data* folder), where the distributions are loaded from the file *dataSimulated.mat*, and the RR error is estimated using the empirical function *calcRRsmoothingerror.m*.

## Step 2 in the method: find optimal biomarker values with parameter estimation
To re-do the optimization to find the best fit to data, go to the Optimizaion folder and run runParamEstimationSeveral_ESS_asNSC. 
The optimization ran at the Swedish NSC. The script runParamEstimationSeveral_ESS_asNSC runs the optimization in a similar same way as it was done at the NSC, but settings might need to be modified and the optimization might need to be re-run several times. On a normal workstation, this could take days or weeks.

The results from the parameter estimation are saved in the *Parameters/ESS* folder for each subject/dataset.

## Step 3 in the method: profile likelihood
To run profile likelihood (PL), the script *Optimization/runPL_asNSC* can be used. Depending on if the PL is applied to replicate the results based on the bootstrapped data or if the method is applied to new data or to the clinical subjects, different PL scripts are used. This is because for the bootstraps, the PL steps were pre-defined from the start, while the classical PL starts from the optimal parameter set found in Step 2. To run the classical PL, use the *EstimatePL_classic.m* function. It needs to be run iteratively or for a long time to make sure that the PL succeeds. For optimal speed, is recommended to use a supercomputing center or a workstation with a lot of CPU and run the estimations in parallell.


Which biomarkers to be estimated is defined either as input to the estimation functions, or, if several biomarkers are to be estiamted at the same time, inside each function.


The results from the profile likelihood are saved in the *Parameters/Pl* folder for each subject/dataset.


To combine all separate folders for each biomarker, run the script *collectAllPLparams.m*.


To plot the results from the profile likelihood, the script *plotPL_realdata.m* can be used (with modification to select the wanted biomarkers and datasets).

# Evaluation of the uncertainty estimation method: The 100 sampled data and defining data uncertainty distributions
The bootstrapped estimation data is saved in *Data/SampledErrors_simRandomSystematic*. The "true" simulated dataset is saved in *Data/dataSimulated.mat*, and can be re-created by running a model simulation. The script used to create the true dasaset is *createTrueData_simulated*, where the error distributions are defined and a model simulation is done using the parameters saved in *trueParameters.mat*.

To create the bootstrap data, measurement errors were sampled using *createSampledMeasurementError.m*, which calls the sampling function *sampleData* several times, and draw errors from the data uncertainty distribution. Each resulting error e is saved in the folder *Data/SampledErrors_simRandomSystematic* as *sampledMeasurementErrorn.mat* where *n* is the bootstrap number 1-100.

To load the bootstraps, the function *Optimization/loadSampledMeasurementError.m* is used. The function loads the *trueData* and the sampled error *e* of bootstrap n, and finally adds the errors using the *sampleData* function and sending in the *trueData* nd the loaded error *e*.

## Smooting error due to RR variation in 4D flow MRI
The empirical function to replicate the RR error is found in *Data/calcRRsmoothingerror_fixed.m*. The script requires input data of the same format as the estimation data. An example of such a dataset can for example be found in *dataP1.mat* which is one of the datasets from Tunedal et al 2023.

The simulations to estimate the size of the RR variation error during bootstraping were done with the script *Optimization/sampleData.m*, calling the simulation function *simulate_avatar_RR*.

# Requirements
The model is created in the [AMICI toolbox](https://doi.org/10.1093/bioinformatics/btab227) in MATLAB. To compile the model, MATLAB 2017b or earlier is needed, but to run the already compiled model any later matlab verison works. The provided compiled model is compiled on Windows, but will not work on Linux or macOS. To re-compile: run GenerateModels in the folder Modelfiles. To compile the model, you need a valid C-compiler (such as xcode on Mac or MinGW on Windows. Run mex -setup to check if you have an installed compiler in matlab) and the MATLAB Symbolic Math Toolbox.

To run the other scripts, the following MATLAB toolboxes are needed: Statistics and Machine Learning Toolbox, Optimization Toolbox, Parallel Computing Toolbox, Symbolic Math Toolbox, Signal Processing Toolbox. The code was created with Matlab R2021a and R2023a. Earlier matlab versions might not be compatible with many of the scripts.

Additionally, the following toolboxes are included in the folder Requirements and are needed to reproduce the results:
* For model simulation, the [AMICI toolbox](https://doi.org/10.1093/bioinformatics/btab227) is needed. 
* For parameter estimation with eSS, the [MEIGO toolbox](https://doi.org/10.1186/1471-2105-15-136) is needed.


## Author
Kajsa Tunedal (kajsa.tunedal@liu.se)

## License
The MIT License (MIT)

Copyright (c) 2024 Kajsa Tunedal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

