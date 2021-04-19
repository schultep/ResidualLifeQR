# ResidualLifeQR
Residual lifetime quantile methods for personalized screening intervals incorporating longitudinal biomarker data

The MATLAB code provided replicates, in part, Tables 1 and 2

Requires MATLAB 9.5 and Statistics and Machine Learning Toolbox 11.4.  Earlier versions may work, but have not been tested.



Included programs and purpose:

    RunSim.m - PRIMARY program. This program will run the simulation study, calling other programs as applicable

    simDat2.m - Generates simulated dataset based on parameters inputed

    QRdata.m - Creates a class object which organizes the data for QRfit_sim1.m and MLEfit_sim1.m

    MLEfit_sim1.m - Generates starting value for quantile regression. The code contins ad-hoc MLE fit

    QRfit_sim1.m - The quantile regression portion of the code


