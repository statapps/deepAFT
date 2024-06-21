
# DeepAFT: AFT model for survival data with deep learning neural networks.

"deepAFT" is a R package for AFT model with deep learning.

### Update June 21, 2024
The "deepAFT" software is now part of the R package "dnn" available from the Comprehensive R Archive Network (https://CRAN.R-project.org/package=dnn) and can be installed using R commands below:

  install.packages("dnn")

  library(dnn)

  help("deepAFT")

References: 

1. Norman, P. and Chen, B. E. (2018). DeepAFAT: A nonparametric accelerated failure time model with artifical neural network. M. Sc. thesis report. Department of Public Health Sciences, Queen's University, Canada.

2. Norman, P. Li, W., Jiang, W. and Chen, B. E. (2024). deepAFT: A nonlinear accelerated failure time model with artificial neural network. Statistics in Medicine. Published June 18, 2024. <https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10152>. 


### The following steps are not longer required.
Please use the steps below to install 'deepAFT' package:

1. First, you need to install the 'devtools' package. You can skip this step if you have 'devtools' inst
alled in your R. Invoke R and then type

  install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "deepAFT" package with R commond

  install_github("statapps/deepAFT")
 
4. Install "keras" and "tensorflow" R package.

   install.packages("keral")
   install.packages("tensorflow")
   
5. Download and install "Python"

6. Install tensorFlow for your compute system (run the following code under commondline of your system (Windows, Linix, or OS X))

   pip install tensorflow
