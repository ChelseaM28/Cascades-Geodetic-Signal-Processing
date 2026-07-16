### covariance_realism.py
# This script will quantify the uncertainty of my velocity estimates using a distribution comparison 
# of Monte-Carlo samples of uncertainty realism metrics against the matching chi-squared distribution. 


### Motivating publications:
# Working Group on Covariance Realism. (n.d.). Covariance and Uncertainty 
#        Realism in Space Surveillance and Tracking (A. B. Poore, J. M. Aristoff, 
#        & J. T. Horwood, Eds.). Air Force Space Command Astrodynamics 
#        Innovation Committee.
# Zaidi, Waqar H. , and Matthew D. Hejduk. Earth Observing System Covariance Realism. American Institute 
#        of Aeronautics and Astronautics, 1 Mar. 2016.
#
# Please NOTE that Poore et al quantifies uncertainty against REAL truth, while my uncertainty can
# only be validated against against *simulated* truth. Claiming true uncertainty realism in the context
# of orbit determination, I'd be making an overstatement. However, in the context of GNSS geodesy, 
# my data is adequately 'realistic' since I'm not assuming white noise. TODO: add source.







