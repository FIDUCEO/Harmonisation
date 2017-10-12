#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author:         Arta Dilo and Peter Harris \ NPL MM
    Date created:   24-01-2017
    Last update:    20-03-2017
    Version:        10.0

Evaluate harmonisation uncertainty for the error structure in har. data 
via Monte Carlo (MC). An MC trial generates data from best estimates of 
harmonisation variables evaluated from regression on simulated data and errors 
generated by respecting the full error structure. Weighted ODR regression is 
performed on the generated data, fit coefficients are used to evaluate 
covariance 
evaluate uncertainty in 
uses full error structure to generate data in each MC trial. """