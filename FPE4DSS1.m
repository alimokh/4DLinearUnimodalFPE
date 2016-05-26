%Copyright 2013 Rice University
%All Rights Reserved
%
%Created by Ali Arya Mokhtarzadeh
%Department of Mechanical Engineering
%

%This program solves the 4-Dimensional Linear System
clc;
clear all;
close all;

tic


FPE4DSS3SetupSym
%Create the Domain
%Create the Boundary driver eqs.
%Create the initial approximation/Boundary Enforcer eqs.
%Create the Exact Solution to compare SOMA results against
%Create the "exact solution", the residual error of the PDE
%Set Up the Loop
%Give the Optimization Information (number of basis functions, max error)
%Set Up Plots

%FPE1DSSnd2PlaceCenters

FPE4DSS5Loop
%set up a loop to run while the error is greater than the tolerance and the
%iteration count is less than the user specified value
%start by stepping to the next basis function (for a while loop)
%Run an optimization routine
%Build in a functionality in case you get a NaN value from the optimization
%routine, such as replacing the parameters you would have found in the
%routine with heuristic values, or by reducing the counter so the iteration
%repeats itself until it arrives at a solution without the NaN values
%Run an update to the approximation

%%Opt

%%Update