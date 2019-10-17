% function create_fig_3A_top
%
% This function creates the top row of Figure 3A of the paper:
%
%   *********************************************************************
%   * Van den Berg, R. & Ma, W.J. (2018). A resource-rational theory of *
%   *   set size effects in human visual working memory. Elife.         *
%   *********************************************************************
%
% Written by Ronald van den Berg, 2018

function create_fig_3A_top
% Fit model and create plot for worst-fitting subject in terms of R^2,
% which is subject 10 in experiment 6:
fit_model(6,10)  