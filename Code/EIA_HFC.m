% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // PaÃ­s Vasco (UPV/EHU), Spain, released under the terms of the GNU 
% // General Public License.
% //
% // Endmember Induction Algorithms Toolbox is free software: you can redistribute 
% // it and/or modify it under the terms of the GNU General Public License 
% // as published by the Free Software Foundation, either version 3 of the 
% // License, or (at your option) any later version.
% //
% // Endmember Induction Algorithms Toolbox is distributed in the hope that it will
% // be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% // of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% // General Public License for more details.
% //
% // You should have received a copy of the GNU General Public License
% // along with Endmember Induction Algorithms Toolbox. 
% // If not, see <http://www.gnu.org/licenses/>.
% // ====================================================================

function [vd] = EIA_HFC(data,alfa)
%% vd = EIA_HFC(data,alfa)
% 
% Manuel Grana <manuel.grana[AT]ehu.es>
% Miguel Angel Veganzones <miguelangel.veganzones[AT]ehu.es>
% Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
% Euskal Herriko Unibertsitatea (UPV/EHU)
% http://www.ehu.es/computationalintelligence
% 
% Copyright (2011) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
%
% Virtual dimensionality by HFC method
% ------------------------------------------------------------------------------
% Input:   data      : column data matrix [bands x pixels]
%          alfa      : vector of false alarm probabilities [1 x p] (default: [10^(-3) 10^(-4) 10^(-5)])
%
% Output:  vd        : vector of virtual dimensionality values [1 x p]
%
% Bibliographical references:
% [1] Chang, C.-I. and Du, Q., â€œEstimation of number of spectrally distinct signal sources in hyperspectral imagery,â€? Geoscience and Remote Sensing, IEEE Transactions on,  vol. 42, 2004, pp. 608-619.
% [2] Wang, J. and Chang, C.-I., â€œApplications of Independent Component Analysis in Endmember Extraction and Abundance Quantification for Hyperspectral Imagery,â€? Geoscience and Remote Sensing, IEEE Transactions on,  vol. 44, 2006, pp. 2601-2616.
% [3] J. Wang and Chein-I Chang, â€œIndependent component analysis-based dimensionality reduction with applications in hyperspectral image analysis,â€? Geoscience and Remote Sensing, IEEE Transactions on,  vol. 44, 2006, pp. 1586-1600.

%% Parameters
if (nargin < 1)
    error('Insufficient parameters');
end
if (nargin < 2)
    alfa = [10^(-3) 10^(-4) 10^(-5)];
end

%% data size
[nvariables, nsamples] = size(data);

%% Correlation and covariance matrix
%R = corr(data); -> ERROR
R = (data*data')/nsamples; % !!
K = cov(data');

%% Eigenvalues
lcorr = sort(eig(R),'descend');
lcov = sort(eig(K),'descend');

%% Differences and variances
diff = lcorr - lcov;
variance = sqrt(2*(power(lcorr,2)+power(lcov,2))/nsamples);

%% Hypothesis Test %%%% nemifahmam
p = size(alfa,2);
vd = zeros(1,p);
for i=1:p
    tau = -norminv(alfa(i),zeros(nvariables,1),variance);
    vd(i) = sum(diff > tau);
end