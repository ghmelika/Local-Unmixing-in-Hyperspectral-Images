% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // País Vasco (UPV/EHU), Spain, released under the terms of the GNU 
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

function [E,C] = EIA_FIPPI(data,p,maxit)
%% [E,C] = EIA_FIPPI(data,p,maxit)
% 
% Manuel Grana <manuel.grana[AT]ehu.es>
% Miguel Angel Veganzones <miguelangel.veganzones[AT]ehu.es>
% Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
% Euskal Herriko Unibertsitatea (UPV/EHU)
% http://www.ehu.es/computationalintelligence
% 
% Copyright (2011) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
%
% Fast Iterative Pixel Purity Index (FIPPI) endmembers induction algorithm.
% ------------------------------------------------------------------------------
% Input:   data      : column data matrix [nvariables x nsamples]
%          p         : number of endmembers to be induced. If not provided it is calculated by HFC method with tol=10^(-5)
%          maxit     : maximum number of iterations. Default = 3*p
%
% Output:  E         : set of induced endmembers [nvariables x p]
%          C         : induced endmembers indexes vector [nsamples] ...
%                    ... with {0,1} values, where '1' indicates that ...
%                    ... the corresponding sample has been identified ...
%                    ... as an endmember.
%
% Bibliographical references:
% [1] Chang, C.-I., “A fast iterative algorithm for implementation of pixel purity index�?, Geoscience and Remote Sensing Letters, IEEE, vol. 3, nº. 1, págs. 63-67, 2006.

%% Parameters
if (nargin < 1)
    error('Insufficient parameters');
end
if (nargin < 2 || p <= 0)
    p = EIA_HFC(data,10^(-5));
end
if nargin < 3 || maxit < 0
    maxit = 0;
end

%% data size
[nvariables,nsamples] = size(data);

%% Dimensionality reduction by PCA

[pc, zscores, pcvars] = pca(data);%[pc, zscores, pcvars] = princomp(data);
data_pca = squeeze(zscores(1:p,:)); % removing singleton dimensions

%% Initialization
C = zeros(1,nsamples);
E = [];
% Initial skewers
skewers = EIA_ATGP(data_pca,p);

%% Algorithm
stop = false; % stop condition
it = 1; % iterations
ne = 0; % number of endmembers
idx = []; %indexes of the induced endmembers
while ~stop
    % Calculate Nppi
    Nppi = zeros(1,nsamples);
    proj = data_pca'*skewers;
    [C1,I1] = min(proj);
    [C2,I2] = max(proj);
    for j=1:p
        Nppi(I1(j)) = Nppi(I1(j)) + 1;
        Nppi(I2(j)) = Nppi(I2(j)) + 1;
    end
    % Check new skewers
    r = find(Nppi);
    [skewers_r, i_r, i_sk] = union(data_pca(:,r)',skewers','rows'); 
    % i_r & i_sk : index vectors
    % skewers_r : the sorted combination of the rows data_pca(i_r,:) and skewers'(i_sk,:)
    dif = size(i_r,1);
    if dif == 0
        stop = true;
	idx = r;
    else
        % new skewers
        skewers = skewers_r';
        % Check iterations
        if maxit > 0 && it == maxit
            stop = true;
            idx = r;
        else
            it = it + 1;
        end
    end
end

% Endmembers
for j=1:size(i_r)
    ne = ne + 1;
    C(r(j)) = 1;
    E(:,ne) = data(:,r(j));
end

