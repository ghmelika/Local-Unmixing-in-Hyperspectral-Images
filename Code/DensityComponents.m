function DensityComponents()

% this code implements equations for inference/learning
% in the following variance component model:

% p(u|B,v) : distribution over the linear features (ICA basis funcs, etc)
%            Gaussian with variance \sigma^2 = exp(B*v)
% p(v)     : Laplacian, variance = 1

% notes:
%   - using MAP estim. during learning favors B growing w/o bounds
%       weight decay on B (or fixed norm as here) can fix this

% yan karklin. jan 2009.


I = 20;  % dimensionality of the data (e.g. image patches)
J = 4;  % number of higher-order basis functions
N = 500;  % number of data samples in a single batch

trueB = setInitB(I,J);

epsv = 0.02;  % step size for v (might want to adjust so it tapers)
epsB = 0.005; % step size for B

% initialize B
B = randn(I,J);
% fix norm of B to the norm of trueB (need to deal with degeneracy of B
%                                     growing without bound)
B = B * diag(sqrt(sum(trueB.^2)./sum(B.^2)));

L = [];
clf; drawDisplay(B, trueB,[],[]);

for i=1:100,

  % get new batch of data here (e.g. resample image patches)
  u = getData(trueB, N);  % u should be a [I x N] matrix
  
  % initialize MAP estimates \hat{v}  
  v = .1*randn(J,N);

  % use gradient ascent to obtain MAP estimate \hat{v}
  Linf = [];
  for j=1:100,
    % d log p(u|B,v)/dv + d log(p(v)/dv
    dv = .5 * B'*(u.^2./ exp(B*v) - 1) - sqrt(2)*sign(v);
    v = v + epsv * dv;
    
    Linf = [Linf  calcL(B, u, v)];
  
  end;

  % update B
  dB = .5 * (u.^2 ./ exp(B*v) - 1) * v';

  B = B + epsB * dB;
  % adjust norm
  B = B * diag(sqrt(sum(trueB.^2)./sum(B.^2)));

  L = [L calcL(B,u,v)];

  drawDisplay(B,trueB,L,Linf);


end;

function drawDisplay(B, trueB, L, Linf)
[I J] = size(B);

subplot(2,2,1); plot(trueB + 2*repmat([1:J],[I 1]));  title('true B');
subplot(2,2,2); plot(    B + 2*repmat([1:J],[I 1]));  title('learned B');
if ~isempty(Linf), subplot(2,2,3); plot(Linf); title('inference  log p(v|u,B)'); end;
if ~isempty(L),    subplot(2,2,4); plot(L); title('param learning  log p(u|B)'); end;
drawnow;



function B = setInitB(I,J)
% make cosine basis set
B =[];
for j=1:J,
  B = [B .5*cos(2*pi*(1:I)/I*j+2*pi*rand)'];
end;


function u = getData(trueB, N);
[I J] = size(trueB);

% Laplacian v
v    = exprnd(1/sqrt(2)*ones(J,N)) .* sign(randn(J,N));
varu = exp(trueB * v);
u    = sqrt(varu) .* randn(I,N);




function L = calcL(B, u, v)
% computes the log-likelihood
% logP(u|B) = logP(u|B,v) + logP(v)

[I N] = size(u);
J     = size(v,1);

Bv  = B * v;

% conditionally gaussian prior P(u|B,v),
% variance  sigma^2 = exp(Bv);  
% log p(u|sigma^2) = -.5*log(2*pi*sigma^2) -.5 * u^2 / sigma^2
lPu = -.5*I*log(2*pi) - .5*sum(sum(Bv + u.^2 ./ exp(Bv)));

% assume Laplacian prior on v
lPv = -.5*J*log(2) - sqrt(2)*sum(sum(abs(v)));

L = lPu/N + lPv/N;