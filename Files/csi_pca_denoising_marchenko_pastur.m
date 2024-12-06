function denoised_patch = csi_pca_denoising_marchenko_pastur(inpMat, svd_method)
% This method uses principal component analysis to remove noise components
% from the signal, which can be automatically be identified using the 
% Marchenco-Pastur distribution.
%
% i.e. denoise by removing eigenvalues equal or below the noise
% eigenvalues.
%
% The returned denoised_patch should be iterativly summated and divided by 
% the numel of the local area to increase SNR.
%
% Optional: svd_method, default Matlab svd-function [0] or fast but lower
% accuracy svdecon-function [1]
%
% Method:
% Calculate eigenvalues for the input matrix, which is signal and noise
% from a local area in a spectrum-matrix. Eigenvalues from solely noise 
% should represent a Marchenko-Pastur distribution. Values outside are 
% components which carry signal information. Therefor, these eigenvalues 
% below the threshold can be set to zero, removing noise-components.
%
% Q. van Houtum, 2023;
% quincyvanhoutum@gmail.com

% NFO Principle Component Analysis:
% PC have a particular ordering – each principal component points in the 
% direction of maximal variance that is orthogonal to each of the previous 
% PC. In this way, each PC accounts for the maximal possible amount of 
% variance, ignoring the variance already accounted for by the previous PC. 
% These directions turn out to be precisely the eigenvectors of  AA'|A'A.
% The singular vector corresponding to the leading singular value usually 
% can be thought as the principal direction in which the data is organized.
% Calculating the SVD consists of finding the eigenvalues and eigenvectors 
% of AAT and ATA. The eigenvectors of A'A make up the columns of V , the 
% eigenvectors of AA'  make up the columns of U. Also, the singular values 
% in S are square roots of eigenvalues from AA'|A'A.
%
% Veraart et al; 2016; 10.1016/j.neuroimage.2016.08.016
% The variance ˜σ2 of the residual noise,
% contained within the P significant components, can be computed as:
% ˜σ2 = σ2 − Pσ
% In particular, we estimate the number ˆP of significant components by 
% incrementally increasing p until
% ∑(over M, incrementing p) = λi ≥ (M − p) ˆσ2(p)

if nargin == 1, svd_method = 0; end

% Dimensions
n = size(inpMat,1); % row
m = size(inpMat,2); % col
r = min(m,n);       % Smallest == #voxels

% Singluar Value Decomposition
% Economy to speed up by removing values from calculations. See nfo of svd.
% S = square roots of eigen values AA'|A'A
if svd_method == 0    
    [U, S, V] = svd(inpMat,'econ'); 
elseif svd_method == 1
    [U, S, V] = svdecon(inpMat); 
end

% Eigenvalues: |A-lambda| = 0, solve for lambda
% S-squared = eigenvalues of A'A or AA'.
eigenvec = flip(diag(S).^2);  % Principle values
% Lambda i.e. eigenvalue of inpMat; solve for n.
lam_r = eigenvec(1) / n; % Principle components

clam = 0;
tolerance = 1;
sigma = NaN;

% Increment over the smallest dimension of input matrix M.
for pp = 1:r
    lam = eigenvec(pp) / n; % Eigenvalue lambda for pp (increment)
    clam = clam + lam; 
    gam = (m-pp) / n; % γ = ˜M/N with ˜M = M − P; see veraart et al.
    varA = clam / pp; % Variance after pp components
    % With λ± = σ2(1 ± √γ)2, from the relation λ+ − λ− = 4 √(γ)σ2
    varB = (lam - lam_r) / (4*sqrt(gam));
    if varB < varA
        % The threshold for setting the MP-distribution is found here.
        % Thus inside the MP-distribution set to 0 after the increment
        % loop.
        sigma = sqrt(varA + varB) / 2; % For noise-map
        threshold = pp - tolerance;
    end
end
% Correct threshold to minimum matrix size.
threshold = r - threshold;
% Undo flip of eigenvec.
eigenvec = flip(eigenvec);


% Modified to avoid local areas without denoising
if threshold <= 1, threshold = 1; end
Snew = zeros(size(S));
Snew(1:threshold, 1:threshold) = diag(sqrt(eigenvec(1:threshold)));

% Calculate input matrix using the new singular value matrix S.
% i.e. A = USV' with S = S-denoised via Marchenco-Pastur distribution.
denoised_patch = reshape(U*Snew*V',size(inpMat));

