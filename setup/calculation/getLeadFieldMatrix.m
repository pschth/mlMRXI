function L = getLeadFieldMatrix( mrxisetup, H, RAM_threshold, verbose )
% calculates the system matrix of the linear forward model of a
% magnetorelaxometry imaging (MRXI) setup defined by 'mrxisetup' and
% magnetic field vectors 'H' in [V*s/(kg*m^2)] (if SI units are provided
% for the input).
% 
% Based on the equations from:
% Liebl, M., et al. "Quantitative imaging of magnetic nanoparticles
% by magnetorelaxometry with multiple excitation coils." Physics in
% Medicine & Biology 59.21 (2014): 6607. -> EQUATION (5) (or rather the
% split into \Delta B = L*X_MNP in EQUATION (9))
% !!ATTENTION!!: the product \chi(t_{mag})*\Delta \kappa_{1,2} equals 1
% [m^3/kg] in the computations of this function
% 
% INPUT:
% mrxisetup - MRXI setup configuration object
% 
% H - (nVoxels*nSequences) x 3 matrix with magnetic field vectors of
% each voxel and each excitation sequence. Call getMagneticField for
% convenient calculation of H.
% 
% RAM_threshold - scalar (optional); RAM-threshold in Gb. If the predicted
% RAM-usage of the calculation is below, the fast version of the algorithm
% is computed. Otherwise, the slower, RAM-efficient version is computed.
% Default: 4 [Gb]
% 
% verbose - boolean (optional); if false, no text messages are generated.
% Default: TRUE
%
% OUTPUT: L - (nSensors*nSequences) x nVoxels, system matrix of the linear
% MRXI forward model
    

    % pass default values if input is missing
    if nargin < 3 || isempty(RAM_threshold)
        RAM_threshold = 4;
    end
    
    if nargin < 4 || isempty(verbose)
        verbose = true;
    end
    
    %%% code
    % get some fixed setup characteristics
    nVoxels = mrxisetup.getNumberOfVoxels;
    nSensors = mrxisetup.getNumberOfSensors;
    nSequences = size(H, 1)/nVoxels;
    
    % reshape matrix of sensor orientations for later computations
    ns = permute(mrxisetup.sensorData.orientations, [1 3 2]);

    % define constants
    mu_0 = 4e-7*pi; % absolute permeability, [Vs/Am]

    % calculate all distance vectors between the sensors and the voxel
    % centers. This corresponds to all combinations of (r-r') in Liebl et
    % al., equation (5), and results in a 3 x nVoxels x nSensors matrix
    r = bsxfun(@minus, ...
            permute(mrxisetup.sensorData.centers, [1 3 2]),... 
            mrxisetup.ROIData.voxels);
    % calculate scalar distances of the above vectors
    r_rssq = rssq(r, 1);
    
    % calculated predicted RAM usage
    predict_RAM = 2* 3*nVoxels*nSensors*nSequences * 8 / 1024^3;

    if predict_RAM < RAM_threshold
        %% faster, RAM-inefficient version
        % print info message regarding which version is used
        if verbose
            fprintf('\tCalculating lead-field matrix - faster, RAM-demanding version. Predicted RAM usage = %.2fGb < RAM threshold = %.2fGb', predict_RAM, RAM_threshold);
        end
        
        % reshape H to a 4 dimensional array for efficient MATLAB matrix
        % calculations
        H = reshape(H', [3 nVoxels 1 nSequences]);

        % calculate r/||r||^5
        r_div_r_norm_pow_5 = bsxfun(@rdivide, r , r_rssq.^5); % r_div_r_norm_pow_5...3 x nVoxels x nSensors

        % calculate the system matrix L based on Liebl et al., equation (5)
        % (without the stuff right of H_{mag}(r'))
        L = sum( bsxfun( @times, ... % make scalar product of the stuff below
            (bsxfun( @times, 3*sum( bsxfun( @times, H, r_div_r_norm_pow_5 ), 1), r)... corrsponds to (3*dot(r,H)*r)/||r||^5
            -bsxfun( @rdivide, H, r_rssq.^3)),... corresponds to -H/||r||^3
            ns), 1)*(mu_0/(4*pi)); % corresponds to *ns*mu_0/(4*pi)
        
        % reshape L to a (nSensors*nSequences) x nVoxels matrix
        L = reshape(L,nVoxels,[])';

    else
        %% slower, RAM-efficient version
        % print info message regarding which version is used
        if verbose
            fprintf('\n\tCalculating lead-field matrix - slower, RAM-efficient version. Predicted RAM usage = %.2fGb > RAM threshold = %.2fGb', predict_RAM, RAM_threshold);
        end
        
        % reshape H to a 3 dimensional array for easier indexing of the
        % activation sequences
        H = reshape(H', [3 nVoxels nSequences]);
        
        % initialize system matrix
        L = zeros(nSensors*nSequences, nVoxels);

        % print status message
        if verbose
            fprintf('\n\tCalculating sequences: '); 
            delCount = 0;
        end
        
        % for all sequences concatenate the sub-system matrices of the
        % individual activation sequences
        for s = 1:nSequences
            
            % print progress
            if verbose
                fprintf(1, repmat('\b',1,delCount)); % delete old progress
                delCount = fprintf('%.0f / %.0f', s, nSequences); % print new progress
                pause(1e-12); % pause for printing
            end
        
            % left term of Liebl et al., equation (5):
            % (3*dot(ns,r)*dot(r,H))/||r||^5
            lt = (3*sum(bsxfun(@times, ns, r), 1) .* sum(bsxfun(@times, r, H(:,:,s)),1))./...
                (r_rssq.^5);

            % left term of Liebl et al., equation (5): 
            % dot(ns,H)/||r||^3
            rt = sum(bsxfun(@times, ns, H(:,:,s)),1)./(r_rssq.^3);

            % get system matrix fill index for current activation sequence
            idx = ((s-1)*nSensors+1) : (s*nSensors);

            % fill current system matrix block
            L(idx,:) = permute( mu_0/(4*pi)* (lt - rt) , [3 2 1]);
        end
    end
    % done
    if verbose
        fprintf('\n\tDone: Computation lead-field matrix.\n');
    end
end
