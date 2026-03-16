classdef RECON_TIKH_ITER_REAL < handle
    % from:
    % Alessandro Buccini, Marco Donatelli, and Lothar Reichel: Iterated Tikhonov regularization with a general penalty term
    % with original source:
    % Engl HW, Hanke M, Neubauer A.Regularization of inverse problems, vol. 375. Springer Science & Business Media,1996.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC MEMBERS
    properties (Access = public)
        
        req; % required algorithm parameters
        aux; % non-mandatory auxiliary algorithm parameters
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE MEMBERS
    properties (Access = private)
        
        int; % internal algorithm parameters
        info; % informational algorithm parameters
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC FUNCTIONS
    methods (Access = public)
        function o = RECON_TIKH_ITER_REAL(L,b)
            % constructor
            %
            % INPUT:
            % SETUP - MRXI setup configuration object
            %
            % L - nActivations*nSensors x nVoxels matrix; leadfield matrix
            %
            % b - nActivations*nSensor vector; measurement vector
            %
            % 
            
            o.info.name       = 'TIKH_ITER_REAL';
            o.info.className       = ['RECON_' o.info.name];
            
            if ~size(L,1)==length(b)
                error('Matrix dimension and measurement dimension do not match.');
            end
            
            o.setSystem(L,b);
            o = o.init();
                        
        end
        
        function o = copySettings(o, RECON, forceCopy)
            % copies the required and auxiliary settings of an object of
            % the same RECON class
            %
            % INPUT:
            % RECON - object; RECON class object
            %
            % forceCopy (optional) - bool; forces the copying of settings although the
            % wrong RECON class is provided; default: false
            
            if nargin < 2
                forceCopy = false;
            end
            
            if isa(RECON, o.info.className) || forceCopy
                o.req = RECON.req;
                o.aux = RECON.aux;
            else
                error('Wrong RECON class.');
            end
            
        end
        
        function o = setSystem(o,L,b)
            % sets leadfield matrix and measurement vector. Use e.g. after
            % copying the object for use of different system with same
            % parameter settings.
            %
            % INPUT:
            % L (optional) - nActivations*nSensors x nVoxels matrix; leadfield matrix
            %
            % b (optional) - nActivations*nSensor vector; measurement vector
            
            if ~isempty(L)
                o.int.L = L;
            end
            if ~isempty(b)
                o.int.b = b;
            end
        end
        
        function printExplanation(o)
            
            % get all fieldnames
            namesReq = fieldnames(o.req);
            namesAux = fieldnames(o.aux);
            namesInfo= fieldnames(o.info);
            
            % get the longest name
            nLongestName = 0;
            for i = 1:length(namesInfo)
                if length(namesInfo{i})>nLongestName
                    nLongestName = length(namesInfo{i});
                end                
            end
            
            fprintf(['\nAlgorithm: ' o.info.name]);
            % print explanations of required variables
            fprintf('\n\nRequired variables:');
            for r = 1:length(namesReq)
                for i = 1:length(namesInfo)
                    if strcmp(namesReq{r}, namesInfo{i})
                        nFill = nLongestName - length(namesReq{r});
                        fillStr = repmat(' ', 1, nFill);
                        fprintf(['\n' namesReq{r} fillStr '\t' o.info.(namesReq{r})]);
                        continue;
                    end
                end
            end
            
            fprintf('\n\nAuxiliary variables:');
            % print explanations of auxiliary variables
            for a = 1:length(namesAux)
                for i = 1:length(namesInfo)
                    if strcmp(namesAux{a}, namesInfo{i})
                        nFill = nLongestName - length(namesAux{a});
                        fillStr = repmat(' ', 1, nFill);
                        fprintf(['\n' namesAux{a} fillStr '\t' o.info.(namesAux{a})]);
                        continue;
                    end
                end
            end
            fprintf('\n');
        end
        
        function int = getInternals(o)
            % returns the internal parameters of the object
            int = o.int;
        end
        
        function o = start(o)
            % starts the reconstruction algorithm
            
            o.int.tic = tic;
            o.reconstruct();
        end
        
        function [result, corrstorage] = getResult(o)
            % returns the result
            
            result = o.int.result;
            if nargout > 1
                corrstorage = o.int.corrstorage;
            end
        end
        
 
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS
    methods (Access = private)
        
        function o = init(o)
            % initializes the reconstruction algorithm.
            
            % set all known internal parameters
            o.initInternals();
            % set informational parameters
            o.initInfo();

            %%% required
            o.req.lambda                = 1e-5;                                 % Tikhonov regularization parameter
            o.req.Gamma                 = ones( o.int.nVox, 1);                 % Tikhonov weighting vector
            o.req.stepsize              = 1e0;                                  % Gradient descent stepsize
            o.req.nonNeg                = true;                                 % Non-negativity requirement
            o.req.nIter                 = 1e4;                                  % Number of iterations
            o.req.deltaMin              = 0;                                    % in %, minimum normalized RMS deviation from one iteration step to the next, otherwise abort
            
            %%% auxiliary
            o.aux.xInit                 = zeros( o.int.nVox, 1);    % Estimation of the result (zeros = cold start, estimation = warm start)
            o.aux.xOrig                 = zeros( o.int.nVox, 1);    % Ground truth for correlation control during computation
            o.aux.OutputInterval        = 100;                      % Defines the verbose interval
            o.aux.verbose               = true;
            o.aux.neglectZeroEntryLimit = Inf;                      % Removes variables from the equations, if they return zero neglectZeroEntryLimit times in a row
        end
        
        function o = initInternals(o)
            % defines all required internal parameters known after calling
            % the constructor
            
            % get setup constants
            o.int.nVox = size(o.int.L,2);
            o.int.rnVox = sqrt(o.int.nVox);
            
            % set algorithm parameters
            o.int.i = 0; % algorithm iteration count
            o.int.vi = 0; % verbose count
            
            % initialize deviation variable
            o.int.dev = inf;
        end  
        
        function o = initInfo(o)
            % defines all informational parameters known after calling
            % the constructor
            
            o.info.runtime = 0; % algorithm runtime
            
            %%% informational (field names must be identical to required
            %%% field names)
            o.info.lambda = 'Scalar - Tikhonov regularization parameter';
            o.info.Gamma = 'Vector(nVoxel) - Tikhonov weighting vector';
            o.info.nonNeg = 'Bool - Non-negativity requirement';
            o.info.nIter = 'Scalar - Number of iterations';
            o.info.deltaMin = 'Scalar - Minimum deviation of one iteration step to the next. Otherwise abort.';
            %%% informational (field names must be identical to auxiliary
            %%% field names)
            o.info.xInit = 'nVoxel vector - Initialisation of the result (zeros = cold start, estimation = warm start)';
            o.info.xOrig = 'nVoxel vector - Ground truth for correlation comparison during computation (zeros = no comparison, otherwise = comparison)';
            o.info.OutputInterval = 'Scalar - Defines the verbose interval';
        end
        
        function o = reconstruct(o)
            % this is the actual reconstruction algorithm.
            
            % initialize zero entry removal variables
            if ~isinf(o.aux.neglectZeroEntryLimit)
                o.int.ZER = true;
                o.int.ZERcounter = zeros(o.int.nVox,1);
            else
                o.int.ZER = false;
            end
            o.int.removeZE = false(o.int.nVox,1);
            
            % if algorithm iteration count = 0, the reconstruction
            % algorithm has not been initialized
            if o.int.i == 0
                % initialize result
                o.int.result = o.aux.xInit;
                
                % reshape weighting vector if necessary
                [o.int.Gamma_s1, o.int.Gamma_s2] = size(o.req.Gamma);
                if o.int.Gamma_s1 == 1 || o.int.Gamma_s2 == 1
                    o.req.Gamma = diag(o.req.Gamma);
                end
                if ~all(size(o.req.Gamma) == o.int.nVox)
                    error('Gamma has the wrong size');
                end
                
                
                % project to non-negativity if necessary
                if o.req.nonNeg
                    o.int.result(o.int.result<0) = 0;
                end
                
                % perform zero entry removal counting
                if o.int.ZER
                    % get all zero entries
                    o.int.zero = o.int.result == 0;
                    % set all zero entries that were already marked to true
                    o.int.zero(o.int.removeZE) = true;
                    % reset all nonzero entries in the counter
                    o.int.ZERcounter(~o.int.zero) = 0;
                    % increment all zero entries in the counter
                    o.int.ZERcounter = o.int.ZERcounter + o.int.zero;
                    % mark zero entries if threshold is hit
                    o.int.removeZE = o.int.removeZE | o.int.ZERcounter == o.aux.neglectZeroEntryLimit;
                end

                % initialize calculation constants
                o.int.LtL = o.int.L'*o.int.L;
                o.int.invLtLw = inv(o.int.LtL + o.req.lambda*o.req.Gamma'*o.req.Gamma);
                o.int.Ltb = o.int.L'*o.int.b;
                o.int.Const1 = o.int.invLtLw*o.int.Ltb;
                o.int.Const2 = eye(o.int.nVox) - o.int.invLtLw*o.int.LtL;
            end
            
            if o.req.deltaMin == 0
                o.int.dev = Inf;
                while o.int.i < o.req.nIter && o.req.deltaMin < o.int.dev
                    o.int.preresult = o.int.result;
                    o.int.keepIdx = ~o.int.removeZE;
                    % perform tikhonov step
                    o.int.result(o.int.keepIdx) = o.int.Const1(o.int.keepIdx) + ...
                        o.int.Const2(o.int.keepIdx, o.int.keepIdx)*o.int.result(o.int.keepIdx);

                    % nonnegativity step
                    if o.req.nonNeg %&& o.int.i < (o.req.nIter-1)
                        o.int.result(o.int.result<0) = 0;
                    end
                    
                    % perform zero entry removal counting
                    if o.int.ZER
                        % get all zero entries
                        o.int.zero(~o.int.removeZE) = o.int.result(~o.int.removeZE) == 0;
                        % reset all nonzero entries in the counter
                        o.int.ZERcounter(~o.int.zero) = 0;
                        % increment all zero entries in the counter
                        o.int.ZERcounter(~o.int.removeZE) = o.int.ZERcounter(~o.int.removeZE) + o.int.zero(~o.int.removeZE);
                        % mark zero entries if threshold is hit
                        o.int.removeZE(~o.int.removeZE) = o.int.removeZE(~o.int.removeZE) | ...
                            o.int.ZERcounter(~o.int.removeZE) == o.aux.neglectZeroEntryLimit;
                    end

                    % increase iteration count
                    o.int.i = o.int.i + 1;

                    % print algorithm progress
                    o.printProgress();
                end
            else
                while o.int.i < o.req.nIter && o.req.deltaMin < o.int.dev
                    o.int.preresult = o.int.result;
                    o.int.keepIdx = ~o.int.removeZE;
                    % perform tikhonov step
                    o.int.result(o.int.keepIdx) = o.int.Const1(o.int.keepIdx) + ...
                        o.int.Const2(o.int.keepIdx, o.int.keepIdx)*o.int.result(o.int.keepIdx);

                    % nonnegativity step
                    if o.req.nonNeg %&& o.int.i < (o.req.nIter-1)
                        o.int.result(o.int.result<0) = 0;
                    end
                    
                    % perform zero entry removal counting
                    if o.int.ZER
                        % get all zero entries
                        o.int.zero(~o.int.removeZE) = o.int.result(~o.int.removeZE) == 0;
                        % reset all nonzero entries in the counter
                        o.int.ZERcounter(~o.int.zero) = 0;
                        % increment all zero entries in the counter
                        o.int.ZERcounter(~o.int.removeZE) = o.int.ZERcounter(~o.int.removeZE) + o.int.zero(~o.int.removeZE);
                        % mark zero entries if threshold is hit
                        o.int.removeZE(~o.int.removeZE) = o.int.removeZE(~o.int.removeZE) | ...
                            o.int.ZERcounter(~o.int.removeZE) == o.aux.neglectZeroEntryLimit;
                    end

                    % calculate deviation to previous step
                    o.int.dev = norm(o.int.result-o.int.preresult)/norm(o.int.preresult)*100;
                    if isnan(o.int.dev)
                        o.int.dev = inf;
                    end

                    % increase iteration count
                    o.int.i = o.int.i + 1;

                    % print algorithm progress
                    o.printProgress();
                end
            end
            if o.aux.verbose
                fprintf('Done!\n');
            end
        end
        
        function o = printProgress(o)
            if o.aux.verbose
                if o.int.i >= o.int.vi || o.int.i == o.req.nIter || o.int.dev < o.req.deltaMin
                    % update runtime
                    o.int.toc = toc(o.int.tic);
                    o.info.runtime = o.info.runtime + o.int.toc;
                    o.int.tic = tic;

                    % print informational stats
                    if isempty(o.aux.xOrig) || all(o.aux.xOrig(:)==0)
                        fprintf([o.info.name '\tIteration %.0f/%.0f\tDeviationPrev: %.2e/%.2e\tRuntime: %.2f\n'],...
                            o.int.i, o.req.nIter, o.int.dev, o.req.deltaMin, o.info.runtime);
                    else
                        currcorr = corr(o.aux.xOrig(:), o.int.result(:));
                        if o.int.i == 1
                            o.int.fillIdx = 1;
                        else
                            o.int.fillIdx = o.int.fillIdx + 1;
                        end
                        o.int.corrstorage(2,o.int.fillIdx) = currcorr;
                        o.int.corrstorage(1,o.int.fillIdx) = o.int.i;
                        if isnan(currcorr)
                            o.int.i = o.req.nIter;
                        end
                        fprintf([o.info.name '\tIteration %.0f/%.0f\tDeviationPrev: %.2e/%.2e\tRuntime: %.2f\tCorrelation: %.5f\n'],...
                            o.int.i, o.req.nIter, o.int.dev, o.req.deltaMin, o.info.runtime, currcorr);
                    end

                    % increase verbose counter
                    o.int.vi = o.int.vi + o.aux.OutputInterval;
                end
            end
        end
        
    end
end

