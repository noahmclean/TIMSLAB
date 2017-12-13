function reducedChiSq = optimizePbStdDeadtime(dt, stnd, runs, MSmethod)
% OPTIMIZE982DEADTIME evaluates fit to NBS982 data for a proposed deadtime dt
% where dt is scalar, in nanoseconds.  
%
% This function for use with fminsearch() to find optimum deadtime given data.
%
%       nbs982 is a structure containing nbs982.r46, nbs982.r76, nbs.r86 'true' values
%       runs is the structure of run data containing dataDT0, standard, skips
%       MSmethod is structure containing mass spec method info for beam interpolation
% 
%% 1. Apply deadtime dt and perform beam interpolation

reducedChiSq = 0; % initialize sum at zero

nruns = size(runs,1);

dtcorr = struct('data', 0, 'ratiosBI', 0);

for irun = 1:nruns
    
    % dead time correction
    dtcorr(irun).data = runs(irun).dataDT0 ./ (1 - (dt*1e-9)*runs(irun).dataDT0); 
    
    % beam interpolation
    switch MSmethod.BImethod
        case 'Dodson'
            dtcorr(irun).ratiosBI = DodsonBI_v1(dtcorr(irun).data, MSmethod);
        case 'Quadrift'
            dtcorr(irun).ratiosBI = QuadDriftCorr_v1(dtcorr(irun).data, MSmethod);
    end % switch beam interpolation method
    
    % eliminate skipped ratios from calculation
    if ~isempty(runs(irun).skips) % if the user brushed some data
        dtcorr(irun).ratiosBI = dtcorr(irun).ratiosBI(~runs(irun).skips,:);
    end
    
end % for runi


%% 2.  For 982: Calculate fractionation, correct 204/206 and 208/206
%      For 981: Calculate fractionation for 204/206, 207/206, 208/206

for irun = 1:nruns
    
    switch runs(irun).standard
        
        case 'NBS981'
            standard = stnd.nbs981;
        case 'NBS982'
            standard = stnd.nbs982;
    end
    
            % exponential law fractionation correction factor using 208/206
            dtcorr(irun).beta = (log(standard.r86)-log(dtcorr(irun).ratiosBI(:,4))) ./ ...
                                log(stnd.massPb206/stnd.massPb208);

            % log(204/206) and log(207/206) corrected using beta
            dtcorr(irun).lr46 = log(dtcorr(irun).ratiosBI(:,1)) + ...
                                dtcorr(irun).beta * log(stnd.massPb206/stnd.massPb204);
            dtcorr(irun).lr76 = log(dtcorr(irun).ratiosBI(:,3)) + ...
                                dtcorr(irun).beta * log(stnd.massPb206/stnd.massPb207);


            % simple initial distance formulation: could be improved
            
            % moving standard deviation calculation, for data weighting
            dtcorr(irun).std46 = movstd(dtcorr(irun).lr46, MSmethod.cyclesPerBlock - 1);
            dtcorr(irun).std76 = movstd(dtcorr(irun).lr76, MSmethod.cyclesPerBlock - 1);
            
            % reduced chi-square values for each measured, dt & fractionation corr ratio
            dtcorr(irun).x2red_46 = (dtcorr(irun).lr46 - log(stnd.nbs982r46)).^2 ./ ...
                                     dtcorr(irun).std46.^2;
            dtcorr(irun).x2red_76 = (dtcorr(irun).lr76 - log(stnd.nbs982r76)).^2 ./ ...
                                     dtcorr(irun).std76.^2;
                                 
            reducedChiSq = reducedChiSq + sum(dtcorr(irun).x2red_46) ...
                                        + sum(dtcorr(irun).x2red_76);
    
    
end % for irun
    
    
