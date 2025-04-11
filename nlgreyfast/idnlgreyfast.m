classdef idnlgreyfast < matlab.mixin.CustomDisplay
% class for nlgreyfast identification results
   properties
      param_est_v %the estimated model parameters
      x0_est_v %the estimated initial states 
      diag_data %contains debugging information
   end
end
