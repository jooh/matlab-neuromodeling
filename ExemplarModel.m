% basic exemplar model based on isotropic gaussian tunings with some
% tuningwidth (expressed in FWHM).
%
% instancing:
% m = ExemplarModel(exemplars,tuningwidth)
classdef ExemplarModel < CompModel
    properties
        tuningwidth
    end

    methods
        function m = ExemplarModel(exemplars,tuningwidth,averaging)
            if ~any(nargin)
                % sub-class and Saveable support
                exemplars = [];
                tuningwidth = [];
                averaging = [];
            end
            m = m@CompModel(exemplars,averaging);
            % store tuning in sigma rather than FWHM
            m.tuningwidth=tuningwidth./(2*sqrt(2*log(2)));
        end

        function resp = getrawresponse(self,x,y)
            % you can just get the product of the response in each feature
            % dimension to recover the response.
            % So it goes like this - responses for each dimension
            diffs = exp(-bsxfun(@rdivide,bsxfun(@minus,x,y).^2,...
                2*self.tuningwidth.^2));
            % and prodded for the total multivariate response
            resp = prod(diffs,2);
        end
    end
end
