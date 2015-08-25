classdef RampModel < CompModel
    properties
        saturation
        offset
    end

    methods
        function m = RampModel(exemplars,saturation,offset,averaging)
            if ~any(nargin)
                exemplars = [];
                saturation = [];
                offset = [];
                averaging = [];
            end
            m = m@CompModel(exemplars,averaging);
            m.saturation = saturation;
            m.offset = offset;
        end

        function xy = model2data(self,weights)
            % return coordinates defining the discriminant dimension for
            % each ramp. Each pair of xy coordinates is interspersed with
            % NaN to support a single plot command. The length of the lines
            % is normalised to self.offset
            %
            % xy = model2data(self,weights)
            assert(size(weights,2)<3,['model2data is only supported ' ...
                'for 2D models']);
            % find slope
            sl = weights(:,2)./weights(:,1);
            x = [-1; 1];
            y = x * sl';
            y(3,:) = NaN;
            x = repmat([-1 1 NaN]',[1 size(y,2)]);
            xy = [x(:) y(:)];
            notnanind = ~any(isnan(xy),2);
            % scale length of each line to offset*2
            xy(notnanind,:) = unitlen(xy(notnanind,:)')' * self.offset;
        end

        function rawresp = getrawresponse(self,weights,coords)
            % project coordinates onto weights. So basically convert the ND
            % coordinate to a 1D vector.
            % This is done either for the populationresponse case (many
            % weights, 1 coordinate) or for the unitresponse case (1 weight
            % set, many coordinates).
            discdata = coords * weights';
            % and then we pass all those values through the sigmoid to
            % obtain the response
            rawresp = sigmoidlike(discdata,self.saturation,self.offset);
        end

    end
end
