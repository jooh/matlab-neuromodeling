classdef CompModel < Saveable
    properties
        exemplars
        averaging
        poprespcachecoords = []
        poprespcache = []
    end

    methods 
        
        function m = CompModel(exemplars,averaging)
            if ~any(nargin)
                exemplars = [];
                averaging = [];
            end
            m.exemplars = exemplars;
            m.averaging = averaging;
        end

        function x = data2model(self,x)
        % return the model representation of the input x, which is in data
        % space.
        %
        % x = data2model(self,x)
            x = x;
        end

        function x = model2data(self,x)
        % return some representation of the exemplar x in data space.
        % For an exemplar model this would be the tuning center.
        %
        % x = model2data(self,x)
            x = x;
        end

        function resp = getrawresponse(self,varargin)
            resp = [];
            error('no getrawresponse method in base class')
        end

        function resp = getresponse(self,weights,coords)
            % get the unnormalised response pattern (either across
            % coordinates or across units)
            rawresp = getrawresponse(self,weights,coords);
            % get the mean activation across the population to the current
            % coords
            popresp = getpopulationresponse(self,coords);
            % interpolate between population mean and pattern.
            resp = ascol((rawresp-popresp)*(1-self.averaging) + popresp);
        end

        function popresp = getpopulationresponse(self,coords)
        % get the mean response to each coordinate in coords across the
        % population of units (weights) in self.
        %
        % popresp = getpopulationresponse(self,coords)
            if isequal(coords,self.poprespcachecoords)
                popresp = self.poprespcache;
                return
            end
            % need to iterate because sometimes we have multiple coordinate
            % inputs (and we can't directly map many coordinates to many
            % weights).
            popresp = NaN([size(coords,1) 1]);
            for n = 1:size(coords,1)
                % get the raw response to this coordinate and average it
                % over the population
                popresp(n) = mean(getrawresponse(self,self.exemplars,coords(n,:)));
            end
            self.poprespcache = popresp;
            self.poprespcachecoords = coords;
        end

        function resp = populationresponse(self,testex)
        % population response (rows) to a single test exemplar.
        %
        % resp = populationresponse(self,testex)
            % response in scaled space
            resp = getresponse(self,self.exemplars,...
                data2model(self,testex));
        end

        function resp = unitresponse(self,exind,testex)
        % single unit response to a set of test exemplars (rows).
        %
        % resp = unitresponse(self,testex,exind)
            % scale the test patterns into the model space
            testex = data2model(self,testex);
            % get responses in model space
            resp = getresponse(self,self.exemplars(exind,:),testex);
        end

        function [h,ph,intmap,cmap] = plotunit(self,exind,varargin)
        % plot the response profile of a single unit.
        %
        % [h,ph,intmap,cmap] = plotunit(self,exind,varargin)
            assert(size(self.exemplars,2)<3,'can only plot 2D models');
            getArgs(varargin,{'xvals',[],'yvals',[],'ax',gca,...
                'cmap',gray(1024),'limits','zerobounded','gridlines',[],...
                'gridcolor',[],'outline',[],'outlinecolor',[],...
                'plotcolor',[1 1 1]});
            if isempty(xvals) && isempty(yvals)
                % auto-lim centered on 0
                l = max(abs(ascol(self.exemplars)));
                xvals = linspace(-l,l,100);
            end
            if isempty(xvals) && ~isempty(yvals)
                xvals = yvals;
            end
            if isempty(yvals) && ~isempty(xvals)
                yvals = xvals;
            end
            [x,y] = meshgrid(xvals,yvals);
            resp = reshape(unitresponse(self,exind,[x(:) y(:)]),size(x));
            [im,intmap,cmap] = intensity2rgb(resp,cmap,limits);
            h = imageplot(ax,im,'gridlines',gridlines,...
                'gridcolor',gridcolor,'outline',outline,...
                'outlinecolor',outlinecolor,'upscale=1',...
                'ylims',[yvals(1) yvals(end)],...
                'xlims',[xvals(1) xvals(end)]);
            % overlay the unit center
            xy = model2data(self,self.exemplars(exind,:));
            hold(ax,'on');
            ph = plot(xy(:,1),xy(:,2),'o-','markeredgecolor',plotcolor,...
                'color',plotcolor);
        end

        function [h,ph,intmap,cmap] = plotdissimilarity(self,xy,varargin)
        % plot the dissimilarity profile relative to a reference coordinate
        % xy.
        %
        % [h,ph,intmap,cmap] = plotdissimilarity(self,xy,varargin)
            assert(size(self.exemplars,2)<3,'can only plot 2D models');
            getArgs(varargin,{'xvals',[],'yvals',[],'ax',gca,...
                'cmap',cmap_wr,'limits','zerobounded','gridlines',[],...
                'gridcolor',[],'outline',[],'outlinecolor',[],...
                'plotcolor',[0 0 0]});
            if isempty(xvals) && isempty(yvals)
                % auto-lim centered on 0
                l = max(abs(ascol(self.exemplars)));
                xvals = linspace(-l,l,100);
            end
            if isempty(xvals) && ~isempty(yvals)
                xvals = yvals;
            end
            if isempty(yvals) && ~isempty(xvals)
                yvals = xvals;
            end
            [x,y] = meshgrid(xvals,yvals);
            % get the responses to all those positions
            responses = arrayfun(@(ind)populationresponse(self,...
                [x(ind) y(ind)]),1:numel(x),'uniformoutput',0);
            refresponse = populationresponse(self,xy);
            % euclidean distances after making the reference xy the norm of
            % the space
            distances = reshape(sqrt(sum(bsxfun(@minus,...
                cat(2,responses{:}),refresponse).^2,1)),size(x));
            [im,intmap,cmap] = intensity2rgb(distances,cmap,limits);
            h = imageplot(ax,im,'gridlines',gridlines,...
                'gridcolor',gridcolor,'outline',outline,...
                'outlinecolor',outlinecolor,'upscale=1',...
                'ylims',[yvals(1) yvals(end)],...
                'xlims',[xvals(1) xvals(end)]);
            % overlay the unit center
            hold(ax,'on');
            ph = plot(xy(1),xy(2),'o','markeredgecolor',plotcolor);
        end
    end
end
