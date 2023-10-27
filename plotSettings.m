function plotSettings(gcf,varargin)
    % set units
    set(gcf, 'Units', 'centimeters'); 
    
    % set position size mode
    switch nargin
        case 2
            % varargin{1} = afFigurePosition
            set(gcf, 'Position', varargin{1}, 'PaperSize', PaperSize, 'PaperPositionMode', 'auto');
        case 3
            % varargin{1} = afFigurePosition
            % varargin{2} = PaperSize
            set(gcf, 'Position', varargin{1}, 'PaperSize', varargin{2}, 'PaperPositionMode', 'auto');
        otherwise
            set(gcf, 'Position', [2 7 19 12], 'PaperSize', [19 11], 'PaperPositionMode', 'auto');
    end
end

