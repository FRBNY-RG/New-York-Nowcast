function spec = load_spec(specfile, disp_flag)
% load_spec Load model specification for a dynamic factor model (DFM)

% Description:
%   Loads model specification 'spec' from excel workbook file given by
%   'specfile'.
% 
% Input Arguments:
%   specfile - 
%   disp_flag - boolean to display summary of model specification
%
% Output Arguments:
%   spec - 1 x 1 structure with the following fields:
%     . series_id
%     . name
%     . frequency
%     . units
%     . transformation
%     . category
%     . blocks
%     . BlockNames


% Preallocate output
spec = struct();

% Read data from file
raw    = readtable(specfile);
header = raw.Properties.VariableNames;
raw    = table2cell(raw);

% Parse fields given by column names in Excel worksheet
field_name = {'SeriesID', ...
              'SeriesName', ...
              'Frequency', ...
              'Units', ...
              'Transformation', ...
              'Category'
              };
   
for i_field = 1:numel(field_name)
    fld   = field_name{i_field};
    i_col = find(strcmpi(fld, header), 1);
    if isempty(i_col)
        error([fld ' column missing from model specification.']);
    else
        spec.(fld) = raw(:, i_col);
    end
end

% Include text for transformed units
spec.UnitsTransformed = spec.Transformation;
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'lin', 'Levels (No Transformation)');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'chg', 'Change (Difference)');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'ch1', 'Year over Year Change (Difference)');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'pch', 'Percent Change');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'pc1', 'Year over Year Percent Change');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'pca', 'Percent Change (Annual Rate)');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'cch', 'Continuously Compounded Rate of Change');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'cca', 'Continuously Compounded Annual Rate of Change');
spec.UnitsTransformed = strrep(spec.UnitsTransformed, 'log', 'Natural Log');

% Parse trend
i_col_trend = strcmpi('Trend', header);
Trend       = cell2mat(raw(:, i_col_trend));
Trend(strcmp(spec.Transformation, 'pch')) = Trend(strcmp(spec.Transformation, 'pch'))/12; % Scale down trend for monthly series
spec.Trend  = Trend;

% Parse blocks
i_col_block = strncmpi('Block', header, length('Block'));
Blocks      = cell2mat(raw(:, i_col_block));
Blocks(Blocks == 1) = NaN; % Set unrestricted loadings to NaN
Blocks(Blocks > 1)  = 1;   % Set normalizing loadings to 1
spec.Blocks = Blocks;

% Parse prior
i_col_prior = strcmpi('Prior', header);
Prior       = cell2mat(raw(:, i_col_prior));
spec.Prior  = Prior;

% Sort all fields of 'Spec' in order of decreasing frequency
frequency   = {'d', 'w', 'm', 'q', 'sa', 'a'};
permutation = [];
for i_freq = 1:numel(frequency)
    permutation = [permutation; find(strcmp(frequency{i_freq}, spec.Frequency))]; %#ok<AGROW>
end
field_name  = fieldnames(spec);
for i_field = 1:numel(field_name)
    fld        = field_name{i_field};
    spec.(fld) = spec.(fld)(permutation, :);
end

% Define block names
spec.BlockNames = regexprep(header(i_col_block), 'Block\d+_', '');

% Create category names
spec.CategoryNames = unique(spec.Category);

% Summarize model specification
if (nargin < 2), disp_flag = false; end
if disp_flag
    fprintf('Table 1: Model specification \n');
    try
        tabular = table(spec.SeriesID, spec.SeriesName, spec.Units, spec.UnitsTransformed, ...
            'VariableNames', {'SeriesID', 'SeriesName', 'Units', 'Transformation'});
        disp(tabular)
    catch
    end

end

end