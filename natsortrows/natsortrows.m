function [B,ndx,dbg] = natsortrows(A,rgx,varargin)
% Natural-order / alphanumeric sort the rows of a cell array, string array, or table.
%
% (c) 2014-2021 Stephen Cobeldick
%
% Sort text by character code and by number value. For a cell/string array:
% - SORTROWS <column> option is supported, selects the columns to sort by.
% - SORTROWS <direction> option is supported, specifies the sort directions.
% For a table/timetable array any categorical, string, or cell-of-character
% variables are sorted alphanumerically, all other types via SORT/SORTROWS.
% For a table/timetable the options above are supported and additionally:
% - SORTROWS 'RowNames' and <rowDimName> options are supported.
% - SORTROWS <vars> option is supported, selects the variables to sort by.
%
%%% Example:
% >> X = {'A2','X';'A10','Y';'A10','X';'A1','X'};
% >> sortrows(X) % SORTROWS for comparison.
% ans =
%     'A1'     'X'
%     'A10'    'X'
%     'A10'    'Y'
%     'A2'     'X'
% >> natsortrows(X)
% ans =
%     'A1'     'X'
%     'A2'     'X'
%     'A10'    'X'
%     'A10'    'Y'
%
%%% Syntax:
%  B = natsortrows(A)
%  B = natsortrows(A,rgx)
%  B = natsortrows(A,rgx,<options>)
% [B,ndx,dbg] = natsortrows(A,...)
%
% To sort the elements of a string/cell array use NATSORT (File Exchange 34464)
% To sort any file-names or folder-names use NATSORTFILES (File Exchange 47434)
%
%% File Dependency %%
%
% NATSORTROWS requires the function NATSORT (File Exchange 34464). Extra
% optional arguments are passed directly to NATSORT. See NATSORT for case-
% sensitivity, sort direction, number format matching, and other options.
%
%% Examples %%
%
% >> A = {'B','2','X';'A','100','X';'B','10','X';'A','2','Y';'A','20','X'};
% >> sortrows(A) % SORTROWS for comparison.
% ans =
%    'A'  '100'  'X'
%    'A'    '2'  'Y'
%    'A'   '20'  'X'
%    'B'   '10'  'X'
%    'B'    '2'  'X'
% >> natsortrows(A)
% ans =
%    'A'    '2'  'Y'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
% >> natsortrows(A,[],'descend')
% ans =
%     'B'  '10'  'X'
%     'B'   '2'  'X'
%     'A' '100'  'X'
%     'A'  '20'  'X'
%     'A'   '2'  'Y'
%
% >> sortrows(A,[2,-3]) % SORTROWS for comparison.
% ans =
%    'B'   '10'  'X'
%    'A'  '100'  'X'
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'A'   '20'  'X'
% >> natsortrows(A,[],[2,-3])
% ans =
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
% >> natsortrows(A,[],[false,true,true],{'ascend','descend'})
% ans =
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
% >> natsortrows(A,[],{'ignore','ascend','descend'})
% ans =
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
%
% T = cell2table(A);
% natsortrows(T,[], [2,-3]) % TABLE
% ans = 
%     A1      A2      A3 
%     ___    _____    ___
%     'A'    '2'      'Y'
%     'B'    '2'      'X'
%     'B'    '10'     'X'
%     'A'    '20'     'X'
%     'A'    '100'    'X'
% natsortrows(T,[], {'A2','A3'},{'ascend','descend'}) % TABLE
% ans = 
%     A1      A2      A3 
%     ___    _____    ___
%     'A'    '2'      'Y'
%     'B'    '2'      'X'
%     'B'    '10'     'X'
%     'A'    '20'     'X'
%     'A'    '100'    'X'
%
% >> B = {'ABCD';'3e45';'67.8';'+Inf';'-12';'+9';'NaN'};
% >> sortrows(B) % SORTROWS for comparison.
% ans =
%    '+9'
%    '+Inf'
%    '-12'
%    '3e45'
%    '67.8'
%    'ABCD'
%    'NaN'
% >> natsortrows(B,'[-+]?(NaN|Inf|\d+\.?\d*(E[-+]?\d+)?)')
% ans =
%    '-12'
%    '+9'
%    '67.8'
%    '3e45'
%    '+Inf'
%    'NaN'
%    'ABCD'
%
% >> C = {'A2', 2; 'A10', 1; 'A2', 1}; % Scalar numerics in a column:
% >> natsortrows(C,[],'sortnum')
% ans =
%     'A2'     [1]
%     'A2'     [2]
%     'A10'    [1]
% >> natsortrows(C,[],'sortnum',{'ascend','descend'})
% ans =
%     'A2'     [2]
%     'A2'     [1]
%     'A10'    [1]
%
%% Input and Output Arguments %%
%
%%% Inputs (**=default):
% A   = Array of size MxN, with atomic rows to be sorted. Can be a string
%       array, or a cell array of character row vectors, or a categorical
%       array, or any other array type supported by NATSORT.
%     = Table of size MxN, with atomic rows to be sorted. Columns/variables
%       that are string, categorical, or cell of character row vectors are
%       sorted using NATSORT, all other column types are sorted using SORT.
% rgx = Optional regular expression to match number substrings.
%     = [] uses the default regular expression '\d+'** to match integers.
% <options> can be supplied in any order:
%     = 'SortNum': If sorting a cell array then it may also contain columns
%       of numeric/logical scalars: each such column is concatenated into
%       a vector (type coersion may occur) which is sorted using SORT.
%     = Logical vector indicating which columns of <A> to sort by.
%     = <column>: numeric vector where each integer specifies which columns
%       of <A> to sort by. Negative integers indicate a descending sort.
%     = <direction>: a cell array containing only the character vectors
%       'ascend', 'descend', and/or 'ignore'. The number of cells must match
%       the number of columns being sorted. The sign of <column> is ignored.
% <options> additionally supported for tables/timetables:
%     = 'RowNames' or <rowDimName> (the name of first dimension of
%       table <A>): sorts table <A> based on its row names.
%     = <vars>: a cell array containing the names (character row vectors)
%       of the timetable/table <A> variables to sort by.
% Any remaining <options> are passed directly to NATSORT.
%
%%% Outputs:
% B   = Array <A> with rows sorted into alphanumeric order.
% ndx = NumericVector Mx1 of row indices such that B = A(ndx,:).
% dbg = CellVectorOfCellArrays, size 1xN. Each cell contains the debug cell array
%       for one column of <A>. Helps debug the regular expression (see NATSORT).
%
% See also SORT SORTROWS NATSORT NATSORTFILES TABLE TIMETABLE CELLSTR REGEXP IREGEXP SSCANF

%% Input Wrangling %%
%
fnh = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
%
varargin = cellfun(@nsr1s2c, varargin, 'UniformOutput',false);
%
assert(ndims(A)<3,...
	'SC:natsortrows:A:NotMatrix',...
	'First input <A> must be a matrix (2D).') %#ok<ISMAT>
%
[nmr,nmc] = size(A);
ndx = 1:nmr;
dbg = cell(1,nmc);
%
cda = {'descend','ascend','ignore'};
iso = false;
ist = isa(A,'table')||isa(A,'timetable');
chk = 'RowNames|SortNum';
%
%% Columns to Sort %%
%
bvc = cellfun(@iscell,varargin); % direction
bvl = cellfun(@islogical,varargin); % column
bvn = cellfun(@isnumeric,varargin); % column
bvv = bvl|bvn; % column
bvt = ~bvc&~bvv; % text
%
txt = varargin(bvt);
assert(all(fnh(txt)),...
	'SC:natsortrows:option:InvalidType',sprintf(...
	'Optional inputs must be character row vectors, or string scalars,\n%s.',...
	'or logical/numeric column indices, or cell arrays of directions/variable names.'))
%
ito = any(strcmpi('descend',txt)|strcmpi('ascend',txt)|strcmpi('ignore',txt));
bt1 = strcmpi('SortNum',txt);
bvt(bvt) = ~bt1;
%
assert(nnz(bvv)<2,...
	'SC:natsortrows:column:Overspecified',...
	'The <column> option is over-specified: one logical or numeric vector is allowed.')
%
if any(bvl)
	axc = find(varargin{bvl});
	assert(max(axc)<=nmc&&isvector(varargin{bvl}),...
		'SC:natsortrows:column:IndexMismatchLogical',...
		'The <column> option must be a vector of logical indices into columns of <A>.')
elseif any(bvn)
	col = varargin{bvn};
	assert(isvector(col)&&isreal(col)&&all(~mod(col,1))&&all(col)&&all(abs(col)<=nmc),...
		'SC:natsortrows:column:IndexMismatchNumeric',...
		'The <column> option must be a vector of subscript indices into columns of <A>.')
	aso = cda((3+sign(col))/2);
	iso = ~ito;
	axc = abs(col);
else
	axc = 1:nmc;
end
%
if ist % table
	prn = A.Properties.RowNames;
	pvn = A.Properties.VariableNames;
	pdn = A.Properties.DimensionNames(1);
	%
	chk = sprintf('|%s','SortNum',prn{:},pvn{:},pdn{:});
	chk = sprintf('RowNames%s',chk);
	%
	btn = ismember(txt,pvn);
	btr = strcmpi(txt,'RowNames')|strcmpi(txt,pdn);
	%
	if any(btr) % sort by table row names
		assert(nnz(btr)<2,...
			'SC:natsortrows:RowNames:Overspecified',...
			'The "RowNames" or <rowDimName> option is over-specified, may be used once.')
		assert(all(bvt) && ~any(btn),...
			'SC:natsortrows:RowNames:NotExclusive',...
			'The "RowNames" or <rowDimName> option cannot be combined with <column> or <var> options.')
		bvt(bvt) = btr;
		varargin(bvt) = [];
		if nargin>1
			nsrChkRgx(rgx,chk)
			varargin = [{rgx},varargin];
		end
		dbg = {[]};
		if numel(prn)
			if nargout<3 % faster:
				[~,ndx] = natsort(prn,varargin{:});
			else % for debugging:
				[~,ndx,dbg] = natsort(prn,varargin{:});
			end
		end
		ndx = ndx(:);
		B = A(ndx,:);
		return
	elseif any(btn) % sort by one variable name
		assert(nnz(btn)<2,...
			'SC:natsortrows:VariableName:Overspecified',...
			'To specify multiple variable names use a cell array of character vectors.')
		assert(all(bvt),...
			'SC:natsortrows:VariableName:NotExclusive',...
			'A variable name cannot be combined with <column> or <var> cell array options.')
		axc = find(strcmp(varargin{btn},pvn));
		bvt(bvt) = ~btn;
	elseif any(bvc) % sort by <direction> and/or <vars>
		prd = false;
		prv = false;
		vnc = find(bvc);
		for k = 1:numel(vnc)
			sbc = varargin{vnc(k)};
			assert(isvector(sbc),...
				'SC:natsortrows:option:CellNotVector',...
				'Optional cell arrays must be vectors.')
			assert(all(fnh(sbc)),...
				'SC:natsortrows:option:CellContentNotChar',...
				'Optional cell arrays must contain only character row vectors.')
			tmp = lower(sbc);
			if all(ismember(tmp,cda)) % <direction>
				assert(~ito && ~prd,...
					'SC:natsortrows:direction:Overspecified',...
					'The <direction> option is over-specified, may be used only once.')
				aso = tmp;
				iso = true;
				prd = true;
			else % <vars>
				assert(~prv,...
					'SC:natsortrows:vars:Overspecified',...
					'The <vars> option is over-specified, may be used only once.')
				assert(~any(bvv),...
					'SC:natsortrows:vars:NotExclusive',...
					'The <vars> option cannot be combined with <column> or RowName options.')
				axc = nan(size(sbc));
				for ii = 1:numel(sbc)
					tmp = strcmp(sbc{ii},pvn);
					assert(nnz(tmp)==1,...
						'SC:natsortrows:vars:UnrecognizedName',...
						'The <vars> option has an unrecognised variable name: "%s".',sbc{ii})
					axc(ii) = find(tmp);
				end
				prv = true;
			end
		end
		if prd
			assert(numel(aso)==numel(axc),...
				'SC:natsortrows:direction:NumberDirections',...
				'The <direction> cell array must specify one direction for each column to be sorted.')
		end
	end
elseif any(bvc) % array
	assert(~ito && nnz(bvc)<2,...
		'SC:natsortrows:direction:Overspecified',...
		'The <direction> option is over-specified, may be used only once.')
	aso = varargin{bvc};
	iso = true;
	assert(isvector(aso),...
		'SC:natsortrows:direction:CellNotVector',...
		'The <direction> option cell array must be a vector.')
	assert(all(fnh(aso)),...
		'SC:natsortrows:direction:CellContentNotChar',...
		'The <direction> cell array must contain only character row vectors.')
	aso = lower(aso);
	assert(all(ismember(aso,cda)),...
		'SC:natsortrows:direction:NotAscendDescend',...
		'The <direction> cell array must contain only ''ascend'' and ''descend''.')
	assert(numel(aso)==numel(axc),...
		'SC:natsortrows:direction:NumberDirections',...
		'The <direction> cell array must specify one direction for each column to be sorted.')
end
%
if any(bt1)
	assert(iscell(A),...
		'SC:natsortrows:SortNum:NotCell',...
		'The "SortNum" option is only supported when <A> is a cell array.')
	ax1 = all(cellfun('prodofsize',A)==1,1) & ... logical/numeric scalars
		(all(cellfun('islogical',A),1) | all(cellfun(@isnumeric,A),1));
else
	ax1 = false(1,nmc);
end
%
if nargin>1
	nsrChkRgx(rgx,chk)
	varargin = [{rgx},cell(1,+iso),varargin(bvt)];
end
%
%% Sort by Column %%
%
for k = numel(axc):-1:1
	axk = axc(k);
	if iso % sort order:
		varargin{2} = aso{k};
	end
	if any(strcmpi(varargin,'ignore'))
		continue
	end
	if ax1(axk) % logical/numeric scalars
		isd = any(strcmpi(varargin,'descend'));
		tmp = vertcat(A{ndx,axk}); % type coersion may occur...
		[~,idx] = sort(tmp,cda{2-isd});
		[~,idb] = sort(ndx);
		dbg{axk} = tmp(idb); %... so we return the numeric vector.
	elseif ist % table
		tmp = A{ndx,axk};
		if isa(tmp,'string')||iscategorical(tmp)||iscell(tmp)&&all(fnh(tmp(:)))
			if size(tmp,2)==1
				if nargout<3 % faster:
					[~,idx] = natsort(tmp,varargin{:});
				else % for debugging:
					[~,idx,trd] = natsort(tmp,varargin{:});
					[~,idb] = sort(ndx);
					dbg{axk} = trd(idb,:);
				end
			else
				[~,idx] = natsortrows(tmp,varargin{:});
			end
		else % numeric, logical, datetime, etc.
			isd = any(strcmpi(varargin,'descend'));
			col = (1-2*isd).*(1:size(tmp,2));
			[~,idx] = sortrows(tmp,col);
		end
	else % char, string, categorical, cell of char vectors, etc.
		if nargout<3 % faster:
			[~,idx] = natsort(A(ndx,axk),varargin{:});
		else % for debugging:
			[~,idx,trd] = natsort(A(ndx,axk),varargin{:});
			[~,idb] = sort(ndx);
			dbg{axk} = trd(idb,:);
		end
	end
	ndx = ndx(idx);
end
%
ndx = ndx(:);
B = A(ndx,:);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortrows
function nsrChkRgx(rgx,chk)
chk = sprintf('^(%s)$',chk);
assert(~ischar(rgx)||isempty(regexpi(rgx,chk,'once')),...
	'SC:natsortrows:rgx:OptionMixUp',...
	['Second input <rgx> must be a regular expression that matches numbers.',...
	'\nThe provided expression "%s" looks like an optional argument (inputs 3+).'],rgx)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsrChkRgx
function arr = nsr1s2c(arr)
% If scalar string then extract the character vector, otherwise data is unchanged.
if isa(arr,'string') && isscalar(arr)
	arr = arr{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsr1s2c