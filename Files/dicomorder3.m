function [imgstr, infstr, dim_labels] = dicomorder3(im, sliceinfo)
%%%% Description:                     Order images DICOM v4.2 - 7T UMCU
%%% Creator: Ir. Q. van Houtum       Version: 3.0        Date: 2016-10
%%% --------------------------------------------------------------------
%%% Output variables:
%%% [image structure, info structure] = dicomorder(.*);
%%%
%%% Input arguments:
%%%         .* No input     - open DCM via file-select GUI
%%%         .* Filepath     - String to dicom file
%%%         .* Dicom images, Slice info
%%%                         - DCM image array: WxHxN with N is # of images.
%%%                           As loaded from matlabs dicomread().
%%%                         - Slice info: info-struct per slice as 
%%%                           cell(1,N) loaded via dicomread-function or
%%%                           dicominfo.
%%%
%%% Dimension ordered: stack - slice - echo - dyn
%%%
%%% ! Cardiac is unfortunately excluded for now !
%%% Tested for UMC Utrecht 7T Philips MRI (R3/R5-Mtx) eDICOM export.-v4.2
%%%
%%% Update notes: v2.0 -> 3.0 changed underlying index mechanism. Increased
%%% robustness as number of image types is now variable at output. Speed
%%% slightly increased, requires use of extractFields function (will be
%%% added to this script).


%% Process input %%
switch nargin 
    case 0, [im, sliceinfo, ~] = dicomread7T({});                          
%         if ~exist('im', 'var') || ~exist('sliceinfo', 'var')
%             load('dicomorder3c.mat');
%         end
    case 1, [im, sliceinfo, ~] = dicomread7T({dcmimg});
    case 2 % Do nothing
end


%% Image type indexing %%

% Size of array and other index dimensions
imsz = size(im); fprintf(['Number of images|info-structs:' ...
    '  %d|%d\nImage types \n'], imsz(end),size(sliceinfo,2));

% Get image type labels of all slices
type_perslice = ...
    extractField(sliceinfo, 'PrivatePerFrameSq.Item_1.MRImageTypeMR');
% Find unique image type labels
type_unique = unique(type_perslice);

% Count # of slices for each found unique image type.
type_count = NaN(1,size(type_unique,2));
for type_iter = 1:size(type_unique,2)
    type_count(1,type_iter) = ...
        sum(strcmp(type_perslice, type_unique(type_iter)));
end

% Create output/storage variables
% Each unique image type will get its own struct-field
imgstr = struct; infstr = struct;

%%                      Loop each unique image-type                     %%
for type_iter = 1:size(type_unique,2)
    % Get img/info index of unique image-type(iter) 
    ind = find(strcmp(type_perslice, type_unique(type_iter)) == 1);
    
    % Get ALL slice information of these type-unique slices
    slinfo  = sliceinfo(ind); 
    
    % n_stacks = size(fieldnames(slinf.Stack),1);
    
    % Get dimension-values of this slice.
    nr_stack  = str2double(extractField(slinfo,...
        'FrameContentSequence.Item_1.StackID'));
    nr_echo   = cell2mat(extractField(slinfo,...
        'PrivatePerFrameSq.Item_1.EchoNumber'));
    nr_slice  = double(cell2mat(extractField(slinfo,...
        'PrivatePerFrameSq.Item_1.ImagePlaneNumber')));                     
    nr_dyn    = double(cell2mat(extractField(slinfo,...
        'FrameContentSequence.Item_1.TemporalPositionIndex'))); 
    
    % N.B if statement above: Increase dyn_nr plus 1 if one if indexed at
    % zero! (Some B1 maps show this behaviour for the dynamics).
    if sum(nr_dyn == 0) ~= 0, nr_dyn = nr_dyn+1; end 
    
    % Save dimension info for all slices of unique_type(iter).
    dim_list     = cat(2,nr_stack', nr_slice', nr_echo', nr_dyn');
    dim_list_uni = unique(dim_list, 'rows', 'stable');
    
            %%           Catch: Multiple stacks           %% 
    % Indexing is impaired, too many empty array entries are made in both 
    % the image info and image data array.
    % Why:
    % Slice nr is a continious count ignoring stack indexing. Thus though a
    % slice is another stack, the slice nr is +1 from the last slice in the
    % previous stack.
     
    % Backup of dimension e.g. index list.
    % dim_list_bu = dim_list;
    
    % Nr of slices per stack
    slices_per_stack = NaN(1,max(nr_stack));
    % Nr of stacks
    nst = max(unique(dim_list(:,1))); 
    
    % The indexing will generate more array space than slices available.
    % This is due stacks!
    % Change slice nr per stack in dim_list col 1 = stack; col 2 = slice
    if prod(max(dim_list)) > size(slinfo,2)
        for sti = 1:nst % For each stack
            % stack #sti at rows r
            [r,~] = find(dim_list(:,1) == sti); 
            % N slices in stack
            % nsl = size(r,1); 

            % Create seperate array with only indexes for stack sti
            tmp = dim_list(r,:); 
            % Get sorting index if slices are sorted from lowest to highest
            [~, sorted_ind] = sort(dim_list(r,2), 1); 

            % Loop each row index pointing to this stack
            for sli  = 1:size(r,1)
                % Find in the sorted and seperated index array the current 
                % index in the index loop: its found index is the new 
                % slice nr
                [r2,~] = find(tmp(sorted_ind,2) == dim_list(r(sli),2));
                % Replace the general slice nr with the stack slice nr.
                dim_list(r(sli),2) = r2;
            end
        end 
    end


    
    
    %%           Catch: Image_type dimensions not unique           %%
    
    % If dimensions/index list is not unique
    if size(dim_list,1) ~= size(dim_list_uni,1)
        % List has duplicates: This overwrites image data in outp var!! 

        % Add dimension by counting each duplicate in de dim-list.
        dim_list_uni_new = NaN(size(dim_list,1),1);
        for iter = 1:size(dim_list_uni,1)           
           % Find (boolean) unique-dim in full dim list
           eqdim_bln = ismember(dim_list,dim_list_uni(iter,:),'rows');
           % Get index of reoccurent dim-list-rows
           eqdim_ind = find(eqdim_bln == 1);
           for ind_iter = 1:size(eqdim_ind,1)
               dim_list_uni_new(eqdim_ind(ind_iter)) = ind_iter;
           end
        end
        
        % Add the new uniquified index-column to the dim list.
        dim_list = cat(2,dim_list,dim_list_uni_new);
        
        
    else
        % List is unique, add blank 5th dimension
        dim_list = cat(2,dim_list, ones(size(dim_list,1),1));
    end
    
    
    %%                  Continue ordering by dimensions             %%
    % Create storage variable
    imgstr.(type_unique{type_iter}) = NaN([imsz(1:2), max(dim_list,[],1)]);
    infstr.(type_unique{type_iter}) = cell([max(dim_list,[],1)]);
    
    % Loop dimension-index, rows in dim_list 
    % - doesnt Matlab let you assign a WxHxN array to mutliple 
    % multidimensional indexes [1 2 3 1; 1 3 3 2]?    
    for dim_iter = 1:size(dim_list,1)
        % Fill image-structure imgstr with field type_unique(iter) and add
        % the images from input array im with the index found for this
        % iteration specific type_uniqe(iter) at the nr_dim indexes found
        % above, stored in dim_list.
        imgstr.(type_unique{type_iter})...
            (:,:,dim_list(dim_iter,1),dim_list(dim_iter,2), ...
                 dim_list(dim_iter,3),dim_list(dim_iter,4), ...
                 dim_list(dim_iter,5)) = im(:,:,ind(dim_iter));
        
        % Same action for info-structure infstr except no width nor height.
        infstr.(type_unique{type_iter})...
            {dim_list(dim_iter,1),dim_list(dim_iter,2), ... 
             dim_list(dim_iter,3),dim_list(dim_iter,4), ...
                 dim_list(dim_iter,5)} = sliceinfo{ind(dim_iter)};
    end
    
    
    %%                         Display info                         %%
    maxdim = max(dim_list,[],1);
    fprintf(['%3s - %d stack(s), %d slice(s),'...
        ' %d echo(es), %d dynamic(s), %d optional, total: %d.\n'],...
        type_unique{type_iter}, maxdim(1),maxdim(2), maxdim(3),maxdim(4),...
        maxdim(5), prod(maxdim));
    
end
dim_labels = {'row','col','stack','slice','echo','dyn','opt'};

%%                          Output check                                %%

% - Img stucture check
% Check for consistency due algorithm. Count all images in output and each
% array element to confirm outp has all elements of input. Infocheck not
% necessary as indexing is linked to both output variables. If image
% structure passes the test, so will info, and vice versa.
outpfields = fields(imgstr); outp_tot = 0; outp_sl = 0;
for field_iter = 1:size(outpfields,1)
   outp_tot = outp_tot + numel(imgstr.(outpfields{field_iter}));
   outp_sz = size(imgstr.(outpfields{field_iter}));
   outp_sl = outp_sl + prod(outp_sz(3:end));
end
fprintf(['Processed: %d.\n'],outp_sl);
if numel(im) ~= outp_tot || imsz(end) ~= outp_sl
   warning('Number of output arrays not equal to input!');
end




function fieldvals = extractField(istruct, ifield)
%%%% Description:                     Extract multiple struct fields 
%%% Creator: Ir. Q. van Houtum       Version: 2.1          Date: 2017-04
%%% --------------------------------------------------------------------
%%%
%%% If istruct is a cell-struct (e.g. istruct{:,:,:,:}), extractField will
%%% get all values from the ifield from each cell-index.
%%%
%%% See also: isfieldfull();
%%%
%%% Supported for up to 5 index-dimensions.
%%% Supported for graphics object in Matlab 2016a. Uses fieldnames and
%%% isfield combined to figure out existance of subfields in struct-cell.

% Get field-depth field of interest in structure!
% 1. Get all subfields in requested ifield.
dotindex = strfind(ifield, '.'); subfields = cell(1,size(dotindex,2));
for qq = 1:size(dotindex,2)+1
    if     qq == 1
        if ~isempty(dotindex)
            subfields{qq} = ifield(1:dotindex(qq)-1);
        else
            subfields{qq} = ifield;
        end
    elseif qq == size(dotindex,2)+1
        subfields{qq} = ifield(dotindex(end)+1:end);
    else
        subfields{qq} = ifield(dotindex(qq-1)+1:dotindex(qq)-1);  
    end
end

% Check if all fields exist in istruct.
tmpstruct = istruct{1,1,1,1,1};
for qq = 1:size(subfields,2)
    % Get fieldnames if possible
    fnames = fieldnames(tmpstruct); 
    fnindex = find(strcmp(fnames,subfields{qq}),1);
    if isfield(tmpstruct, subfields{qq}) || ~isempty(fnindex) 
        tmpstruct = tmpstruct.(subfields{qq});
    else warning('RealWorldSys:NoExtractField',...
        'Extracted field(s) not in structure - values not extracted'); 
        fieldvals = []; return;
    end
end

% Get size of field-values ---> use tmpstruct
valsz = size(tmpstruct);
% Storage variable
fieldvals = cell(size(istruct)); % 1 x N values - stslecdy

% Loop through expected 4 dimensions
for st = 1:size(istruct,1)
    for sl = 1:size(istruct,2)
        for ech = 1:size(istruct,3)
            for dyn = 1:size(istruct,4)
                for dum = 1:size(istruct,5) % Dummy dim.
                    tstruct = istruct{st,sl,ech,dyn,dum};
                    for qq = 1:size(subfields,2)
                        tstruct = tstruct.(subfields{qq});
                    end
                    fieldvals{st,sl,ech,dyn,dum} = tstruct; 
                end
            end
        end
    end
end



