function method = parseTIMSAM(filename)
% PARSEXML Convert XML file to a MATLAB structure.
% Code from Mathworks website, xmlread documentation
% 
% NOTE 4/4/23 - can replace all with x = readstruct(methodName, 'FileType', 'xml')
% and then change some field names for compatibility
% e.g. use https://www.mathworks.com/matlabcentral/fileexchange/28516-renamefield

try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

try 
    removeIndentNodes(tree.getChildNodes)
catch 
    error('Failed to remove pesky #text nodes from excess whitespace')
end


% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end

% now collapse child nodes
theStruct = theStruct.Children; % discard top level node

% method name
method.methodName = extractBefore(filename, ".");

% parse names/specifiers for each xml block
nodeNames = string({theStruct.Name});

% create header
thisNode = theStruct(nodeNames == "HEADER").Children;
nFields = length(thisNode);
for iField = 1:nFields
    nameOfField = thisNode(iField).Name;
    method.header.(nameOfField) = thisNode(iField).Children.Data;
end % for iField of header

% take care of settings
thisNode = theStruct(nodeNames == "SETTINGS").Children;
nFields = length(thisNode);
for iField = 1:nFields

    nameOfField = thisNode(iField).Name;
    if ~isempty(thisNode(iField).Children) % if data exists
        method.settings.(nameOfField) = thisNode(iField).Children.Data;
    else % no data
        method.settings.(nameOfField) = [];
    end

end % for iField for settings

% baselines
BLindices = find(nodeNames == "BASELINE");
nBLs = length(BLindices);

for iBL = 1:nBLs
    BLindex = BLindices(iBL);

    thisNode = theStruct(BLindex).Children;
    nFields = length(thisNode);
    for iField = 1:nFields

        nameOfField = thisNode(iField).Name;
        if ~isempty(thisNode(iField).Children) % if data exists
            method.baselines(iBL).(nameOfField) = thisNode(iField).Children.Data;
        else % no data
            method.baselines(iBL).(nameOfField) = [];
        end

    end % for iField in iBL

    % use MassID as "Name" field, e.g. "BL1"
    method.baselines(iBL).Name = method.baselines.MassID;

end % for iBL

% on peaks
OPindices = find(nodeNames == "ONPEAK");
nOPs = length(OPindices);

for iOP = 1:nOPs
    OPindex = OPindices(iOP);

    thisNode = theStruct(OPindex).Children;
    nFields = length(thisNode);
    for iField = 1:nFields

        nameOfField = thisNode(iField).Name;
        if ~isempty(thisNode(iField).Children) % if data exists
            method.onpeaks(iOP).(nameOfField) = thisNode(iField).Children.Data;
        else % no data
            method.onpeaks(iOP).(nameOfField) = [];
        end

    end % for iField in iBL

    % use S + 'Sequence' as "Name" field, e.g. "S1"
    method.onpeaks(iOP).Name = "S" + method.onpeaks(iOP).Sequence;

end % for iBL


% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end

% ------ Local function removeIndentNodes 
% from https://stackoverflow.com/a/11552839
% removes (some) whitespace to prevent #text nodes

function removeIndentNodes( childNodes )

numNodes = childNodes.getLength;
remList = [];
for i = numNodes:-1:1
   theChild = childNodes.item(i-1);
   if (theChild.hasChildNodes)
      removeIndentNodes(theChild.getChildNodes);
   else
      if ( theChild.getNodeType == theChild.TEXT_NODE && ...
           ~isempty(char(theChild.getData()))         && ...
           all(isspace(char(theChild.getData()))))
         remList(end+1) = i-1; % java indexing
      end
   end
end
for i = 1:length(remList)
   childNodes.removeChild(childNodes.item(remList(i)));
end

