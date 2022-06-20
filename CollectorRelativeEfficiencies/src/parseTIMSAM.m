function method = parseTIMSAM(filename)
% PARSEXML Convert XML file to a MATLAB structure.
% Code from Mathworks website, xmlread documentation
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

% take care of header
header = theStruct(1);
method.header(1).Name = 'Filename';
method.header(2).Name = 'DateModified';
method.header(3).Name = 'DateCreated';
method.header(4).Name = 'CreatedBy';
method.header(5).Name = 'ModifiedBy';

for iField = 1:5
    method.header(iField).Value = header.Children(iField).Children.Data;
end % for iField of header

% take care of settings
setting = theStruct(2);
for iField = 1:20
    method.settings(iField).Name = setting.Children(iField).Name;
    
    if ~isempty(setting.Children(iField).Children)
        method.settings(iField).Value = setting.Children(iField).Children.Data;
    end % if ~isempty

end % for iField for settings

BLcount = 0; OPcount = 0;
% baselines and on-peaks
for iField = 2 : size(theStruct,2)-1 % last is 'equilibration'
    
nodeType = string(theStruct(iField).Name);
thisNode = theStruct(iField).Children;

switch nodeType

    case "BASELINE"
    BLcount = BLcount + 1;
    BLname = string(thisNode(3).Children.Data);
    method.baselines(BLcount).Name = BLname;
    BL = struct('Name', [], 'Value', []);
    for iProperty = 1:12
        BL(iProperty).Name = string(thisNode(iProperty).Name);
        if ~isempty(thisNode(iProperty).Children)
            BL(iProperty).Value = thisNode(iProperty).Children.Data;
        end % if
    end % for iProperty
    method.baselines(BLcount).Info = BL;

    case "ONPEAK"
    OPcount = OPcount + 1;
    OPname = "S" + string(thisNode(1).Children.Data);
    method.onpeaks(OPcount).Name = OPname;
    OP = struct('Name', [], 'Value', []);
    for iProperty = 1:17
        OP(iProperty).Name = string(thisNode(iProperty).Name);
        if ~isempty(thisNode(iProperty).Children)
            OP(iProperty).Value = thisNode(iProperty).Children.Data;
        end % if
    end % for iProperty
    method.onpeaks(OPcount).Info = OP;

end % switch nodeType 


end % for iField

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

