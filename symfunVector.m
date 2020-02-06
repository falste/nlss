function res = symfunVector(name, size, symvar)
% res = SYMFUNVECTOR(name, size, symvar)
%	Initializes a vector with symfun elements 'namei(symvar)'

	res = sym(name, [size, 1]);
	for i = 1:size
		
		res(i) = str2sym([name num2str(i) '(' convertStringsToChars(string(symvar)) ')']);
	end
end