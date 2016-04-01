function res = isoctave()
	res = exist('OCTAVE_VERSION', 'builtin') ~= 0;