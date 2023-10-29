function svcm_fromfile(infile, outfile)
load(infile);
[beta, pval, omega, tscore] = svcm(imgData, xMatrix);
save(outfile, 'beta', 'pval', 'omega', 'tscore');
