function S = phononShearSuite(pos, bonds, masses, latticeVectors, basisFrac, k, gammaMax, nSteps, alphaMorse)
%PHONONSHEARSUITE Run three shear-sweep models on the lattice exactly as the app built it (built-in or any CIF).
%  S = phononShearSuite(pos, bonds, masses, latticeVectors, basisFrac, k, gammaMax, nSteps, alphaMorse)   (alphaMorse default 4)
%  cases{1} harmonic, affine (= applyShearForce / sliders); cases{2} + pre-stress + relaxation;
%  cases{3} Morse + pre-stress + relaxation + bond rupture (main model, S.main = 3).
%  Shared by the app buttons and runPhononShearAnalysis.m so both give identical results.

    if nargin < 7 || isempty(gammaMax), gammaMax = 0.5; end
    if nargin < 8 || isempty(nSteps),   nSteps = 40;   end
    if nargin < 9 || isempty(alphaMorse), alphaMorse = 4;  end
    base = struct('latticeVectors', latticeVectors, 'basisFrac', basisFrac, ...
        'gammaMax', gammaMax, 'nSteps', nSteps, 'k', k);

    oA = base; oA.model = 'harmonic'; oA.prestress = false; oA.allowBreak = false; oA.relax = false;
    oB = base; oB.model = 'harmonic'; oB.prestress = true;  oB.allowBreak = false; oB.relax = true;
    oC = base; oC.model = 'morse';    oC.prestress = true;  oC.allowBreak = true;  oC.relax = true;  oC.alphaMorse = alphaMorse;

    S.cases  = {phononShearSweep(pos, bonds, masses, oA), ...
                phononShearSweep(pos, bonds, masses, oB), ...
                phononShearSweep(pos, bonds, masses, oC)};
    S.labels = {'harmonic, affine (as in sliders)', 'harmonic + pre-stress + relaxation', 'Morse + pre-stress + relaxation + rupture'};
    S.main = 3;  S.gammaMax = gammaMax;  S.k = k;  S.alphaMorse = alphaMorse;
end
