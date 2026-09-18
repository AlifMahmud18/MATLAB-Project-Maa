function [colors, elementNames] = getElementDetails(masses)
%GETELEMENTDETAILS Map atomic masses to element names and CPK colours (all 118 elements).
%  [colors, elementNames] = getElementDetails(masses); names look like 'Calcium (Ca)'.
%  The nearest table mass is accepted within max(0.3, 1% of its mass);
%  otherwise the atom is 'Generic Atom (m = ...)' with a palette colour.
%  Used for atom colours in the views and for the hover label.

    Natoms = length(masses);
    colors = zeros(Natoms, 3);
    elementNames = cell(Natoms, 1);

    elements = {
        1.008,   [0.90 0.90 0.90], 'H',  'Hydrogen';
        4.003,   [0.85 1.00 1.00], 'He', 'Helium';
        6.94,    [0.80 0.50 1.00], 'Li', 'Lithium';
        9.0122,  [0.76 1.00 0.00], 'Be', 'Beryllium';
        10.81,   [1.00 0.70 0.70], 'B',  'Boron';
        12.011,  [0.35 0.35 0.35], 'C',  'Carbon';
        14.007,  [0.10 0.30 0.90], 'N',  'Nitrogen';
        15.999,  [0.90 0.10 0.10], 'O',  'Oxygen';
        18.998,  [0.60 0.90 0.60], 'F',  'Fluorine';
        20.180,  [0.70 1.00 1.00], 'Ne', 'Neon';
        22.990,  [0.60 0.10 0.90], 'Na', 'Sodium';
        24.305,  [0.50 1.00 0.00], 'Mg', 'Magnesium';
        26.982,  [0.80 0.80 0.80], 'Al', 'Aluminum';
        28.085,  [0.90 0.75 0.20], 'Si', 'Silicon';
        30.974,  [1.00 0.50 0.00], 'P',  'Phosphorus';
        32.06,   [1.00 1.00 0.20], 'S',  'Sulfur';
        35.45,   [0.10 0.90 0.10], 'Cl', 'Chlorine';
        39.948,  [0.50 0.82 0.89], 'Ar', 'Argon';
        39.098,  [0.50 0.20 0.80], 'K',  'Potassium';
        40.078,  [0.20 0.80 0.20], 'Ca', 'Calcium';
        44.956,  [0.90 0.90 0.90], 'Sc', 'Scandium';
        47.867,  [0.75 0.76 0.78], 'Ti', 'Titanium';
        50.942,  [0.65 0.65 0.67], 'V',  'Vanadium';
        51.996,  [0.54 0.60 0.78], 'Cr', 'Chromium';
        54.938,  [0.61 0.48 0.78], 'Mn', 'Manganese';
        55.845,  [0.88 0.40 0.20], 'Fe', 'Iron';
        58.933,  [0.94 0.56 0.63], 'Co', 'Cobalt';
        58.693,  [0.30 0.80 0.30], 'Ni', 'Nickel';
        63.546,  [0.78 0.50 0.20], 'Cu', 'Copper';
        65.38,   [0.50 0.50 0.65], 'Zn', 'Zinc';
        69.723,  [0.76 0.56 0.56], 'Ga', 'Gallium';
        72.630,  [0.40 0.56 0.56], 'Ge', 'Germanium';
        74.922,  [0.74 0.50 0.89], 'As', 'Arsenic';
        78.971,  [1.00 0.63 0.00], 'Se', 'Selenium';
        79.904,  [0.65 0.16 0.16], 'Br', 'Bromine';
        83.798,  [0.36 0.72 0.78], 'Kr', 'Krypton';
        85.468,  [0.44 0.18 0.69], 'Rb', 'Rubidium';
        87.62,   [0.00 1.00 0.00], 'Sr', 'Strontium';
        88.906,  [0.58 1.00 1.00], 'Y',  'Yttrium';
        91.224,  [0.58 0.88 0.88], 'Zr', 'Zirconium';
        92.906,  [0.45 0.76 0.79], 'Nb', 'Niobium';
        95.95,   [0.33 0.71 0.71], 'Mo', 'Molybdenum';
        98,      [0.23 0.62 0.62], 'Tc', 'Technetium';
        101.07,  [0.14 0.56 0.56], 'Ru', 'Ruthenium';
        102.91,  [0.04 0.49 0.55], 'Rh', 'Rhodium';
        106.42,  [0.00 0.41 0.52], 'Pd', 'Palladium';
        107.87,  [0.75 0.75 0.75], 'Ag', 'Silver';
        112.41,  [1.00 0.85 0.56], 'Cd', 'Cadmium';
        114.82,  [0.65 0.46 0.45], 'In', 'Indium';
        118.71,  [0.40 0.50 0.50], 'Sn', 'Tin';
        121.76,  [0.62 0.39 0.71], 'Sb', 'Antimony';
        127.60,  [0.83 0.48 0.00], 'Te', 'Tellurium';
        126.90,  [0.58 0.00 0.58], 'I',  'Iodine';
        131.29,  [0.26 0.62 0.69], 'Xe', 'Xenon';
        132.91,  [0.34 0.09 0.56], 'Cs', 'Cesium';
        137.33,  [0.00 0.79 0.00], 'Ba', 'Barium';
        138.91,  [0.44 0.83 1.00], 'La', 'Lanthanum';
        140.12,  [1.00 1.00 0.78], 'Ce', 'Cerium';
        140.91,  [0.85 1.00 0.78], 'Pr', 'Praseodymium';
        144.24,  [0.78 1.00 0.78], 'Nd', 'Neodymium';
        145,     [0.64 1.00 0.78], 'Pm', 'Promethium';
        150.36,  [0.56 1.00 0.78], 'Sm', 'Samarium';
        151.96,  [0.38 1.00 0.78], 'Eu', 'Europium';
        157.25,  [0.27 1.00 0.78], 'Gd', 'Gadolinium';
        158.93,  [0.19 1.00 0.71], 'Tb', 'Terbium';
        162.50,  [0.12 1.00 0.65], 'Dy', 'Dysprosium';
        164.93,  [0.00 1.00 0.61], 'Ho', 'Holmium';
        167.26,  [0.00 0.90 0.46], 'Er', 'Erbium';
        168.93,  [0.00 0.83 0.32], 'Tm', 'Thulium';
        173.05,  [0.00 0.75 0.22], 'Yb', 'Ytterbium';
        174.97,  [0.00 0.67 0.14], 'Lu', 'Lutetium';
        178.49,  [0.30 0.76 1.00], 'Hf', 'Hafnium';
        180.95,  [0.30 0.65 1.00], 'Ta', 'Tantalum';
        183.84,  [0.13 0.58 0.84], 'W',  'Tungsten';
        186.21,  [0.15 0.49 0.67], 'Re', 'Rhenium';
        190.23,  [0.15 0.40 0.59], 'Os', 'Osmium';
        192.22,  [0.09 0.33 0.53], 'Ir', 'Iridium';
        195.08,  [0.82 0.82 0.88], 'Pt', 'Platinum';
        196.97,  [1.00 0.82 0.14], 'Au', 'Gold';
        200.59,  [0.72 0.72 0.82], 'Hg', 'Mercury';
        204.38,  [0.65 0.33 0.30], 'Tl', 'Thallium';
        207.2,   [0.34 0.35 0.38], 'Pb', 'Lead';
        208.98,  [0.62 0.31 0.71], 'Bi', 'Bismuth';
        209,     [0.67 0.36 0.00], 'Po', 'Polonium';
        210,     [0.46 0.31 0.27], 'At', 'Astatine';
        222,     [0.26 0.51 0.59], 'Rn', 'Radon';
        223,     [0.26 0.00 0.40], 'Fr', 'Francium';
        226,     [0.00 0.49 0.00], 'Ra', 'Radium';
        227,     [0.44 0.67 0.98], 'Ac', 'Actinium';
        232.04,  [0.00 0.73 1.00], 'Th', 'Thorium';
        231.04,  [0.00 0.63 1.00], 'Pa', 'Protactinium';
        238.03,  [0.00 0.56 1.00], 'U',  'Uranium';
        237,     [0.00 0.50 1.00], 'Np', 'Neptunium';
        244,     [0.00 0.42 1.00], 'Pu', 'Plutonium';
        243,     [0.33 0.36 0.95], 'Am', 'Americium';
        247,     [0.47 0.36 0.89], 'Cm', 'Curium';
        247.1,   [0.54 0.31 0.89], 'Bk', 'Berkelium';
        251,     [0.63 0.21 0.83], 'Cf', 'Californium';
        252,     [0.70 0.12 0.83], 'Es', 'Einsteinium';
        257,     [0.70 0.12 0.73], 'Fm', 'Fermium';
        258,     [0.70 0.05 0.65], 'Md', 'Mendelevium';
        259,     [0.74 0.05 0.53], 'No', 'Nobelium';
        262,     [0.78 0.00 0.40], 'Lr', 'Lawrencium';
        267,     [0.80 0.00 0.35], 'Rf', 'Rutherfordium';
        270,     [0.82 0.00 0.31], 'Db', 'Dubnium';
        271,     [0.85 0.00 0.27], 'Sg', 'Seaborgium';
        270.1,   [0.88 0.00 0.22], 'Bh', 'Bohrium';
        277,     [0.90 0.00 0.18], 'Hs', 'Hassium';
        278,     [0.92 0.00 0.15], 'Mt', 'Meitnerium';
        281,     [0.94 0.00 0.13], 'Ds', 'Darmstadtium';
        282,     [0.90 0.00 0.20], 'Rg', 'Roentgenium';
        285,     [0.55 0.55 0.55], 'Cn', 'Copernicium';
        286,     [0.60 0.30 0.30], 'Nh', 'Nihonium';
        289,     [0.60 0.30 0.35], 'Fl', 'Flerovium';
        290,     [0.60 0.30 0.40], 'Mc', 'Moscovium';
        293,     [0.60 0.30 0.45], 'Lv', 'Livermorium';
        294,     [0.60 0.30 0.50], 'Ts', 'Tennessine';
        294.1,   [0.30 0.60 0.65], 'Og', 'Oganesson'
    };

    tableMasses = cell2mat(elements(:,1));

    uniqueMasses = unique(masses);
    palette = lines(length(uniqueMasses));

    for i = 1:Natoms
        m = masses(i);

        diffs = abs(tableMasses - m);
        [minDiff, idx] = min(diffs);
        refMass = tableMasses(idx);

        absTol = 0.3;
        relTol = 0.01;
        tol = max(absTol, relTol * refMass);

        if minDiff <= tol
            colors(i, :) = elements{idx, 2};
            elementNames{i} = sprintf('%s (%s)', elements{idx, 4}, elements{idx, 3});
        else
            idx2 = find(uniqueMasses == m, 1);
            colors(i, :) = palette(idx2, :);
            elementNames{i} = sprintf('Generic Atom (m = %.2f)', m);
        end
    end
end
