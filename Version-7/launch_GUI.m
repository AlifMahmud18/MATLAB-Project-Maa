%% launch_GUI.m
% Run this to open the Crystal Vibration Explorer GUI.
%
% This is the interactive deliverable from the project proposal (the
% "Interactive MATLAB-Based Crystalline Structure Viewer and Vibrational
% Mode Animator"). MainProject.m remains a separate script - it validates
% the underlying physics/math engine (Jacobi solver vs. eig(), theGFG
% Maxwell rigidity check, etc.) and is worth keeping and running too,
% but this GUI is what you actually interact with and demo.

app = CrystalVibrationApp; %#ok<NASGU>
