This repository contains the custom simulation pipeline and metric calculation scripts for the paper: Ultimate Intrinsic Temporal Performance in Magnetic Resonance Imaging Using a Realistic Human Head Model (Currently under review / Published in IEEE AWPL).

📌 Overview
This project provides a full-wave electromagnetic (EM) simulation framework to evaluate the fundamental temporal stability limits of RF receive coils in MRI (at 1.5 T, 3 T, and 7 T) under head motion.

Unlike conventional image-domain motion correction, this pipeline isolates the pure electrodynamic effects (coil-sample interaction) by calculating the intrinsic sensitivity variability ($\lambda^*$) and intrinsic thermal-noise variability ($\alpha^*$) using a generalized Huygens-surface approach.

⚙️ Prerequisites & Dependencies
To run the scripts, you will need:

MATLAB (Tested on R2022a or newer)
