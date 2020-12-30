# Pole_Figure_Inversion
This repository contains ADMM based inversion of Crystallographic Pole Figures

### Requirement
1.	Subroutines in this repository can be run using either Matlab or Octave 
2.	Install packages statistics & image if you opt to run calculations using Octave

### Main Files
1. main_Sf.m inverts simulated polefigures from simulated Santa Fe ODF
2. main_Cu.m inverts experimental polefigures from Cu sample

### Results
<img src="/results/Sf_complete/TV/pf_err_figs/noiselevel2.png" alt="Input PF fig"/>
<img src="/results/Sf_complete/TV/pf_rec_px_figs/noiselevel2.png" alt="Output Pf fig"/>

### Citation
    @article{po5145,
        author={S. Singh and P. Kc and S. Kashyap and M. D. Graef}, 
        journal={Journal of Applied Crystallography}, 
        title={An iterative reconstruction algorithm for pole figure inversion using total variation regularization}, 
        year={2019}, 
        volume={52}, 
        number={}, 
        pages={}, 
        keywords={pole figure inversion; orientation distribution function; total variation regularization; texture}, 
        doi={https://doi.org/10.1107/S1600576719013529},
        ISSN={1600-5767}, 
        month={Oct}
    }

### Contact
Prabhat KC
prabhat.kc077@gmail.com<br>
Saransh Singh
saransh1@llnl.gov<br>
Marc De Graef
mdg@andrew.cmu.edu
