# gpmaterials
GPMaterials implements a rich set of finite element formulations for finite strain hyperelasticity, damage, multi-scale solid mechanics in both continuum and discrete level (truss-like network of fibres). This library was developed during my master and PhD studies at LNCC, under guidance of Prof. Pablo Blanco. Particularly this version is more related with my PhD and indeed it is the last stable version from late 2018, before starting writing up the thesis. 

The library should be compiled in a .so file and then used by SolverGP (see below reference). SolverGP itself is not opensource. Consider using "Piola", developed by myself mimeting the interface of SolverGP, using C++. If you didn't manage to couple them please drop me a e-mail at f.rocha.felipe@gmail.com. 

Instructions

i) cd ./build

ii) cmake .

iii) make 

iv) check if libgpmaterials.so and/or libgpmaterials_static.a have been created. 

PS : check if you have a valid fortran compiler, I'm using gfortran of the gcc version 7.5.0.


SolverGP reference : Urquiza, S. A. and Venère, M. J. 2002. An application framework architecture for fem and
other related solvers. In Mecánica Computacional, S.R.Idelsohn, V.E.Sonzogni, and A.Cardona,
Eds. Vol. 22. Santa Fé - Paran´a , Argentina, 3099–3109.


You can find my thesis at
https://sites.google.com/view/feliperocha/phdmaster-dissertations?authuser=0
and some pieces of my research at
https://sites.google.com/view/feliperocha/research?authuser=0

Please consider citing at least one of the following documents this library has been useful to your research:

@article{Rocha2018,
	author = {Rocha, Felipe Figueredo and Blanco, Pablo Javier and S{\'{a}}nchez, Pablo Javier and Feij{\'{o}}o, Ra{\'{u}}l Antonino},
	doi = {10.1016/j.cma.2018.06.031},
	file = {:home/felipe/thesis/bibliography/Group/MTS{\_}2018{\_}Rocha{\_}MultiscaleFibres.pdf:pdf},
	issn = {00457825},
	journal = {Computer Methods in Applied Mechanics and Engineering},
	keywords = {biological tissues,fibre network,multi-scale modelling,non-affinity,representative volume element,virtual power},
	pages = {740--787},
	publisher = {Elsevier B.V.},
	title = {{Multi-scale modelling of arterial tissue: Linking networks of fibres to continua}},
	url = {https://linkinghub.elsevier.com/retrieve/pii/S0045782518303281},
	volume = {341},
	year = {2018}
}

@article{Rocha2021,
title = {Damage-driven strain localisation in networks of fibres: A computational homogenisation approach},
journal = {Computers & Structures},
volume = {255},
pages = {106635},
year = {2021},
issn = {0045-7949},
doi = {https://doi.org/10.1016/j.compstruc.2021.106635},
url = {https://www.sciencedirect.com/science/article/pii/S0045794921001577},
author = {Felipe Figueredo Rocha and Pablo Javier Blanco and Pablo Javier Sánchez and Eduardo {de Souza Neto} and Raúl Antonino Feijóo},
}


@MastersThesis{mythesisMaster,
author = {F.F Rocha},
title  = {Basics Aspects of Multi-Scale Modelling in Biological Tissues (In Portuguese)},
school = {National Laboratory of Scientific Computing},
year   = {2014}
}

@PhDThesis{mythesisPhD,
	author = {F.F Rocha},
	title  = {Multiscale Modelling of Fibrous Materials: from
	the elastic regime to failure detection in soft tissues},
	school = {National Laboratory of Scientific Computing},
	year   = {2019}
}

@article{Blanco_Rocha_2019,
	author = {Blanco, Pablo Javier and S{\'{a}}nchez, Pablo Javier and Rocha, Felipe Figueredo and Toro, Sebastian  and Feij{\'{o}}o, Ra{\'{u}}l Antonino},
	journal = {Computer Methods in Applied Mechanics and Engineering},
	publisher = {Elsevier B.V.},
	title = {{A consistent multiscale mechanical formulation for media with randomly distributed voids (in submission process)}},
	year = {2019}
}

@inproceedings{Rocha2015,
	author =    {F.F. Rocha and P.J. Blanco and R.A. Feij{\'{o}}o   and P.J. S{\'{a}}nchez and A.E. Huespe},
	title =     {A multi-scale approach to model arterial tissue},
	publisher = {XXXVI Congresso Ibero-Latino-Americano de M{\'e}todos Computacionais em Engenharia},
	booktitle = {Ibero-Latin American Congress on Computational Methods in Engineering (CILAMCE)},
	address =   {Rio de Janeiro},
	year =      {2015},
	volume =    {},
	pages =     {},
	note =      {}
}

