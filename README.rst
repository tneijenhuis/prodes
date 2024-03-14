Prodes
===========
Prodes is calculates features from 3D protein structures with a general focus on the protein surface.

Requirements
~~~~~~~~~~~~~~

Prodes requires Python 3.7 or higher. Additional requirements are specified in the `requirements.txt` file and automatically installec when installing with pip.

Installation
~~~~~~~~~~~~~

Installation of ProDes is currently only available using the source code by cloning repository and running

``pip install .``

or if an editable version is required by running

``pip install -e .``

in the source directory.

It is advised to perform the installation in a virtual enviroment like 
`venv <https://docs.python.org/3/library/venv.html#:~:text=A%20virtual%20environment%20is%20a,part%20of%20your%20operating%20system>`_.

Getting started
~~~~~~~~~~~~~~~~~~
ProDes can be used as a module or via the installed script by running

``python -m prodes [infile] [outfile] [options]``

The infile requires to be a PDB file and the outfile a csv.

``python -m prodes --help``

will return the list of options available, which include:

* Specific pH
* Specific surface probe size
* Costum pKa values

Prodes as a module examples
-------------------------------
To run Prodes and perform a similar action as the comandline implementation simply use

.. code-block:: python
    
    import prodes
    prodes.run("./tests/data/1GDW.pdb", "example.csv")
    
Alternatively, the surface area of a protein can be calculated by running a script similar to 

.. code-block:: python

    from prodes.io.parser import PDBparser
    from prodes.calculations import grid_wizard, sasa
    
    structure = PDBparser().parse("./tests/data/1GDW.pdb")
    grid = grid_wizard.Grid(10)
    grid.construct_cells(structure.heavy_atoms)
    grid.fill_cells(structure.heavy_atoms)
    
    sasa.shrake_rupley(grid)
    
    print(structure.surface_area())
    

How to cite
~~~~~~~~~~~~~~~

If this package is useful for you, please cite Prodes in your publications:

Neijenhuis, T., Le Bussy, O., Geldhof, G., Klijn, M. E., & Ottens, M. (2024). Predicting protein retention in ion-exchange chromatography using an open source QSPR workflow. Biotechnology Journal, 19, e2300708. https://doi.org/10.1002/biot.202300708