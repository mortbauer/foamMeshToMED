This is a OpenFOAM mesh to MED mesh (Salome) translator
#######################################################

In order to use the mesh of the excellent `cfMesh`_ produced mesh with
`code_saturne`_ I needed a mesh converter. Actually OpenFOAM can write the mesh to
Ensight Gold format which is readable by code_saturne, but the patches weren't
translated. This translater is able to do it.

Instructions
************

Download it best into `$FOAM_UTILITIES/mesh/converter/foamToMED`, but probably
the path doesn't matter, didn't check it. Then compile it in the usual OpenFOAM
way, by first sourcing your shell environment file::

    source $HOME/OpenFOAM/OpenFOAM-X.Y.Z/etc/bashrc

and then cd into the source path of the downloaded meshconverter and type::

    wmake

this should compile the conervert and place the executable "foamMeshToMED" into
your `$FOAM_APPBIN` dir.

If it doesn't work it probably doesn't find the header or the library of `med`
which is a prerequisite_ for this converter.

.. _prerequisite:

Prerquisites
************

Additionally to a working OpenFOAM setup you need also a version of med which
you probably have anyways, the source can also be found at the `salome
download`_.

.. _cfMesh: http://www.c-fields.com/solutions/products/meshing
.. _salome download: http://www.salome-platform.org/downloads/current-version
.. _code_saturne: http://code-saturne.org/cms/
