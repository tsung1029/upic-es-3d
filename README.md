# upic-es-3d
a 3D electrostatic code for LPI and other fundamental plasma physics phenomenon

UPIC-ES-3D is based on Viktor Decyk's UPIC framework, and use OSIRIS's diagnostic routines to provide a uniform post-processing environment.  

to compile, type

make -f (makefile) production

to produce the production version of the code, to compile the production code along with generic upic, type

make -f (makefile) all


Currently available makefile include:

macosx-gnu.make  --> for OS X using the GNU compiler.

cori-haswell.make --> for Cori HASWELL CPU's.  (using intel compiler).  It requires the following modules

<ul>
    <li>PrgEnv-intel 
    <li>cray-hdf5-parallel 
    <li>craype-haswell 
</ul>
