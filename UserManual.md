# GMBCG cluster finder user manual #

Before using this package, you need to have IDL installed and a sdssidl package (http://code.google.com/p/sdssidl/) needs to be installed too.

The following tells how to run the gmbcg cluster finder. Jiangang Hao 12/14/10

0. The input catalog should be in fit file format. Exactly specified as here:

  * ID: A unique id # for that galaxy
  * OMAG: The observed DES magnitudes. Should be an array of g,r,i,z,Y magnitudes
  * OMAGERR: Estimated photometric errors for each band.
  * RA: Galaxy position.
  * DEC: Galaxy position.
  * PHOTOZ: Estimated photo-z based using other methods

#---------the following is the directions about how to run it--------------

1. in the config file, change the input directory and output directory. Please do not add any space when you do the changes. Just replace the directory part.

2. Specify a radius in the config file. This is the searching aperture, in unit of Mpc. I put 1.0 Mpc for example.

3. In command line, type: idl rungmbcg

4. It is done!