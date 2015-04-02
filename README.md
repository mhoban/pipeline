#Analysis pipeline
Analysis pipeline for next-generation (Illumina) sequence data consisting of pooled samples

###Usage notes
Before doing anything with your data or any of the scripts in here, modify the pipeline.sh script
so that $PIPEDIR points to the appropriate location. Then **source** the pipeline.sh
script located in the "root" of the pipeline tree:

**$ . pipeline.sh**  
or  
**$ source pipeline.sh**  

###Installation notes
There is a script in the src/ directory that gets source packages you might need. Run the script and attempt
to compile and install all the stuff. In general, it's a good idea to ./configure each source package with --prefix=$PIPEDIR
To get cutadapt to install, first install Cython, the following will get you there, for each package:

**$ cd cython-src**  
**$ python setup.py install --user**
**$ cd ../cutadapt-src**  
**$ python setup.py install --user**

To get fastx_toolkit to install properly, you'll need to first build libgtextutils. Again, configure it with a --prefix=$PIPEDIR.
When you go to ./configure fastx_toolkit, make sure you fist set the environment variable PKG_CONFIG_PATH to include the location
of the pkgconfig file spat out by building libgtextutils (usually $PREFIX/lib/pkgconfig).

There is some other stuff you may need that's not automatically downloaded:

    -Latest zlib (needed by PEAR)
        -This will need to be compiled and installed to the same prefix as PEAR
         or else PEAR won't compile properly. 
         .i.e. the ./configure script needs the same --prefix argument to both
         zlib and PEAR
        -You'll know you need this if PEAR compilation keeps breaking with unresolved _gzoffset symbols
    -Latest autoconf (also needed by PEAR)
        -Build and install this before running the autogen.sh script
