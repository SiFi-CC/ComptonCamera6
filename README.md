# ComptonCamera6

Prerequisites
------------------------------------------------
* Required ROOT version: 6.13/01
* Required MathMore (optional ROOT module that requires GSL)
* Doxygen version 1.8.15
* Compiler supporting c++14

Sources
-------

Sources repository:
```
https://github.com/SiFi-CC/ComptonCamera6/
```

To get sources run:
```
git clone https://github.com/SiFi-CC/ComptonCamera6/
```

Building and installation
-------------------------
```
in the ComptonCamera6 directory run:
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/dir
make
make install
```

Documentation
-------------
By default documentation is built and installed.
Add following options to cmake:

* `-DBUILD_DOC=Off`  -- to prevent building of documentation
* `-DINSTALL_DOC=Off` -- will prevent installation of documentation to final location


Installation directories
------------------------

Following directory structure will be created
```
<location>/bin
          /lib or /lib64            # depends on system
          /include
          /share/ComptonCamera6/
          /share/ComptonCamera6/doc # if doc installed
```

Scripts and macros
------------------

Following scripts and macros will be also generated and installed into `<location>/share/ComptonCamera6/`:

* `profile.sh` - set up all paths to run the CC6 executables
* `rootlogon.C` - startup macro for ROOT to load CC6 libraries

To make ROOT use this `rootlogon.C` set following value in your `~/.rootrc`:

```
Rint.Logon:              $(ROOTLOGON)
```
`profile.sh` adds location of CC6 lib and bin directories to system variables 'LD_LIBRARY_PATH' and 'PATH', respectively. If you wish to have it done automatically, insert the following line in your `~\.bashrc` file:
```
source <location>/share/ComptonCamera6/profile.sh
```
