# ComptonCamera6

Prerequisites
------------------------------------------------
* Required ROOT version: 6.13/01
* Doxygen version 1.8.15

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
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_DIR=/path/to/install/dir
make
make install
```

Documentation
-------------

Add following options to cmake:

* `-DBUILD_DOC=On`  -- to build documentation, will be found in build/doc
* `-DINSTALL_DOC=On` -- will install documentation to final location


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

Following scripts and macros will be also generated:

* `profile.sh` - set up all paths to run the CC6 executables
* `rootlogin.C` - startup macro for ROOT to load CC6 libraries

To make ROOT using the `rootlogin.C` set following value in your `~/.rootrc`:

```
Rint.Logon:              $(ROOTLOGON)
```