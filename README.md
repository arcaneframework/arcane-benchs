Benchs for Arcane
=================

To get sources files for the benchs:

~~~{.sh}
git clone --recurse-submodules /path/to/git
~~~

or

~~~{.sh}
git clone /path/to/git
cd framework && git submodule update --init --recursive
~~~

Compilation
-----------

You need to compile and install
[Arcane](https://github.com/arcaneframework/framework) in Release
configuration following the instructions in the github repository. We
suppose Arcane is installed in the directory `${ARCANE_INSTALL_PATH}`

You need to compile MaHyCo using the following command:

~~~{.sh}
# To configure
cmake -S mahyco -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_PATH} -B build_mahyco
# To compile
cmake --build build_mahyco
~~~

To launch MicroHydro tests:

...

To launch MaHyCo tests:

...
