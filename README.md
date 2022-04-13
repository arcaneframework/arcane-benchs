# Benchs for Arcane

This repository contains benchs for evaluating performance of some
message passing patterns.

To get sources files for the benchs:

~~~{.sh}
git clone https://github.com/arcaneframework/arcane-benchs
~~~

## Compilation

### Arcane

You need to compile and install
[Arcane](https://github.com/arcaneframework/framework) in Release
configuration following the instructions in the github
repository. Version 3.5.7+ of Arcane is required. We
suppose Arcane will be installed in the directory
`${ARCANE_INSTALL_PATH}`. Instead of using a separate distribution
of Arcane You can use the sources of Arcane provided in this
repository if you use git submodules:

~~~{.sh}
git submodule update --init --recursive
cmake -S framework -B /tmp/arcane_build_path -DCMAKE_INSTALL_PREFIX=${ARCANE_INSTALL_PATH}
cmake --build /tmp/arcane_build_path --target install
~~~

### Bench

A recent (3.19+) version of [CMake](https://cmake.org/) is needed to compile. In-source
compilation is forbidden. You can configure and compile the benchs using the following command:

~~~{.sh}
# To configure if you are in the source directory
cmake -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_PATH} -B build_bench
# To compile
cmake --build build_bench
~~~

## Benchs

Two benchs are availables:

- MicroHydro for load balance
- MaHyCo for synchronisation

The benchs are launched by shell scripts generated in the build
directory. You can specifiy the environment variable ` MPI_ARGS` if you
want to add some arguments to the MPI launch command.

It is important to note that by default only the processus with rank 0
will display information in the listing. You may set the environment
variable `ARCANE_PARALLEL_OUTPUT` to `1` to allow all other ranks to print
their corresponding listing in a file `output%R` where `%R` is the rank.

### MicroHydro

This bench is used to evaluate cell migration during
load-balancing. During cell migration, a subdomain send cells and its
associated variables to other subdomains. All the informations for a
subdomain are sent in one message. In this benchmark, the size of one
message is in the range 10Mo to 100Mo.

Two implementation are available:

- The default implementation use MPI_Irecv/MPI_Isend/MPI_Wait.
- A collective implementation using 'MPI_Alltoallv'. This
implementation is available if you set the
`ARCANE_MESH_EXCHANGE_USE_COLLECTIVE` environment variable to
`TRUE`.

The important information is the time spent to send cell
information. The following line in the listing gives the information:

~~~{txt]
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end)
total_time=0.0058741569519043
~~~

You can launch the test and get the results with the following
commands:

~~~{sh}
sh launch_micro_hydro_lb8.1.sh > result1.txt
grep 'total_time=' result1.txt | grep Cell
~~~

The results will be similar to this listing:

~~~{txt}
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.0058741569519043
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00423097610473633
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00581598281860352
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00305461883544922
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00515961647033691
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00465488433837891
*I-Mesh       ParallelExchanger Cell: ProcessExchange (end) total_time=0.00334906578063965
~~~

There is two tests using 512 sub-domains:

~~~{txt}
launch_micro_hydro_lb512.1.sh
launch_micro_hydro_lb512.2.sh
~~~

If you want to test if anything is OK you can launch smaller tests:

~~~{txt}
launch_micro_hydro_lb32.1.sh
launch_micro_hydro_lb8.1.sh
~~~

### MaHyCo

This bench is used to evaluate weak scaling of synchronisations. In
this configuration, a subdomain is the owner of approximately 18000
cells with 3 layers of ghost cells. The scripts are provided to launch
the test from 128 subdomains to 2048. The important information is the
total time of the execution loop. If the scalability of
synchronisations is perfect, that time should be constant. You can
launch a test and get the results with the following commands:

~~~{sh}
sh launch_mahyco_512.sh > result2.txt
grep '  loop:' result2.txt
~~~

The output is:

~~~{txt}
*I-Internal   TotalReel = 15.5567035675049 secondes (init: 0.841777086257935  loop: 14.7149264812469 )
~~~
