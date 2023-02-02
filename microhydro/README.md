# MicroHydro

This bench is used to test the basics functionalities of Arcane and the load-balancing

The source files for this bench are located in the directory [`../framework/arcane/extras/benchs/microhydro`](https://github.com/arcaneframework/framework/tree/main/arcane/extras/benchs/microhydro)


## Output results in CSV file

Like QAMA, MicroHydro uses SimpleCsvOutput to provide an easy method for outputting results
to CSV files.

In the .arc file, you have the following options:
- `tableDir`  : directory for the output CSV files
- `tableName` : name of the output CSV files

In these file names, you can use `@proc_id@` (which will be replaced by the rank of the process)
and/or `@n_proc@` (which will be replaced by the number of processes) for customizing the file names.

Example:
```xml
<tableDir>ExampleFull</tableDir>
<tableName>Results_P@proc_id@</tableName>
```
With this example, the CSV files will be located in `./output/csv/ExampleFull/Results_P0.csv`.


## Compare results with CSV reference file

MicroHydro uses SimpleCsvComparator to provide a method for comparing the current
results with a reference result. Reference files can be created by you,
or you can use the reference files available in this repository.

You have one option in .arc: `check-numerical-result`.
This option determines whether or not to check the results.
It uses SimpleCsvOutput options to determine the exact path of the reference files.

For example:
```xml
<check-numerical-result>true</check-numerical-result>
<tableDir>ExampleFull</tableDir>
<tableName>Results_P@proc_id@</tableName>
```
The path determined by SimpleCsvComparator is :
`./output/csv_refs/ExampleFull/Results_P0.csv`

To write or overwrite reference files, command line arguments must be used.
___
In MicroHydro, you have two command line arguments:
- `ReferenceDirectory`
- `OverwriteReference`

The `ReferenceDirectory` argument specifies the directory that contains all the reference files.

Note: `ReferenceDirectory` overrrides `<check-numerical-result>` option.
If `ReferenceDirectory` if set, `<check-numerical-result>` option is considered
as `true`.

___

For example:

`ReferenceDirectory="/tmp/my_tmp_refs"`
```xml
<tableDir>Example123</tableDir>
<tableName>Results456_P@proc_id@</tableName>
```
The path determined by SimpleCsvComparator is:
`tmp/my_tmp_refs/Example123/Results456_P0.csv`

___


`ReferenceDirectory="/tmp/my_tmp_refs"`
```xml
<!-- <tableDir></tableDir> -->
<tableName>Results456_P@proc_id@</tableName>
```
The path determined by SimpleCsvComparator is :
`tmp/my_tmp_refs/Results456_P0.csv`

___

`ReferenceDirectory="/tmp/my_tmp_refs"`
```xml
<!-- <tableDir></tableDir> -->
<!-- <tableName></tableName> -->
```
The path determined by SimpleCsvComparator is :
`tmp/my_tmp_refs/Table_P0.csv`
___

```xml
<check-numerical-result>false</check-numerical-result>
```
No comparison will be made.

___

To write or overwrite reference files, use `OverwriteReference`
argument.

For example:

`ReferenceDirectory="/tmp/my_tmp_refs"`
`OverwriteReference=true`
```xml
<tableDir>Example123</tableDir>
<tableName>Results456_P@proc_id@</tableName>
```
The path determined by SimpleCsvComparator is:
`tmp/my_tmp_refs/Example123/Results456_P0.csv`

___

Note:
In order to write new reference files, the destination directory
must be created beforehand.

___

Examples of usage:

To write new reference files in the `/tmp/my_tmp_refs` directory:
```bash
mkdir -p "/tmp/my_tmp_refs"
./MicroHydro -A,MaxIteration=10,ReferenceDirectory="/tmp/my_tmp_refs",OverwriteReference=true ./src/microhydro/data/MicroHydro.1.1.arc
ls "/tmp/my_tmp_refs/MicroHydro.1.1/Results_P0.csv"
```

To verify results using the reference files from the repository:
```bash
./MicroHydro -A,MaxIteration=10,ReferenceDirectory="./src/microhydro/reference_files" ./src/microhydro/data/MicroHydro.1.1.arc
```
