I/O
===

Gridded output
--------------

Variables defined on the fixed Eulerian grid can be output to netcdf files. Within the code grid variables are defined as two dimensional Kokkos arrays:

.. code::

   Kokkos::View<double**> griddedVariableExample;

The variable must be registered for potential output:

.. code::

   grid->register_field_for_write(&griddedVariableExample);

To output the variable to a netcdf file during the simulation a gridded output stream should be defined in the configs file:

.. code::

  <!-- Gridded output -->
  <ConfigGroup name="griddedOutput">
    <Array name="griddedOutputStreams">
      <GriddedOutputStream name="griddedOutput">
	<Config name="griddedWriteFilenameTemplate" type="String" value="gridded_out.$Y-$M-$D_$h.nc" units="none"
		description="Name of NetCDF particle output filename template."
		possible_values='File name, or "none"'/>
	<Config name="griddedWriteOutputDirectory" type="String" value="./output/" units="none"
		description="Name of NetCDF particle output directory."
		possible_values='Directory name, or "none"'/>
	<Config name="griddedWriteInterval" type="String" value="HOURS:6" units="none"
		description="Interval for particle output."
		possible_values="Time interval string format."/>
	<Config name="griddedWriteClobber" type="Boolean" value="true" units="none"
		description="If true clobber particle output files."
		possible_values="'True', or 'False'"/>
	<Array name="OutputFields">
	  <OutputField>griddedVariableExample</OutputField>
	</Array>
      </GriddedOutputStream>
    </Array>
  </ConfigGroup>

Several options controlling the output are defined:

+------------------------------+----------------------------------------------------------+
| Option name                  | Description                                              |
+==============================+==========================================================+
| griddedWriteFilenameTemplate | Filename template for the output file                    |
+------------------------------+----------------------------------------------------------+
| griddedWriteOutputDirectory  | Directory to write the output file to                    |
+------------------------------+----------------------------------------------------------+
| griddedWriteInterval         | Time interval for writing the output gridded data        |
+------------------------------+----------------------------------------------------------+
| griddedWriteClobber          | Whether to over write the file if it already exists      |
+------------------------------+----------------------------------------------------------+
