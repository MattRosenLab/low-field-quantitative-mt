# TNMRfileloadScripts

Helper functions for parsing Tecmag `.tnt` data files. Used internally by the processing scripts in `MATLAB_processingScripts/`; these do not need to be called directly.

| Script | Description |
|---|---|
| `Read_Tecmag_DK.m` | Loads raw FID data and header from a Tecmag file (interactive file picker). |
| `Read_Tecmag.m` / `Read_Tecmag_hdr.m` | Lower-level file readers for Tecmag binary data and header blocks. |
| `Read_Tecmag_Header.m` / `Read_Tecmag_Header_Table.m` / `Read_Tecmag_Header_Var.m` | Routines for extracting specific variables and parameter tables from the Tecmag file header. |
| `Read_Raw_Field.m` / `Parse_Tecmag_Variable.m` | Primitive parsing utilities used by the header readers. |
