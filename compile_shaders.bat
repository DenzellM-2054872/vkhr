@REM forfiles /P share\shaders /S /C "glslc @RELPATH/@file -o @RELPATH/@file.spv"
For /R .\share\shaders %%i in ("*.frag")do glslc %%i -o %%i.spv
For /R .\share\shaders %%i in ("*.geom")do glslc %%i -o %%i.spv
For /R .\share\shaders %%i in ("*.vert")do glslc %%i -o %%i.spv
For /R .\share\shaders %%i in ("*.comp")do glslc %%i -o %%i.spv