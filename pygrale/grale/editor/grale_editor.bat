@echo off
for /f %%i in ('python -c "import os;import grale.editor;print(os.path.dirname(grale.editor.__file__))" ') do set X=%%i
python "%X%\mainwindow.py" %*
