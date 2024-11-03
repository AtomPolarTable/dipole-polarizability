# Generates API html for horton_part while ignoring the test and data folders.
# Stores it in pyapi/

sphinx-apidoc -a -o pyapi/ ../src/pydipole ../src/pydipole/tests/ ../src/pydipole/test/ ../src/pydipole/data/ --separate
