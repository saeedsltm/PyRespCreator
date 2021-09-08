# PyRespCreator

Create your response file from scratch and save it in varoius standard formats!

- Create a new folder in "instruments" directory, then put pole, zero and other required info into the same files names in created folder,
- Update "request.json" file according to your new station info,
- Run the "PyCloneFromPZ.py" script to make new RESPONSE files,

# Note

- Inside "PyCloneFromPZ.py" script file, there are two helper functions for creating "request.json" (dlsv2json) and "STATION0.HYP" (dlsv2sta0) from dataless files,
- The output RESPONSE files will be stored in dataless format. If you need to create RESP format file too, you need to uncomment the lines "163-166" from "PyCloneFromPZ.py" script file,
- If you need to have SCML RESPONSE format (SeisComP xml format), uncomment the lines "169-172" from "PyCloneFromPZ.py" script file,
