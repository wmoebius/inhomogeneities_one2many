To run everything on the FAS RC cluster odyssey, one needs to:

* copy this directory somewhere and rename it to ManyFeatures_template
* incorporate fmmfronthelper.py into ManyFeatures_template

* cp/mv createallfoldersrandnumseeds.sh to the parent directory and execute the script from there,
  argument will be appended to names of directories
* if not automatically submitted by createallfoldersrandnumseeds.sh,
  go into each directory and run submit.sh
  or copy submitonebyone.sh to parent directory and automate this process somewhat

* cp/mv analyseall.sh to the parent directory and execute the script from there
