from flask import Flask
import subprocess
import shutil
import os
import glob
from hs_restclient import HydroShare, HydroShareAuthBasic
import urllib2
import zipfile
from time import sleep
import json

app = Flask(__name__)


@app.route('/id/<id>')
def runScript(id):
    auth = HydroShareAuthBasic(username='', password='')
    hs = HydroShare(auth=auth)

    hs.getResource(id, destination='/home/ubuntu/', unzip=True)
    # Locate the file with the .nam extension
    os.chdir('/home/ubuntu/' + str(id) +"/" + str(id) + '/data/contents')
    for file in glob.glob("*.nam"):
        filename = file

    fpath = '/home/ubuntu/' + str(id) +"/" + str(id) + '/data/contents/'+ filename.split(".")[0]+'.lst'

    # Run the model
    subprocess.call("sudo docker run --rm -v /home/ubuntu/" + str(id) +"/" + str(id) +
                    "/data/contents:/workspace mjstealey/docker-modflow mf2005 " + filename.split(".")[0], shell=True)

    try:
        hs.deleteResourceFile(id,filename.split(".")[0]+'.lst')
    except:
        pass
    #Upload to hydroshare
    hs.addResourceFile(id, fpath)
    #metatdata
#    metadata = {
#    "title": "Metadata",
#        "modeloutput": [
#        {"includes_output": True}
#        ]
#    }
#    science_metadata_json = hs.updateScienceMetadata(
#        str(id), metadata=metadata)
    return json.dumps(hs.getScienceMetadata(id))

if __name__ == "__main__":
    app.run(debug=True)
