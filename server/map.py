import requests
import json

def generateNewMap(markers):

	with open('static/json/map_style.json') as dataFile:    
	    data = json.load(dataFile)

	data["style"] = ( ("|".join( key + ':' + str(style[key]) for key in sorted(style.keys()) )) for style in data["style"])
	data["markers"] = "size:tiny|color:green|University of Cambridge|icon:http://s33.postimg.org/dz4w2oavj/maps_icon2.png"

	r = requests.get("http://maps.googleapis.com/maps/api/staticmap", params = data, stream=True)

	print r.url

	chunk_size = 10

	if r.status_code == 200:
		with open("map.png", 'wb') as fd:
			for chunk in r.iter_content(chunk_size):
				fd.write(chunk)