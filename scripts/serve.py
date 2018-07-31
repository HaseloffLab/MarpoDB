from marpodb.server import app

if __name__ == '__main__':
	print app.instance_path
	app.config['INTERPRO_PATH'] = "/Users/md/MarpoDB4/data/marpodb4/marpodb4_cam_tak_2018_06_27_interpro"
	app.run(debug=True,host='0.0.0.0', port = 8084)