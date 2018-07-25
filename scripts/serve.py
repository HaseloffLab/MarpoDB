from marpodb.server import app

if __name__ == '__main__':
	print app.instance_path
	app.run(debug=True,host='0.0.0.0', port = 8084)