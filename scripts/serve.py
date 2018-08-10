from marpodb.server import app

if __name__ == '__main__':
	print app.instance_path
	app.config['INTERPRO_PATH'] = "/Users/md/MarpoDB4/data/marpodb4/marpodb4_cam_tak_2018_06_27_interpro"
	app.config['BLASTDB_PATH']	= "/Users/md/MarpoDB4/data/marpodb4/blastdb"
	app.config['BLAST_PATH']	= "/usr/local/ncbi/blast/bin/"
	app.config["FASTA_PATH"]	= "/Users/md/MarpoDB4/data/marpodb4/fasta"
	app.config["TEMP_PATH"]	= "/Users/md/MarpoDB4/data/temp"
	app.run(debug=True,host='0.0.0.0', port = 8084)