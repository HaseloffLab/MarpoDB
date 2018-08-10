from marpodb.server import app

if __name__ == '__main__':
	app.config['INTERPRO_PATH'] = ""
	app.config['BLASTDB_PATH']	= ""
	app.config['BLAST_PATH']	= ""
	app.config["FASTA_PATH"]	= ""
	app.config["TEMP_PATH"]	= ""
	app.run(debug=True,host='0.0.0.0', port = 8084)