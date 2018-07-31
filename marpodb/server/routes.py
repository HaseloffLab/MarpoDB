from ..core import *
from ..server import app

from .user import *
from .backend import *

import os

from flask import redirect, render_template, url_for, request, abort, safe_join, flash

@app.route('/')
def index():
	return redirect(url_for('query'))

@app.route('/query', methods=['GET', 'POST'])
def query():
	loginForm = LoginForm()
	if request.method == 'POST':
		scope = request.form.getlist("searchScope")
		scope = '|'.join(scope)
		term = request.form["term"]

		dataset = '|'.join( request.form.getlist("dataset") )

		if scope !="" and term !="":
			return redirect(url_for('result', term=request.form['term'], scope = scope, dataset = dataset))
		else:
			flash('Empty search term or scope')

	return render_template("query.html", loginForm = loginForm)

@app.route('/results', methods = ['GET'])
def result():
	nHits = 5
	displayColumns = ['origin', 'description', 'eVal']

	session = marpodb.Session()

	term = request.args.get('term','')
	scope = request.args.get('scope','').split('|')
	dataset = request.args.get('dataset', '').split('|')

	if not (term and scope and dataset):
		 flash('Search term, scope and dataset should not be empty!')
		 return redirect(url_for('query'))

	table = processQuery(session, scope, term, dataset, displayColumns, nHits)

	session.close()

	if len(table['data']) != 0:
		return render_template('results.html', table = table, title='Query results')
	else:
		flash("No entries found")
		return redirect(url_for('index'))

@app.route('/interprohtml/<cdsdbid>')
def interproHTML(cdsdbid):
	return send_from_directory(app.config['INTERPRO_PATH'], "{0}.html".format(cdsdbid))

@app.route('/details')
def details():
	dbid = request.args.get('dbid','')

	if dbid:
		dbidParts = dbid.split('.')
		if len(dbidParts) == 3:
			tCls = getClassByTablename( dbidParts[1] )
			if tCls:
				session = marpodb.Session()

				gene	= None
				cds		= None

				if tCls == Gene:
					gene = session.query(Gene).filter(Gene.dbid == dbid).first()
				
				else:
					part = session.query(tCls).filter(tCls.dbid == dbid).first()
					if part:
						gene = part.gene[0]

				if gene:
					cds = gene.cds
					
					geneDetails = getGeneDetails(session, gene.dbid)
					annotation = getCDSDetails(session, cds.dbid)

					interproAnnotation = ""
					interproFileName = safe_join(app.config["INTERPRO_PATH"], "{0}.html".format(cds.dbid))

					if os.path.exists(interproFileName):
						interproFile = open(interproFileName)
						interproAnnotation = parseInterProHTML( interproFile.read() )

					print annotation["dbxref"]
					
					return render_template('details.html', geneDBID = gene.dbid, cdsDBID = cds.dbid, geneCoordinates = geneDetails["parts"], seq = geneDetails["seq"],  title = "Details for {0}".format(gene.dbid), blastp=annotation["blastp"], alias = annotation["dbxref"], stared = False, interpro = interproAnnotation)

@app.route('/blast', methods=['GET', 'POST'] )
def blast():
	return render_template('blast.html', title='BLAST to MarpoDB')