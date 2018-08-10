from ..core import *
from ..server import app
from ..server import LoginForm, RegisterForm
from ..server import loginManager, User, StarGene, bcrypt, userDB

from .user import *
from .backend import *

from flask_login import logout_user, login_user, current_user, login_required

import os
import md5
import subprocess

from flask import redirect, render_template, make_response, url_for, request, abort, safe_join, flash, session

############################## USER/LOGIN ROUTES ##############################
@app.context_processor
def user_data():
	return {"user_data" : current_user.data}

@app.route('/register', methods=['GET', 'POST'])
def register():
	form = RegisterForm(request.form)
	if request.method == 'POST' and form.validate():
			if User.query.filter(User.username == form.username.data).first():
				form.username.errors.append('Username is already registred')
				return render_template('registration.html', title = 'Registration', form = form)

			if User.query.filter(User.email == form.email.data).first():
				form.email.errors.append('E-mail is already registred')
				return render_template('registration.html', title = 'Registration', form = form)

			newUser = User(username=form.username.data, email=form.email.data, password=bcrypt.generate_password_hash(form.password.data, 12), affiliation=form.affiliation.data, active = True)

			userDB.session.add(newUser)
			userDB.session.commit()
		
			try:
				generateNewMap(User, app.static_folder)
			except:
				pass

			flash('You have been successfully registred!')
			return redirect(url_for('map'))

	return render_template('registration.html', title = 'Registration', form = form)

@app.route('/login', methods=['POST'])
def login():

	form = LoginForm(request.form)
	user = User.query.filter(User.username == form.username.data).first()
	
	if user and bcrypt.check_password_hash(user.password, form.password.data):
		login_user(user)
		session.pop("starGenes", None)
	else:
		flash("Invalid username or/and password")
	return redirect(url_for('index'))

@app.route("/logout")
@login_required
def logout():
	logout_user()
	return redirect(url_for('index'))


@app.route('/changestarname')
def changeStarName():
	cdsdbid = request.args.get('cdsdbid','')
	newName = request.args.get('newname','')
	marpodbSession = marpodb.Session()

	cds = marpodbSession.query(CDS).filter(CDS.dbid == cdsdbid).first()
	if cds and current_user.is_authenticated:
		star = StarGene.query.filter(StarGene.cdsdbid == cdsdbid, StarGene.userid == current_user.id).first()
		star.genename = newName
		userDB.session.add(star)
		userDB.session.commit()
		return "OK", 200, {'Content-Type': 'text/plain'}
	abort(500)

@app.route('/stargene')
def starGene():
	cdsdbid = request.args.get('cdsdbid','')
	marpodbSession = marpodb.Session()

	cds = marpodbSession.query(CDS).filter(CDS.dbid == cdsdbid).first()
		
	mpdbSession = marpodb.Session()
	annotation = getCDSDetails(mpdbSession, cds.dbid)
	mpdbSession.close()

	geneName = cdsdbid

	if 'MarpolBase' in annotation["dbxref"]:
		geneName = annotation["dbxref"]["MarpolBase"][0]

	if cds:
		if current_user.is_authenticated:
			star = StarGene.query.filter(StarGene.cdsdbid == cdsdbid, StarGene.userid == current_user.id).first()

			if star:
				userDB.session.delete(star)
				userDB.session.commit()
			
			else:
				newStar = StarGene(current_user.id, cdsdbid, geneName)
				userDB.session.add(newStar)
				userDB.session.commit()
		else:
			starGenes = session.get("starGenes", {})
			
			if cdsdbid in starGenes:
				starGenes.pop(cdsdbid)
			else:
				starGenes[cdsdbid] = cdsdbid

			session["starGenes"] = starGenes

		return "OK", 200, {'Content-Type': 'text/plain'}
	else:
		 abort(500)

@app.route('/map')
def map():
	if not os.path.isfile( os.path.join( app.static_folder, "img/map.png" ) ):
		generateNewMap(User, app.static_folder)
		print "Map not found"
	return render_template('map.html', title='Community map')

@app.route('/top')
def top():

	marpodbSession = marpodb.Session()
	cdsTop = getTopGenes(marpodbSession, StarGene, 10)
	marpodbSession.close()

	return render_template('top.html', list = cdsTop, title='Top genes')

############################## CORE ROUTES ##############################
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

					if cds.dbid in current_user.data['starGenes']:
						titleEx = '<img src="static/img/star.png" onclick="starGene()" id="star_img"/>'
					else:
						titleEx = '<img src="static/img/star_na.png" onclick="starGene()" id="star_img"/>'
					
					return render_template('details.html', geneDBID = gene.dbid, cdsDBID = cds.dbid, geneCoordinates = geneDetails["parts"], seq = geneDetails["seq"],  title = "Details for {0}".format(gene.dbid), titleEx = titleEx, blastp=annotation["blastp"], alias = annotation["dbxref"], stared = False, interpro = interproAnnotation)
	abort(404)

@app.route('/blast', methods=['GET', 'POST'] )
def blast():
	if request.method == 'POST':
		query	= request.form["query"]
		evalue	= request.form["evalue"]
		program	= request.form["program"]
		matrix	= request.form["matrix"]
		perc	= request.form["perc"]
		dataset	= request.form["dataset"]

		print dataset

		if query == "":
			flash('Please provide a query sequence')

		if not (query and evalue and program and matrix and perc and dataset):
			flash('Request error')

		if program.startswith('blastn'):
			dbRoute = 'MarpoDB_{0}_Genes'.format(dataset)
			idType = 'Gene'
		else:
			dbRoute = 'MarpoDB_{0}_Proteins'.format(dataset)
			idType = 'CDS'

		dbRoute  = os.path.join( app.config["BLASTDB_PATH"], dbRoute )
		blastDir = os.path.join( app.config["BLAST_PATH"] )

		blastRoute = os.path.join(blastDir, program)

		print blastRoute, dbRoute

		if (program.startswith('blastn') ):
			cmd = subprocess.Popen( blastRoute.split() + ['-db', dbRoute, '-evalue', evalue,'-num_alignments', '10', '-num_threads', '1', '-outfmt', '15', '-perc_identity' , perc, '-strand', 'plus'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		else:
			cmd = subprocess.Popen( blastRoute.split() + ['-db', dbRoute, '-evalue', evalue,'-num_alignments', '10', '-num_threads', '1', '-outfmt', '15', '-matrix', matrix, '-strand', 'plus'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		
		out,error = cmd.communicate(query)

		print out, error

		if not out:
			flash('BLAST output error')

		results = parseBlastResult(out)

		if not results:
			flash('No hits found')

		return render_template('blast_result.html', title='BLAST result', result = results, idType = idType, maxLen = max([ row["len"] for row in results ]) )	

	return render_template('blast.html', title='BLAST to MarpoDB')

@app.route('/hmmer', methods=['GET', 'POST'])
def hmmer():
	title = 'HMMER search'

	if request.method == 'POST':
			smaFile = request.files['file']
			dataset = request.form['dataset']

			if smaFile.filename:
				smaContent = smaFile.read()

				hexer = md5.new()
				hexer.update(smaContent)

				smaFileName = os.path.join(app.config["TEMP_PATH"], hexer.hexdigest() )

				smaFile = open(smaFileName, 'w')
				smaFile.write(smaContent)
				smaFile.close()

				hmmFileName = smaFileName + '.hmm'

				cmd = subprocess.Popen( ['hmmbuild', '--informat', 'STOCKHOLM', hmmFileName, smaFileName], stderr=subprocess.PIPE,  stdin=subprocess.PIPE, stdout=subprocess.PIPE )
				out, err = cmd.communicate(smaContent)

				if not err:
					tableFileName = hmmFileName + '.tblout'
					
					protDB = os.path.join(app.config["FASTA_PATH"], "Prot_{0}.fa".format(dataset))

					cmd = subprocess.Popen( ['hmmsearch', '--tblout', tableFileName, hmmFileName, protDB], stderr=subprocess.PIPE, stdin=subprocess.PIPE, stdout=subprocess.PIPE )
					out, err = cmd.communicate()
					
					print out
					print err

					if not err:
						try:
							errorMessage = ''
							results = parseHMMResult(tableFileName)
						except Exception as e:
							print e
							errorMessage = 'Failed to run hmmsearch. Please report if happens again.'

						if not errorMessage:
							if results:
								return render_template('hmmer_results.html', title='HMMER result', result = results )
							else:
								errorMessage = 'No hits found' 
					else:
						errorMessage = 'Error running hmmsearch. Check your input file.'

				else:
					errorMessage = 'Error running hmmbuild. Check your input file.'

			else:
				errorMessage = 'No file selected'

			try:
				os.remove(hmmFileName)
				os.remove(smaFileName)
				os.remove(tableFileName)
			except:
				pass

			flash(errorMessage)

	return render_template('hmmer.html', title=title)

@app.route('/export/gene')
def exportGene():
	dbid = request.args.get('dbid','')

	if not dbid:
		return ('', 204)

	marpodbSession = marpodb.Session()
	gene = marpodbSession.query(Gene).filter(Gene.dbid == dbid).first()

	if not gene:
		return ('', 204)

	response = make_response(gene.record.format("gb"))
	response.headers["Content-Type"] = "application/octet-stream"
	response.headers["Content-Disposition"] = "attachement; filename={0}".format(dbid+'.gb')
	return response

@app.route('/recode')
def recode():
	dbid 	= request.args.get('dbid', '')
	seqType = request.args.get('seqType', '')

	if not dbid:
		abort(404)

	session = marpodb.Session()
	gene = session.query(Gene).filter(Gene.dbid == dbid).first()

	if not gene:
		abort(404)

	print seqType

	if seqType == 'cds':
		seq = gene.cds.seq

	elif seqType == 'promoter':
		seq = gene.promoter.seq

	elif seqType == "promoter5":
		seq = gene.promoter.seq
		if gene.utr5:
			seq += gene.utr5.seq
	else:
		abort(404)

	print "Seq: ", gene.promoter.seq

	seq = Seq(seq)

	BsaI_F = "GGTCTC"
	BsaI_F_replace = ["GGaCTC","GGTCcC","GGTaTC"]
	BsaI_R = "GAGACC"
	BsaI_R_replace = ["GAaACC","GAGAtC" ,"GAGgCC"]
	SapI_F = "GCTCTTC"
	SapI_F_replace = ["GCcCTTC","GCTCgTC","GCTgTTC"]
	SapI_R = "GAAGAGC"
	SapI_R_replace = ["GAgGAGC","GAAGgGC","GAAaAGC"]
	BsaI_sitesF = []
	BsaI_sitesR = []
	SapI_sitesF = []
	SapI_sitesR = []
	
	BsaI_sitesF = BsaI_sitesF + recfind(BsaI_F, seq)
	BsaI_sitesR = BsaI_sitesR + recfind(BsaI_R, seq)
	SapI_sitesF = SapI_sitesF + recfind(SapI_F, seq)
	SapI_sitesR = SapI_sitesR + recfind(SapI_R, seq)

	orgseq = seq
	for site in BsaI_sitesF:
		for i in range(3):
			if (site % 3 == i):
				seq = seq[:site] + BsaI_F_replace[i] + seq[site + 6:]
	
	for site in BsaI_sitesR:
		for i in range(3):
			if (site % 3 == i):
				seq = seq[:site] + BsaI_R_replace[i] + seq[site + 6:]
	
	for site in SapI_sitesF:
		for i in range(3):
			if (site % 3 == i):
				seq = seq[:site] + SapI_F_replace[i] + seq[site + 7:]
				
	for site in SapI_sitesR:
		for i in range(3):
			if (site % 3 == i):
				seq = seq[:site] + SapI_R_replace[i] + seq[site + 7:]
					
	BsaI_sites = BsaI_sitesF+BsaI_sitesR
	SapI_sites = SapI_sitesF+SapI_sitesR

	session.close()

	return render_template('recode.html', dbid = dbid, seq = orgseq, seqType = seqType, BsaI = BsaI_sites, SapI=SapI_sites, newseq=seq, title='Recode sequence')

@app.route('/export/recode')
def exportrecode():
	geneName = request.args.get('geneName','')
	seq = request.args.get('seq','')
	seqType = request.args.get('type','')
	
	if not geneName:
		return ('', 204)

	record = recodeExport( geneName, seq, seqType )

	response = make_response(record.format("gb"))
	response.headers["Content-Type"] = "application/octet-stream"
	response.headers["Content-Disposition"] = "attachement; filename={0}".format(geneName+'.gb')
	return response

@app.route('/help')
def help():
	return render_template('help.html', title = "Help")

@app.route('/about')
def about():
	return render_template('about.html', title = "About")

@app.route('/whatsnew')
def whatsnew():
	return render_template('whatsnew.html', title = "What's new?")

