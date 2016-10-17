# MarpoDB - Server 
#
# MIT License
#
# Copyright (c) 2016 HaseloffLab - Bernardo Pollak, Mihails Delmans & Jim Haseloff
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import itertools
import collections

from flask import Flask, session, redirect, url_for, escape, request, render_template, make_response, flash, abort
import subprocess

from marpodb import processQuery, getGeneCoordinates, getCDSDetails, exportGB, generateNewMap, getGeneHomolog, parseBlastResult, getTopGenes, getUserData

import psycopg2

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC

from flask_sqlalchemy import SQLAlchemy, sqlalchemy
from flask_user import UserMixin, SQLAlchemyAdapter, UserManager, LoginManager
from flask_login import login_user, login_required, logout_user, current_user, user_logged_in
from user import RegisterForm, LoginForm

import os.path

app = Flask(__name__)
app.secret_key = 'xxx'
app.config['SQLALCHEMY_DATABASE_URI'] = "postgresql:///userdb"
app.debug = True
userDB = SQLAlchemy(app)

loginManager = LoginManager()
loginManager.init_app(app)

class User(userDB.Model, UserMixin):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	username = userDB.Column(userDB.String(50), nullable = False, unique = True)
	password = userDB.Column(userDB.String(50), nullable = False, server_default = '')
	email = userDB.Column(userDB.String(120), nullable = False, unique=True)
	affiliation = userDB.Column(userDB.Text())
	active = userDB.Column('is_active', userDB.Boolean(), nullable=False, server_default='0')

class StarGene(userDB.Model):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	userid = userDB.Column(userDB.Integer, userDB.ForeignKey('user.id'))
	geneid = userDB.Column(userDB.Integer)
	def __init__(self, userid = None, geneid = None):
		self.userid = userid
		self.geneid = geneid

db_adapter = SQLAlchemyAdapter(userDB,  User)
user_manager = UserManager(db_adapter, app)

@app.context_processor
def user_data():
	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()
	
	userData = getUserData(cur, StarGene, current_user, session)
	
	print "userData:", userData

	cur.close()
	conn.close()

	return dict(user_data = userData)

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html', title='404'), 404

@app.errorhandler(500)
def internal_error(e):
	return render_template('500.html', title='500'), 500

@loginManager.user_loader
def load_user(username):
	User.query.filter(User.username == username).first()

@app.route('/geneName')
def geneName():
	geneid = request.args.get('geneid','')
	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()

	geneName = getGeneHomolog(cur, geneid)

	return geneName

@app.route('/')
def index():
	return redirect(url_for('query'))

@app.route('/about')
def about():
	return render_template("about.html", title='About')

@app.route('/help')
def help():
	return render_template("help.html", title='Help')

@app.route('/query', methods=['GET', 'POST'])
def query():
	loginForm = LoginForm()
	if request.method == 'POST':
		scope = request.form.getlist("searchScope")
		scope = '|'.join(scope)
		term = request.form["term"]

		if scope !="" and term !="":
			return redirect(url_for('result', term=request.form['term'], scope = scope))
		else:
			flash('Empty search term or scope')
	return render_template("query.html", loginForm = loginForm)

@app.route('/results')
def result():
	nHits = 5

	columns = {'pfam_hit' : {'name', 'acc', 'e_val', 'descr'}, 'blastp_hit' : ['protein_name', 'gene_name', 'origin', 'e_val'], 'blastm_hit' : ['protein_name', 'gene_name', 'origin', 'e_val'], 'cds': ['name'] }

	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()

	term = request.args.get('term','')
	scope = request.args.get('scope','').split('|')

	if not (term and scope):
		 flash('Empyt search term or scope')
		 return redirect(url_for('query'))

	table = processQuery(cur, scope, term, columns, nHits)

	conn.commit()
	cur.close()
	conn.close()

	if len(table['data']) != 0:
		return render_template('results.html', table = table, title='Query results')

	else:
		flash('No entries matching the query were found.')
		return redirect(url_for('query'))

@app.route('/details')
def details():
	name = request.args.get('name','')

	geneName = name.split('.')[0]
	
	if len( name.split('|') ) == 2:
		cdsName = name
	else:
		cdsName = ""

	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()

	cur.execute("SELECT id, seq, alias FROM gene WHERE name=%s", (geneName, ))

	gene = cur.fetchone()
	geneid = gene[0]
	geneSeq = gene[1]
	geneAlias = gene[2]
		
	if not geneSeq:
		abort(404)

	geneCoordinates = getGeneCoordinates(cur, geneName)

	if not (geneCoordinates['mrnas'] and geneCoordinates['cdss'] and geneCoordinates['gene']):
		abort(500)

	if not cdsName:
		cdsName = geneCoordinates['cdss'].keys()[0]

	cdsDetails = getCDSDetails(cur, cdsName)

	stared = False
	if current_user.is_authenticated:
		if StarGene.query.filter(StarGene.userid == current_user.id, StarGene.geneid == geneid).first():
			stared = True
	else:
		
		if not "stars" in session:
			session["stars"] = ""
		print "DETAILS STARS: ", session["stars"]
	 	print "DETAILS ID: ", str(geneid), session["stars"].find(str(geneid))
		if session["stars"].find(str(geneid)) > -1:
			stared = True
		print "STARRED: ", stared
	return render_template('details.html', cdsName = cdsName, geneName = geneName, geneSeq = geneSeq, gene = geneCoordinates['gene'], mrnas = geneCoordinates['mrnas'], cdss = geneCoordinates['cdss'], blastm = cdsDetails['blastm'], blastp = cdsDetails['blastp'], stared = stared, alias = geneAlias)

@app.route('/user_<cmd>')
def addToFavourite(cmd=None):
	print cmd
	if cmd == 'stargene':
		geneName = request.args.get('genename','')

		conn = psycopg2.connect("dbname=marpodb")
		cur = conn.cursor()

		cur.execute("SELECT id FROM gene WHERE name=%s", (geneName, ))
		geneid = cur.fetchone()[0]

		if current_user.is_authenticated:

			starGene = StarGene.query.filter(StarGene.userid == current_user.id, StarGene.geneid == geneid).first()

			if geneid:
				if starGene:
					userDB.session.delete(starGene)
					userDB.session.commit()
				else:
					newStar = StarGene(current_user.id, geneid)
					userDB.session.add(newStar)
					userDB.session.commit()
				return "OK", 200, {'Content-Type': 'text/plain'}
			else:
				abort(500)
		else:

			if geneid:
				if not "stars" in session:
					session["stars"] = ""
                                
				if str(geneid) in session["stars"]:
                                        session["stars"] = session["stars"].replace(":{0}:".format(geneid), "")
                                else:
                                        session["stars"] += ":{0}:".format(geneid)
				print "STARS: ", session["stars"]
				return "OK", 200, {'Content-Type': 'text/plain'}
                        else:
                                abort(500)
			
# Helper function for recursive search of restriction sites
			
def recfind(pattern, string, where_should_I_start=0):
    # Save the result in a variable to avoid doing the same thing twice
    pos = string.find(pattern, where_should_I_start)
    if pos == -1:
        # Not found!
        return []
    # No need for an else statement
    return [pos] + recfind(pattern, string, pos + len(pattern))
			
@app.route('/recode', methods=['GET'])
def recode():
	seq = request.args.get('seq','')
	strand = request.args.get('strand','')
	seq_type = request.args.get('type','')
	cdsName = request.args.get('cdsName','')

	seq = Seq(seq)

	if strand == '-':
		seq = seq.reverse_complement()
	seq = str(seq)
		
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
#	for search in SapI:
#		SapI_sites = SapI_sites + recfind(search, seq)
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
	
		
	return render_template('recode.html', seq=orgseq, strand=strand, seq_type=seq_type, cdsName=cdsName, BsaI=BsaI_sites, SapI=SapI_sites, newseq=seq, title='Recode sequence')

@app.route('/blast', methods=['GET', 'POST'])
def blast():
	
	if request.method == 'POST':
		query = request.form["query"]
		evalue = request.form["evalue"]
		program = request.form["program"]
		matrix = request.form["matrix"]
		perc = request.form["perc"]
		
		if query !="":
			return redirect(url_for('blast_result', query=query, evalue=evalue, program=program,matrix=matrix, perc=perc))
		else:
			flash('Empty query')
	
	return render_template('blast.html', title='BLAST to MarpoDB')

@app.route('/blast_result')
def blast_result():
	query = request.args.get('query','')
	evalue = request.args.get('evalue','')
	program = request.args.get('program','')
	matrix = request.args.get('matrix','')
	perc = request.args.get('perc','')

	# Validate fasta PENDING!
	
	if (program.startswith('blastn') or program == 'tblastx'):
		route = 'blast/MarpoDB_Genes'
	else:
		route = 'blast/MarpoDB_Proteins'

	if (program.startswith('blastn') ):
		cmd = subprocess.Popen( program.split() + ['-db', route, '-evalue', evalue,'-num_alignments', '10', '-num_threads', '1', '-outfmt', '13', '-perc_identity' , perc], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	else:	
		cmd = subprocess.Popen( program.split() + ['-db', route, '-evalue', evalue,'-num_alignments', '10', '-num_threads', '1', '-outfmt', '13', '-matrix', matrix], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	
	out,error = cmd.communicate(query)

	print "CMD output: ", out
	print "CMD error:", error
	
	return render_template('blast_result.html', title='BLAST result', result = parseBlastResult(out))

@app.route('/login', methods=['POST'])
def login():

	session.pop("stars", None)

	form = LoginForm(request.form)
	user = User.query.filter(User.username == form.username.data).first()
	if user and user_manager.verify_password(form.password.data, user):
		login_user(user)
		#flash("Logged in as " + user.username)
	else:
		flash("Invalid username or/and password")
	return redirect(url_for('index'))

@app.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for('index'))

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

			newUser = User(username=form.username.data, email=form.email.data, password=user_manager.hash_password(form.password.data), affiliation=form.affiliation.data, active = True)

			userDB.session.add(newUser)
			userDB.session.commit()
		
			try:
				generateNewMap(User)
			except:
				pass

			flash('You have been successfully registred!')
			return redirect(url_for('map'))

	return render_template('registration.html', title = 'Registration', form = form)

@app.route('/map')
def map():
	if not os.path.isfile('server/static/img/map.png'):
		generateNewMap(User)
		print "Map not found"
	return render_template('map.html', title='Community map')

@app.route('/top')
def top():

	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()

	geneTop = getTopGenes(cur, StarGene, 10)

	cur.close()
	conn.close()

	return render_template('top.html', list = geneTop, title='Top genes')

@app.route('/export/cds')
def exportCds():
	name = request.args.get('cds','')

	if not name:
		return ('', 204)

	conn = psycopg2.connect("dbname=marpodb")
	cur = conn.cursor()

	record = exportGB(cur, name)

	if not record:
		return ('', 204)

	response = make_response(record.format("gb"))
	response.headers["Content-Type"] = "application/octet-stream"
	response.headers["Content-Disposition"] = "attachement; filename={0}".format(name+'.gb')
	return response

@app.route('/export/recode')
def exportrecode():
	name = request.args.get('cdsName','')
	seq = request.args.get('seq','')
	seq_type = request.args.get('type','')
	
	if not name:
		return ('', 204)

	record = SeqRecord( Seq( seq, IUPAC.unambiguous_dna ), name = name, description = "Generated by MarpoDB", id="M. polymorpha")
	
	strand = 1;
	location = FeatureLocation(11, len(seq)-11, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type=seq_type, id="Part" ) )
	
	location = FeatureLocation(0, 6, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="BsaI_F", id="BsaI") )
	
	location = FeatureLocation(7, 11, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="misc_recomb", id="Front" ) )

	strand = -1;
	location = FeatureLocation(len(seq)-6, len(seq), strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="BsaI_R", id="BsaI" ) )
	

	location = FeatureLocation(len(seq)-11, len(seq)-7, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="misc_recomb", id="Front" ) )


	response = make_response(record.format("gb"))
	response.headers["Content-Type"] = "application/octet-stream"
	response.headers["Content-Disposition"] = "attachement; filename={0}".format(name+'.gb')
	return response

if __name__ == '__main__':
	app.run(debug=True,host='0.0.0.0', port = 8081)

