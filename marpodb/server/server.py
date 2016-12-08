from partsdb.partsdb import PartsDB
from tables import *

from flask import Flask, session, redirect, url_for, escape, request, render_template, make_response, flash, abort
from flask.ext.sqlalchemy import SQLAlchemy
from flask_user import UserMixin, SQLAlchemyAdapter, UserManager, LoginManager
from flask_login import login_user, login_required, logout_user, current_user, user_logged_in

from user import RegisterForm, LoginForm
from system import getUserData, generateNewMap, getTopGenes, processQuery, getGeneCoordinates

import os

app = Flask(__name__)
app.secret_key = 'HJKDGSA&^D%HJKN.zczxcoasdk2194uru'
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

@loginManager.user_loader
def load_user(username):
	User.query.filter(User.username == username).first()

@app.context_processor
def user_data():
	marpodbSession = marpodb.Session()
	userData = getUserData(StarGene, current_user, session, marpodbSession)
	
	marpodbSession.close()
	print "userData:", userData

	return dict(user_data = userData)

marpodb = PartsDB('postgresql:///testdb', Base = Base)

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

		if scope !="" and term !="":
			return redirect(url_for('result', term=request.form['term'], scope = scope))
		else:
			flash('Empty search term or scope')
	return render_template("query.html", loginForm = loginForm)


@app.route('/login', methods=['POST'])
def login():

	session.pop("stars", None)

	form = LoginForm(request.form)
	user = User.query.filter(User.username == form.username.data).first()
	if user and user_manager.verify_password(form.password.data, user):
		login_user(user)
	else:
		flash("Invalid username or/and password")
	return redirect(url_for('index'))

@app.route('/results')
def result():
	nHits = 5
   	columns = { 'pfamhit'	: [PfamHit.name, PfamHit.acc, PfamHit.eVal, PfamHit.description],\
   				'blastphit': [BlastpHit.proteinName, BlastpHit.geneName, BlastpHit.origin, BlastpHit.eVal],\
   				'cds'		: [CDS.dbid]\
   	}

	marpodbSession = marpodb.Session()

	term = request.args.get('term','')
	scope = request.args.get('scope','').split('|')

	if not (term and scope):
		 flash('Empyt search term or scope')
		 return redirect(url_for('query'))

	table = processQuery(marpodbSession, scope, term, columns, nHits)

	marpodbSession.close()

	if len(table['data']) != 0:
		return render_template('results.html', table = table, title='Query results')
	else:
		return "ERROR: no data"
@app.route('/details')
def details():
	dbid = request.args.get('dbid','')

	marpodbSession = marpodb.Session()

	locus = None
	gene  = None
	cds   = None

	locus = marpodbSession.query(Locus).filter(Locus.dbid == dbid).first()

	if not locus:
		
		gene = marpodbSession.query(Gene).filter(Gene.dbid == dbid).first()
		
		if not gene:
			cds = marpodbSession.query(CDS).filter(CDS.dbid == dbid).first()
			if cds:
				gene  = marpodbSession.query(Gene).filter(Gene.cdsID == cds.id).first()
				locus = marpodbSession.query(Locus).filter(Locus.id == Gene.locusID).\
						filter(Gene.cdsID == cds.id).first()
				cdsdbid = cds.dbid

		else:
			locus = marpodbSession.query(Locus).filter(Locus.id == gene.locusID).first()
			cdsdbid = marpodbSession.query(CDS.dbid).filter(CDS.id == gene.cdsID).first()[0]
	else:
		cdsdbid = marpodbSession.query(CDS.dbid).filter(CDS.id == Gene.cdsID).\
					filter(Gene.locusID == locus.id).first()[0]
	response = getGeneCoordinates(marpodbSession, locus.id)

	print response['genes']

	return render_template('details2.html', cdsDBID = cdsdbid, geneCoordinates = response['genes'], seq = response['seq'],  title = "Details for {0}".format(cdsdbid))

	# if not (geneCoordinates['mrnas'] and geneCoordinates['cdss'] and geneCoordinates['gene']):
	# 	abort(500)

	# if not cdsName:
	# 	cdsName = geneCoordinates['cdss'].keys()[0]

	# cdsDetails = getCDSDetails(cur, cdsName)

	# stared = False
	# if current_user.is_authenticated:
	# 	if StarGene.query.filter(StarGene.userid == current_user.id, StarGene.geneid == geneid).first():
	# 		stared = True
	# else:
		
	# 	if not "stars" in session:
	# 		session["stars"] = ""
	# 	print "DETAILS STARS: ", session["stars"]
	#  	print "DETAILS ID: ", str(geneid), session["stars"].find(str(geneid))
	# 	if session["stars"].find(str(geneid)) > -1:
	# 		stared = True
	# 	print "STARRED: ", stared
	# return render_template('details.html', cdsName = cdsName, geneName = geneName, geneSeq = geneSeq, gene = geneCoordinates['gene'], mrnas = geneCoordinates['mrnas'], cdss = geneCoordinates['cdss'], blastm = cdsDetails['blastm'], blastp = cdsDetails['blastp'], stared = stared, alias = geneAlias)

@app.route('/map')
def map():
	if not os.path.isfile('static/img/map.png'):
		generateNewMap(User)
		print "Map not found"
	return render_template('map.html', title='Community map')

@app.route('/top')
def top():

	marpodbSession = marpodb.Session()
	geneTop = getTopGenes(marpodbSession, StarGene, 10)
	marpodbSession.close()

	return render_template('top.html', list = geneTop, title='Top genes')

@app.route('/about')
def about():
	return render_template("about.html", title='About')

@app.route('/help')
def help():
	return render_template("help.html", title='Help')

if __name__ == '__main__':
	app.run(debug=True,host='0.0.0.0', port = 8081)








