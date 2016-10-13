from partsdb.partsdb import PartsDB
from tables import *

from flask import Flask, session, redirect, url_for, escape, request, render_template, make_response, flash, abort
from flask.ext.sqlalchemy import SQLAlchemy
from flask_user import UserMixin, SQLAlchemyAdapter, UserManager, LoginManager
from flask_login import login_user, login_required, logout_user, current_user, user_logged_in

from user import RegisterForm, LoginForm
from system import getUserData, generateNewMap, getTopGenes, processQuery

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
	app.run(debug=True,host='0.0.0.0', port = 8082)








