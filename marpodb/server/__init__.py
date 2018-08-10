from flask import Flask, session
from flask_login import LoginManager, AnonymousUserMixin, UserMixin, current_user
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt

from sqlalchemy.orm import scoped_session
from sqlalchemy.ext.hybrid import hybrid_property

from wtforms import Form, StringField, PasswordField, validators

from ..core import *

app = Flask(__name__)
app.secret_key = 'HJKDGSA&^D%HJKN.zczxcoasdk2194uru'
app.debug = True

############################## LOGIN ##############################

class RegisterForm(Form):
	username = StringField("Username",[validators.Length(min=2, max=20)])
	email = StringField("E-mail", [validators.Length(min=6, max=120), validators.Email()])
	password = PasswordField("Password", [validators.Length(min=8, max=20)])
	affiliation = StringField("Affiliation", [validators.Length(min=1)])

class LoginForm(Form):
	username =  StringField("Username")
	password = PasswordField("Password")

class AnonymousUser(AnonymousUserMixin):
	@property
	def data(self):
		starGenes = {}
		if "starGenes" in session:
			starGenes = session["starGenes"]

		return {"starGenes" : starGenes}

loginManager = LoginManager()

loginManager.anonymous_user = AnonymousUser

loginManager.init_app(app)

############################## USERDB ##############################

app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = "postgresql:///userdb3"

userDB = SQLAlchemy(app)

class User(userDB.Model, UserMixin):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	username = userDB.Column(userDB.String(50), nullable = False, unique = True)
	password = userDB.Column(userDB.String(500), nullable = False, server_default = '')
	email = userDB.Column(userDB.String(120), nullable = False, unique=True)
	affiliation = userDB.Column(userDB.Text())
	active = userDB.Column('is_active', userDB.Boolean(), nullable=False, server_default='0')
	stars = relationship("StarGene", backref = "user")

	@hybrid_property
	def data(self):
		return {"starGenes" : { star.cdsdbid: star.genename for star in self.stars } }

@loginManager.user_loader
def load_user(userID):
	user =  User.query.filter(User.id== userID).first()
	return user

class StarGene(userDB.Model):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	userid = userDB.Column(userDB.Integer, userDB.ForeignKey('user.id'))
	cdsdbid = userDB.Column(userDB.Text)
	genename = userDB.Column(userDB.Text)

	def __init__(self, userid = None, cdsdbid = None, genename = None):
		self.userid = userid
		self.cdsdbid = cdsdbid
		self.genename = genename

userDB.create_all()

bcrypt = Bcrypt(app)

import routes
