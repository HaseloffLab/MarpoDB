from ..server import app

from flask_user import UserMixin, UserManager, LoginManager, SQLAlchemyAdapter, current_user
from flask_login import AnonymousUserMixin

from flask_sqlalchemy import SQLAlchemy

from sqlalchemy.ext.hybrid import hybrid_property

from wtforms import Form, StringField, PasswordField, validators

app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = "postgresql:///userdb4"

# UserDB
userDB = SQLAlchemy(app)

class User(userDB.Model, UserMixin):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	username = userDB.Column(userDB.String(50), nullable = False, unique = True)
	password = userDB.Column(userDB.String(500), nullable = False, server_default = '')
	email = userDB.Column(userDB.String(120), nullable = False, unique=True)
	affiliation = userDB.Column(userDB.Text())
	active = userDB.Column('is_active', userDB.Boolean(), nullable=False, server_default='0')
	
	@hybrid_property
	def data(self):
		return {"starGenes" : {}}

class StarGene(userDB.Model):
	id = userDB.Column(userDB.Integer(), primary_key = True)
	userid = userDB.Column(userDB.Integer, userDB.ForeignKey('user.id'))
	cdsdbid = userDB.Column(userDB.Text)
	genename = userDB.Column(userDB.Text)
	def __init__(self, userid = None, cdsdbid = None, genename = None):
		self.userid = userid
		self.cdsdbid = cdsdbid
		self.genename = genename

class AnonymousUser(AnonymousUserMixin):
	def __init__(self):
		self.data = {"starGenes" : {}}

db_adapter = SQLAlchemyAdapter(userDB,  User)
user_manager = UserManager(db_adapter, app)
userDB.create_all()

# Login Manager
loginManager = LoginManager()
loginManager.init_app(app)
loginManager.anonymous_user = AnonymousUser

@loginManager.user_loader
def load_user(username):
	User.query.filter(User.username == username).first()

@app.context_processor
def user_data():
	return {"user_data" : current_user.data}

# Forms
class RegisterForm(Form):
	username = StringField("Username",[validators.Length(min=2, max=20)])
	email = StringField("E-mail", [validators.Length(min=6, max=120), validators.Email()])
	password = PasswordField("Password", [validators.Length(min=8, max=20)])
	affiliation = StringField("Affiliation", [validators.Length(min=1)])

class LoginForm(Form):
	username =  StringField("Username")
	password = PasswordField("Password")