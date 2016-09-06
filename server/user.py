from wtforms import Form, StringField, PasswordField, validators


class RegisterForm(Form):
	username = StringField("Username",[validators.Length(min=2, max=20)])
	email = StringField("E-mail", [validators.Length(min=6, max=120), validators.Email()])
	password = PasswordField("Password", [validators.Length(min=8, max=20)])
	affiliation = StringField("Affiliation", [validators.Length(min=1)])

class LoginForm(Form):
	username =  StringField("Username")
	password = PasswordField("Password")
