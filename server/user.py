# MarpoDB - User Tables Definition
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

from wtforms import Form, StringField, PasswordField, validators


class RegisterForm(Form):
	username = StringField("Username",[validators.Length(min=2, max=20)])
	email = StringField("E-mail", [validators.Length(min=6, max=120), validators.Email()])
	password = PasswordField("Password", [validators.Length(min=8, max=20)])
	affiliation = StringField("Affiliation", [validators.Length(min=1)])

class LoginForm(Form):
	username =  StringField("Username")
	password = PasswordField("Password")
