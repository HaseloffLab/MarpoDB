from ..core import *
from ..server import app

from .user import *
from .backend import *

from flask import redirect, render_template, url_for, request

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

@app.route('/results', methods = ['GET'])
def result():
	nHits = 5
	displayColumns = ['origin', 'description', 'eVal']

	session = marpodb.Session()

	term = request.args.get('term','')
	scope = request.args.get('scope','').split('|')

	if not (term and scope):
		 flash('Search term and scope should not be empty!')
		 return redirect(url_for('query'))

	table = processQuery(session, scope, term, displayColumns, nHits)

	session.close()

	if len(table['data']) != 0:
		return render_template('results.html', table = table, title='Query results')
	else:
		flash("No entries found")
		return redirect(url_for('index'))