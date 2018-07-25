from flask import Flask
from flask_restless import APIManager
from sqlalchemy.orm import scoped_session

from ..core import *

app = Flask(__name__)
app.secret_key = 'HJKDGSA&^D%HJKN.zczxcoasdk2194uru'
app.debug = True

# manager = APIManager(app, session = scoped_session(marpodb.Session) )
# manager.create_api(Gene, methods = ['GET'], include_columns = ["dbid", "id", "seq"])

import routes
