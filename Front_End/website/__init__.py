#Import statements
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from os import path
import sqlite3
import pandas as pd

db = SQLAlchemy()
db_name = "data.db"

#Define application factory
def create_app():
    app = Flask(__name__)
    app.config['SECRET_KEY'] = 'vegan'
    app.config ['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + db_name
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True

    #initialise database with flask app
    db.init_app(app)

    #register blueprints with application factory
    from . import home
    app.register_blueprint(home.bp)

    return app
