#Import packages and modules from website folder
import sqlite3
import pandas as pd
from website import create_app
from website.db_setup import create_the_database
from os import path

#If database does not exist, create database 
db_name = 'data.db'
if not path.exists('website/' + db_name):
        create_the_database()
        print('Database created.')


app = create_app()

if __name__ == '__main__':
    app.run(debug=True)
