#!flask/bin/python
from flask import Flask, jsonify, request, send_from_directory, render_template, redirect, url_for
import sys
import swarmCaller
import time
import json
import logging
from logging.handlers import RotatingFileHandler

app = Flask(__name__, static_folder='static')

@app.route('/api/v1/benchop', methods=['GET'])
def benchop():
	result = swarmCaller.deleteResult()
	result = swarmCaller.benchop()
	return 'HTTP status 200 (OK) Please wait while the data is processed \n The data can be found by sending a GET request to /api/v1/result in about 5 minutes'

@app.route('/', methods=['GET'])
def sayHi():
	htmlstr = "<!DOCTYPE html> \n \
		<html> \n \
		<head> \n \
		<title>BENCHOP</title> \n \
		</head> \n \
		<body> \n \
		<h1>BENCHOP API</h1> \n \
		<p>Welcome to the BENCHOP API</p> \n \
		</body> \n \
		</html>"
	return htmlstr

@app.route('/api/v1/benchop', methods=['POST'])
def benchop_post():
	result = swarmCaller.deleteResult()
	problems = request.form['problem'].split()
	methods = request.form['method'].split()
	result = swarmCaller.benchop(problems, methods)
	return 'HTTP OK 200 Please wait while the data is processed \n The data can be found by sending a GET request to /api/v1/result in about 5 minutes'
	
@app.route('/api/v1/result', methods=['GET'])
def result_get():
	result = swarmCaller.getResult()
	return json.dumps(result)

@app.route('/api/v1/result', methods=['POST'])
def result_post():
	app.logger.info(request.form)
	resdict = {}
	resdict[request.form['problem'].encode('UTF8')+request.form['method'].encode('UTF8')] = {}
	for word in request.form:
		resdict[request.form['problem'].encode('UTF8')+request.form['method'].encode('UTF8')][word.encode('UTF8')] = request.form[word.encode('UTF8')].encode('UTF8')
	app.logger.info(resdict)
	result = swarmCaller.postResult(resdict)
	return '0'

@app.route('/api/v1/result', methods=['DELETE'])
def apiDeleteResult():
	result = swarmCaller.deleteResult()
	return 'HTTP status 200 (OK)'

if __name__ == '__main__':
	handler = RotatingFileHandler('foo.log', maxBytes=10000, backupCount=1)
	handler.setLevel(logging.INFO)
	app.logger.addHandler(handler)
	app.run(host='0.0.0.0',debug=True)