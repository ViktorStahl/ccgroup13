import json
import subprocess
import os
import thread
#['MC','MC-S','QMC-S','MLMC','MLMC-A','FFT','FGL','COS','FD','FD-NU','FD-AD','RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']
def benchop(problems=[1,2,3,4],methods=['COS','FD', 'RBF-FD']):
	thread.start_new_thread(startServices, (problems, methods))
	return '0'

def startServices(problems, methods):
	for problem in problems:
		for method in methods:
			subprocess.call("docker service create --replicas 1 --name benchop"+ str(problem) + method +" maxwatson142/worker python benchop "+str(problem)+" "+method, shell=True)
	return '0'

def getResult():
	try:
		with open('result.txt', 'r') as file:
			result = json.load(file)
			file.close()
			return result
	except:
		return {}

def postResult(result):
	feeds = {}
	try:
		with open('result.txt', mode='r') as file:
			feeds={}
			try:
				feeds = json.load(file)
			except:
				pass
	except:
		pass
	with open('result.txt', mode='w+') as file:
		feeds.update(result)
		json.dump(feeds, file)
		file.close()
		return '0'

def deleteResult():
	with open('result.txt', mode='a+') as f:
		json.dump({}, f)
		f.close()
		return '0'