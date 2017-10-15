import json
import subprocess
import os

def benchop(problems=[1,2,3,4],methods=['MC','MC-S','QMC-S','MLMC','MLMC-A','FFT','FGL','COS','FD','FD-NU','FD-AD','RBF','RBF-FD','RBF-PUM','RBF-LSML','RBF-AD','RBF-MLT']):
	for problem in problems:
		for method in methods:
			subprocess.call("docker service create --replicas 1 --name benchop"+ problem + method +" maxwatson142/worker python benchop "+problem+" "+method)
	return 0

def getResult():
	try:
		with open('/home/ubuntu/result.txt', 'r') as file:
			result = json.load(file)
			file.close()
			return result
	except:
		return {}

def postResult(result):
	with open('/home/ubuntu/result.txt', mode='a+') as file:
		feeds = json.load(file)
		feeds.append(result)
		json.dump(feeds, file)
		file.close()
		return 0

def deleteResult():
	with open('home/ubuntu/result.txt', mode='w') as f:
		json.dump([], f)
		f.close()