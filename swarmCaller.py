import json
import subprocess
import os

def benchop(problems=[],methods=[]):
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

def putResult(result):
	with open('/home/ubuntu/result.txt', 'w') as file:
		json.dump(result, file)
		file.close()
		return 0

def deleteResult():
	try:
		os.remove('/home/ubuntu/result.txt')
	except:
		pass