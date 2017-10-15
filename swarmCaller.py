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
		with open('/home/result.txt', 'r') as file:
			result = json.load(file)
			file.close()
			return result
	except:
		return {}

def postResult(result):
	with open('/home/result.txt', mode='a+', encoding='utf-8') as file:
		feeds = json.load(file)
		feeds.append(result)
		json.dump(feeds, file)
		file.close()
		return 0

def deleteResult():
	with open('home/result.txt', mode='w', encoding='utf-8') as f:
		json.dump([], f)
		f.close()