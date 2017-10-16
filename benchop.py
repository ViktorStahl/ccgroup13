import sys
import subprocess
import oct2py
import requests
import time

#def benchop(problem, method, apiIP):
def benchop(problem,method):
		
	oc.addpath('/home/ccgroup13/')
	start = time.time()
	res = oc.benchop(problem, method)
	timeTaken = time.time()-start
	#res = subprocess.call(["Ocatave -r benchop "+problem+" "+method],shell=True) <- IF oct2py doesn't work
	r = requests.post("http://"+apiIP+":5000/api/v1/result", data={str(problem)+method:{'problem': problem, 'method': method, 'time': timeTaken, 'result': res}})
	if (r.status_code//100)%10==2:
		return
	else:
		#ErrorHandling
		print("Something went wrong") 
		return
if __name__ == '__main__':
	benchop(sys.argv[1],sys.argv[2])
	#assert(len(sys.argv)>=2)
	#benchop(sys.argv[1],sys.argv[2])
	#benchop(sys.argv[1],sys.argv[2], sys.argv[3])
