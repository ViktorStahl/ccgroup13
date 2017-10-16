from  flask import Flask, jsonify
import sys
import subprocess
import oct2py
import requests
from oct2py import octave
app = Flask(__name__)
app.route('/benchops/api/v1.0/results', methods=['GET'])

#def benchop(problem, method, apiIP):
def benchop(problem,method):
    oc = oct2py.Oct2Py()
    oc.addpath('./')
    oc = oct2py.Oct2Py()
    x = oc.feval("benchop",problem,method);
    print(x);
    return x
    #timeTaken = time.time()-start
    #res = subprocess.call(["Ocatave -r benchop "+problem+" "+method],shell=True) <- IF oct2py doesn't work
    #r = requests.post("http://"+apiIP+":5000/api/v1/result", data={str(problem)+method:{'problem': problem, 'method': method, 'time': timeTaken, 'result': res})
	#if (r.status_code//100)%10==2:
	#	return
	#else:
		#ErrorHandling
	#	print("Something went wrong") 
	#	return
if __name__ == "__main__":
    benchop(sys.argv[1],sys.argv[2])
