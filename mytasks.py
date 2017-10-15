from  flask import Flask, jsonify
import sys
import subprocess
import oct2py
import requests
from oct2py import octave
app = Flask(__name__)
app.route('/benchops/api/v1.0/results', methods=['GET'])

result = []
#def benchop(problem, method, apiIP):
def benchop(problem,method):
    octave.addpath('/home/sujata/Documents/school/YearII/AppliedCloudComputing/myproject/ccgroup13/');
    res = octave.benchop(problem,method, nout=1);
    oc = oct2py.Oct2Py()
    x = oc.feval("benchop",problem,method);
    print(x);
    #result.append(x)

#print (result)
        
	#res = subprocess.call(["Ocatave -r benchop "+problem+" "+method],shell=True) <- IF oct2py doesn't work
	#r = requests.post(apiIP, data={'problem': problem, 'method': method, 'time': res.time, 'result': toJSONstr(res)})
	#if (r.status_code//100)%10==2:
	#	return
	#else:
		#ErrorHandling
	#	print("Something went wrong")
	#	return
if __name__ == '__main__':
        benchop(sys.argv[1],sys.argv[2])
	#assert(len(sys.argv)>=2)
        #benchop(sys.argv[1],sys.argv[2])
	#benchop(sys.argv[1],sys.argv[2], sys.argv[3])
        #print (result)
