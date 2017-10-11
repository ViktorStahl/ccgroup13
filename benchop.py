import sys
import subprocess
import oct2py
import requests

def benchop(problem, method, apiIP):
	oc = oct2py.Oct2Py()
	res = oc.benchop(problem,method)
	#res = subprocess.call(["Ocatave -r benchop "+problem+" "+method],shell=True) <- IF oct2py doesn't work
	r = requests.post(apiIP, data={'problem': problem, 'method': method, 'time': res.time, 'result': toJSONstr(res)})
	if (r.status_code//100)%10==2:
		return
	else:
		#ErrorHandling
		print("Something went wrong")
		return
if __name__ == '__main__':
	assert(len(sys.argv)>=4)
	benchop(sys.argv[1],sys.argv[2], sys.argv[3])