import requests
def post(apiIP,problem, method, res, timeTaken):
	r = requests.post("http://"+apiIP+":5000/api/v1/result", data={str(problem)+method:{'problem': problem, 'method': method, 'time': timeTaken, 'result': res}})

if __name__ == '__main__':
	post('130.239.81.130', 1, 'COS', 0.923, 2.41)